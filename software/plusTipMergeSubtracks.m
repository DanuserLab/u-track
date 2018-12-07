function [dataMatMerge,dataMatReclass,dataMatCrpSecMic,projData]=plusTipMergeSubtracks(projData,dataMat,varargin)
% plusTipMergeSubtracks merges growth fgaps with the flanking growth phases
%
%
% 
% SYNOPSIS  : [dataMatMerge,dataMatReclass,dataMatCrpSecMic, projData]=...
%                   plusTipMergeSubtracks(projData,dataMat)
%
% INPUT
% projData  : structure containing frame rate, pix size info, etc.
%             if dataMatrix isn't given, then projData should have the
%             field nTrack_sF_eF_vMicPerMin_trackType_lifetime_totalDispPix
%             which contains reclassified fgaps and bgaps.
% dataMat   : mainly to be used in plusTipPostTracking, where it is the
%             data matrix prior to reclassification.
%
% fgapReclassScheme : a scalar specifyinh the scheme for reclassifying
%                      fgaps into undetected growth
%                      1 - using growth velocity of the last 2-3 frames
%                      before
%                      2 - using full growth subtrack velocity
%                      3 - unimodal reclassification of fgaps per project
%                      4 - skip reclassification
%
% bgapReclassScheme : a scalar specifying the scheme for reclassifying
%                      fgaps into undetected growth
%                      1 - if speed is lower than the 95th percentile of fgap speeds
%                      2 - using unimodal fgap thresholding
%                      3 - using unimodal fgap thresholding corrected for 
%                      comet latency
%                      4 - using the fluctuation radius
%
% OUTPUT
% dataMatMerge        : matrix where fgaps that should be reclassified (as
%                       growth, because their speeds are
%                       >=70% the speed at the end of the growth phase
%                       prior) are consolidated with flanking growth
%                       subtracks. also bgaps that should be reclassified
%                       as pause, because their speeds are lower than the
%                       95th percentile of fgap speeds (after reclass),
%                       have their indices changed to 2.
% dataMatReclass      : matrix where fgaps that should be reclassified are
%                       not consolidated but the track type is changed to
%                       5, and bgaps that should be reclassified have their
%                       track type changed to 6.
% dataMatCrpSecMic    : like dataMatMerge except the growth phases and
%                       linked fgaps/bgaps from the first and last frames
%                       of the movie have been removed. also, column 6
%                       represents lifetime in seconds and 7 represents
%                       displacement in microns. all speeds and displacements
%                       are positive (which makes more sense anyway for
%                       these measurements)
% percentFgapsReclass : percentage of fgaps that get reclassified as
%                       continuation of growth
% percentBgapsReclass : percentage of bgaps that get reclassified as
%                       pause
%
% 
% this function is called by:
% plusTipPostTracking
% plusTipPoolGroupData
% plusTipParamPlot
% plusTipGetSubtrackCoords
% plusTipSpeedMovie
%
% Copyright (C) 2018, Danuser Lab - UTSouthwestern 
%
% This file is part of u-track.
% 
% u-track is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% u-track is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with u-track.  If not, see <http://www.gnu.org/licenses/>.
% 
% 

% Check input
ip = inputParser;
ip.addOptional('fgapReclassScheme', 1, @(x) ismember(x,1:4));
ip.addOptional('bgapReclassScheme', 1, @(x) ismember(x,1:4));
ip.parse(varargin{:});
fgapReclassScheme = ip.Results.fgapReclassScheme;
bgapReclassScheme = ip.Results.bgapReclassScheme;

%% Specify Reclassification Schemes (eventually put into input of function)
% could make this a switch instead of binary and if statements
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OPTIONS FOR FGAP --> UNDETECTED GROWTH RECLASSIFICATION

localFramesBeforeGap = (fgapReclassScheme == 1); 
% Kathyryn's original scheme: Compares gap 
% velocity to the velocity of the last 2-3 frames of the growth just before
% Problem with this scheme is that 
% if the velocity of the comet significantly slows right before pausing
% (will especially happen if there is a delay in the dissociation of
% the comet upon a pause) potentially viable pause information may
% become reclassified as an undetected growth event.  
%NOTE it appears the old scheme seems to have a bit of a bug in that if a 
    % bgap proceeds a short growth subtrack it will result in those frame
    % to frame velocities being included in the beforeFgapSpeed estimation
    % making this estimate negative and thus the fgap will automatically be
    % reclassified. (MB: 03/11)


localFullGrowthSubtrack = (fgapReclassScheme == 2); 
% easiest without artifacts to just use full 
% growth subtrack for now, in the end would like to eliminate the last 2-3 frames 
% before the pause event (the latency time associated with dissociation of
% the comet)

unimodalReclassSingleProj = (fgapReclassScheme == 3);  % if 1, perform unimodalReclassification 
% of fgaps per Project


unimodalReclassPool = (fgapReclassScheme == 4); % if 1, skip reclassification step entirely until 
% have proceeded with the tracking of all projects. 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OPTIONS FOR BGAP-->PAUSE RECLASSIFICATION

bgap95thPercFGapSpeeds = (bgapReclassScheme == 1); % old reclassification scheme

bgapUniModeThreshCorrect = (bgapReclassScheme == 2); % if 1 will correct for comet latency 
% typically more shrinkage events are maintained 

bgapUniModeThreshNoCorrect = (bgapReclassScheme == 3); % use the unimodal thresh as the thresh 
% for the bgap reclass 

bgapFluctRadius = (bgapReclassScheme == 4); % base the bgap based on the fluct radius (get from 
% estimate of fgap displacement) used 2um(could use max value from pause
% data)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   




%% Check Inputs To See If Data is Reclassified Yet
mergeTracks=0;
if nargin<2
    
    % dataMat is the output matrix in projData, where the gaps to
    % consolidate are labeled as trackType=5
    dataMat=projData.nTrack_sF_eF_vMicPerMin_trackType_lifetime_totalDispPix;
    
    % these are the fgaps to consolidate
    growthFgapIdx=find(dataMat(:,5)==5);
    %% Find all fgaps and bgaps 

fgapIdx=find(dataMat(:,5)==2 | dataMat(:,5)==5);
bgapAllIdx = find(dataMat(:,5) == 3 | dataMat(:,5) == 6);
mergeTracks = 1; 
dataMatReclass = dataMat; % same because already reclassified

else % dataMat input and need to perform reclassification with one of the schemes below
    
%% Local FGap Reclassifications 
%% Find all fgaps and bgaps 

fgapIdx=find(dataMat(:,5)==2 | dataMat(:,5)==5);
bgapAllIdx = find(dataMat(:,5) == 3 | dataMat(:,5) == 6);

if (localFullGrowthSubtrack == 1 || localFramesBeforeGap == 1)
    
    beforeFgapIdx=fgapIdx-1; % get index of the growth subtrack before fgap
    eF=dataMat(beforeFgapIdx,3);% get the end frame of subtrack right before pause
    beforeFgapSpeed=zeros(length(fgapIdx),1); % initiate matrix  
   
  if localFullGrowthSubtrack == 1
      
        for iGap = 1:length(fgapIdx)
           beforeFgapSpeed(iGap) = dataMat(beforeFgapIdx(iGap),4); % just use the average 
           % frame-to-frame velocity of the entire growth subtrack for now
        end % end for iGap
        
        projData.fgapReclassScheme = 'Local Scheme: FullGrowth Subtrack Velocity';
  end 
       
 if localFramesBeforeGap == 1
    
        for iGap=1:length(fgapIdx)
            idxRange=[max(1,eF(iGap)-3):eF(iGap)-1]; % not sure why she uses max here
            beforeFgapSpeed(iGap)=nanmean(projData.frame2frameVel_micPerMin(dataMat(beforeFgapIdx(iGap),1),idxRange)); 
        end 
        
        projData.fgapReclassScheme = 'Local Scheme: Velocity 2-3 Frames Before Pause';
       
        
 end % if localFramesBeforeGap
 
    % these are the fgaps to consolidate
    growthFgapIdx=fgapIdx(dataMat(fgapIdx,4)>0.7.*beforeFgapSpeed);
    projData.cutOffValueFGap = NaN; % no global thresh because local scheme
    
end %  if either local scheme


%% Unimodal Thesholding Reclassification Scheme for FGAPS 
% Calculates UniModal Threshold Based on Data In Current Project:
% Note for some cells the sampling may not be sufficient and the data
% should be pooled among all cells of the same condition for that day
% check histogram

if unimodalReclassSingleProj == 1  
    
    fgapSpeeds = dataMat(fgapIdx,4);
           
    % perform unimodal thresholding using all fgap speeds to obtain 
    % maximum fgap speed corresponding to a pause event (above this thresh
    % are likely undetected growth events)
    [cutoffIdx, cutoffValueFGap, sp, axesH,maxBinValue] = cutFirstHistMode(fgapSpeeds,1);
     
   
    % Mark the fgaps for reclassification to growth under unimodal scheme
    
    growthFgapIdx = find(dataMat(:,4) > cutoffValueFGap & dataMat(:,5) == 2);
   
    projData.cutOffValueFGap_VelMicPerMin = cutoffValueFGap;
    projData.maxBinValue_VelMicPerMin = maxBinValue;
    projData.fgapReclassScheme = 'Unimodal Thresholding Per Project';
    
end

%% Unimodal Reclass Pool: Do not perform reclassification just document scheme
%  we will perform reclassification after all data has been pooled 

if unimodalReclassPool == 1
    projData.fgapReclassScheme = 'No Reclassification Until Pool Data';
    projData.bgapReclassScheme = 'No Reclassification Until Pool Data'; 
   
end

%% Record Fgap Reclassifications
dataMatReclass = dataMat;

if unimodalReclassPool ~= 1
dataMatReclass(growthFgapIdx,5) = 5; 
end

%% Old Reclassification Scheme for BGaps Based on fGap Velocity Distribution: Makes Assumption that Comet Latency Negligible
 % Use by default bgapUniModeThreshCorrectForLatwith with either local reclassification scheme
 
if bgap95thPercFGapSpeeds == 1  
    
    fgapMaxSpeed=prctile(dataMat(dataMat(:,5)==2,4),95);

    bgap2pauseIdx=find((dataMat(:,5)==3 | dataMat(:,5)==6) & abs(dataMat(:,4))< ...
        fgapMaxSpeed) ;

    projData.bgapReclassScheme = 'bGapThresh = 95th Percentile of fGapVel';
    projData.cutOffValueBGap_VelMicPerMin = fgapMaxSpeed;
end % if bgap95thPercFGapSpeeds


%% Bgap Reclassification Using Unimodal Fgap Thresh

if bgapUniModeThreshNoCorrect == 1
    
    bgap2pauseIdx = find((dataMat(:,5) == 3 | dataMat(:,5) == 6)...
        & abs(dataMat(:,4)) < cutoffValueFGap);
    
    projData.bgapReclassScheme = 'Unimodal Fgap Thresh: No Correct'; 
    projData.cutOffValueBGap_VelMicPerMin = cutoffValueFGap;
    
end     
    
%%  Bgap Reclassification Using "Corrected" Unimodal Fgap Thesh

if bgapUniModeThreshCorrect == 1 
    
    widthDist = cutoffValueFGap - maxBinValue;
    
    bgapThresh = widthDist - maxBinValue;
    
    if bgapThresh < 0 
        bgapThresh = 0 ; % if the width of the distribution is greater than 
                         % than the width from 0 to max peak value the 
                         % entire distribution of fgaps is in the positive 
                         % range therefore there should be NO bgap
                         % reclassifications
    end
    
    bgap2pauseIdx = find((dataMat(:,5) == 3 | dataMat(:,5) ==6)...
        & abs(dataMatReclass(:,4)) < bgapThresh);
    
    projData.bgapReclassScheme = 'Unimodal Fgap Thresh: Correct For Comet Latency';
    projData.cutOffValueBGap_VelMicPerMin = bgapThresh;
    projData.maxBinValue_VelMicPerMin = maxBinValue;
    projData.fgapDispWidth_VelMicPerMin = widthDist;

end 
    
%% Bgap Reclassification Using Fluctuation Radius

if bgapFluctRadius == 1
    
    bgap2pauseIdx = find((dataMat(:,5) == 3 | dataMat(:,5) == 6)...
        & abs(dataMat(:,7)) < projData.trackingParameters.fluctRadius);
    
   projData.bgapReclassScheme = 'Fluctuation Radius Displacement Cut-off';
   projData.cutOffValueBGapDispPixels = projData.trackingParameters.fluctRadius;
   
end 
    
    
%% Record Bgap Reclassifications
if unimodalReclassPool  ~=  1
dataMat(bgap2pauseIdx,5)=2; % reassign dataMat to have type 2 
dataMatReclass(bgap2pauseIdx,5)=6; % reassign reclassified matrix to have new type of 6
end


%% Calculate Fraction of SubTracks Reclassified 
% fraction of bgaps that are changed to pause, because their speeds are
% slower than the cutoff for fgap pauses

if unimodalReclassPool ~= 1
    
if isempty(bgapAllIdx)
    projData.percentBgapsReclass=NaN;
else
    projData.percentBgapsReclass=100*length(bgap2pauseIdx)/length(bgapAllIdx);
end

% fraction of fgaps that were consolidated into growth, since their speeds
% were more than 50% of the growth speed just prior to the gap
if isempty(fgapIdx)
    projData.percentFgapsReclass=NaN;
else
    projData.percentFgapsReclass=100*length(growthFgapIdx)/length(fgapIdx);
end

else 
    projData.percentFgapsReclass = 0 ;
    projData.percentBgapsReclass = 0;
end









end % if nargin < 2

%% Merge Reclassifications (ie subTrack ID 5 --> 1 need to incorporate reclassified pauses in growth stats)
% merging will be performed regardless of reclassification scheme unless we
% are pooling data

if unimodalReclassPool ~= 1 || mergeTracks == 1
    
    
% these are the affected track numbers for growth fgaps
tracks2check=unique(dataMat(growthFgapIdx,1));
rows2remove=[];

for i=1:length(tracks2check)
    % rows of dataMat corresponding to i track
    subIdx=find(dataMat(:,1)==tracks2check(i));
    % rows of dataMat corresponding to fgaps in track to consolidate
    fgap2remIdx=intersect(growthFgapIdx,subIdx);
    % rows of dataMat corresonding to bgaps or real pauses in track
    sepIdx=union(subIdx(dataMat(subIdx,5)==3),setdiff(intersect(subIdx,fgapIdx),fgap2remIdx))';
    % split the track based on bgaps, so that all fgaps that
    % should be consolidated together can be done dataMat the same time
    sIdx=[subIdx(1); sepIdx(:)];
    eIdx=[sepIdx(:); subIdx(end)];
    % loop through groups of subtracks (split by bgaps)
    for j=1:length(sIdx)
        % pTemp contains the ones to consolidate in this section
        pTemp=intersect(fgap2remIdx,sIdx(j):eIdx(j));
        if ~isempty(pTemp)
            fIdx=min(pTemp)-1; % first row - prior to first fgap
            lIdx=max(pTemp)+1; % last row - after final fgap

            dataMat(fIdx,3)=dataMat(lIdx,3); % change end frame
            dataMat(fIdx,6)=sum(dataMat(fIdx:lIdx,6)); % sum lifetimes
            dataMat(fIdx,7)=sum(dataMat(fIdx:lIdx,7)); % sum total displacements
            dataMat(fIdx,4)= pixPerFrame2umPerMin(dataMat(fIdx,7)/dataMat(fIdx,6),projData.secPerFrame,projData.pixSizeNm); % find new average velocity

            % keep track of which are the extra rows
            rows2remove=[rows2remove fIdx+1:lIdx];
        end
    end
end


% remove the extra rows
dataMat(rows2remove,:)=[];

% NOTE: dataMat now has all fgaps that are reclassified merged
% with their preceding growth subtrack.  Therefore, the indexing of 
% dataMat and dataMatReclass (before merging) will be different!

end % if unimodalReclassPool


%% Conversions and Calculation of Comet Latency
% calculate necessary avg values from above merged data (where all 
% reclassified pauses have been merged with preceding growth 
% subtrack
% Data with pauses reclassified as growth merged in the growth velocity 
% stats (for output) Save this before conversion
dataMatMerge=dataMat;

% perform lifetime and displacement unit conversions
dataMat(:,6)=dataMat(:,6).* projData.secPerFrame; % convert lifetimes to seconds
dataMat(:,7)=dataMat(:,7).*(projData.pixSizeNm/1000); % convert displacements to microns



%avgVelGrowth = mean(dataMat(dataMat(:,5) == 1,4));

%avgDispPause = mean(dataMat(dataMat(:,5) == 2,7)); % in microns
%avgDispPauseBeforeBgapReclass = mean(dataMatReclass((dataMatReclass(:,5) == 2 | dataMatReclass(:,5) == 5),7).*(projData.pixSizeNm/1000));
%absAvgDispPause = mean(abs(dataMat(dataMat(:,5) == 2,7))); % in microns consider bgaps positive


% calc avg latency of comet formation  
%avgLat = avgDispPause/avgVelGrowth; % in minutes
%projData.avgComLatSec = avgLat*60; % in seconds

%avgLatCalcBeforeBgapReclass = avgDispPauseBeforeBgapReclass/avgVelGrowth;
%projData.avgComLatCalcBeforeBgapReclass = avgLatCalcBeforeBgapReclass*60;

%absAvgLat = absAvgDispPause/avgVelGrowth;
%projData.avgComLatSecAbs = absAvgLat*60;
 

%% Remove Growths Initiated in First Frame or Ending in Last Frame From Stats
% do this so one does not bias growth lifetime/displacement data (might not
% be what we want for 

subIdx2rem=[];
% get index of growth and following fgap or bgap (if it exists) that
% begin in the first frame

sF = projData.detectionFrameRange(1,1);
eF = projData.detectionFrameRange(1,2);

%sF=min(dataMat(:,2));

% compound track IDs of all subtracks with nminimum starting frame number

fullIdx2rem=unique(dataMat(dataMat(:,2)==sF,1)); 
for iTr=1:length(fullIdx2rem)
    subIdx=find(dataMat(:,1)==fullIdx2rem(iTr));
    if (length(subIdx)>1) && (dataMat(subIdx(2),5) > 1) % if there is a forward backward gap linked to growth in first frame
        subIdx2rem=[subIdx2rem; subIdx(1:2)]; % Don't remove entire compound track only the fgap or bgap it's linked to
    else
        subIdx2rem=[subIdx2rem; subIdx(1)];
    end
end

% get index of growth and preceeding fgap or bgap
% (if it exists) that end in the last frame

%eF=max(dataMat(:,3));
fullIdx2rem=unique(dataMat(dataMat(:,3)==eF,1));
for iTr=1:length(fullIdx2rem)
    subIdx=find(dataMat(:,1)==fullIdx2rem(iTr));
    if (length(subIdx)>1) && (dataMat(subIdx(end-1),5) > 1)
        subIdx2rem=[subIdx2rem; subIdx(end-1:end)]; % take out the last two of list  
    else
        subIdx2rem=[subIdx2rem; subIdx(end)];
    end
end
% remove both classes for statistics
dataMat(subIdx2rem,:)=[];
dataMatCrpSecMic=dataMat; % NOTE: Kathyrn makes these all absolute values 
% I think for stats it is better to keep sign (MB) 


