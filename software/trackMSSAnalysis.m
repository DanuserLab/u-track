function [trackClass,mssSlope,genDiffCoef,scalingPower,normDiffCoef] ...
    = trackMSSAnalysis(tracks,probDim,momentOrders,alphaMSS)
%TRACKMSSANALYSIS classifies trajectories based on their moment scaling spectrum
%
% Copyright (C) 2020, Danuser Lab - UTSouthwestern 
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

%SYNPOSIS [trackClass,mssSlope,genDiffCoef,scalingPower,normDiffCoef] ...
%    = trackMSSAnalysis(tracks,probDim,momentOrders,alphaMSS)
%
%INPUT  tracks      : Matrix indicating the positions and amplitudes of the
%                     tracked features. Number of rows = number of tracks,
%                     number of columns = 8*number of time points.
%                     Each row consists of
%                     [x1 y1 z1 a1 dx1 dy1 dz1 da1 x2 y2 z2 a2 dx2 dy2 dz2 da2 ...].
%                     NaN is used to indicate time points where the track
%                     does not exist.
%       probDim     : Problem dimensionality. Optional. Default: 2.
%       momentOrders: Orders of moments to be calculated.
%                     Optional. Default: 0 through 6.
%       alphaMSS    : Alpha-value for classification. Can take the
%                     following values:
%                     OLD SCHEME: 0.2, 0.1, 0.05 or 0.01. In this scheme,
%                     the alpha-value is the false positive rate from a
%                     free diffusion-centric perspective.
%                     NEW SCHEME (RECOMMENDED): -0.2, -0.1, -0.05 or -0.01.
%                     See Remark #5 below for explanation and more details.
%                     Optional. Default: 0.1.
%
%OUTPUT trackClass  : # tracks x 1 vector of track classification.
%                     Values mean the following ...
%                     0 = Immobile.
%                     1 = Confined Brownian.
%                     2 = Pure Brownian.
%                     3 = Brownian with drift (directed).
%                     NaN = not classified.
%       mssSlope    : # tracks x 1 vector of each track's slope of the line
%                     representing moment scaling power vs. moment order.
%                     NaN indicates tracks that could not be analyzed.
%       genDiffCoef : # tracks x # orders array of generalized diffusion
%                     coefficients for every moment order considered.
%                     NaN indicates tracks that could not be analyzed.
%       scalingPower: # tracks x # orders array of powers with which moment
%                     values scale with time.
%                     NaN indicates tracks that could not be analyzed.
%       normDiffCoef: # tracks x 1 vector of each track's "normal"
%                     diffusion coefficient.
%                     NaN indicates tracks that could not be analyzed.
%
%REMARKS
%
%(1) Algorithm is based on Ewers et al. 2005. PNAS 102: 15110-15115 and
%Ferrari et al. 2001. Physica D 154: 111-137.
%
%(2) Analysis assumes that there are no kinks in the moment scaling
%spectrum curve, i.e. that the motion is strongly self-similar. Weakly
%self-similar processes will generate an MSS which is piece-wise
%continuous, hence before fitting to estimate the slope the curve must be
%chopped into smaller straight-line pieces (but this is not done).
%
%(3) MSS slope thresholds for confined and directed:
%Those corresponding to alphaMSS = 0.1, 0.05 and 0.01 in 2D case are
%calculated using a smoothing spline fit to the corresponding percentiles,
%which were in turn derived from a sample of 35500 simulations. 
%For all other conditions, the threshold is calculated using
%simple piece-wise line fit to the corresponding percentiles which were
%derived from a sample of 2000 simulations.
%
%(4) MSS slope threshold for immobile:
%Relevant only for OLD SCHEME (free-diffusion centric error rate (Remark 3 above).
%Taken as the threshold for confined divided by 4. Validated
%visually on 2D data only. 
%
%(5) NEW CLASSIFICATION SCHEME (Tony Vega): 
% New scheme considers not only MSS value distributions resulting from
% simulations of freely diffusing tracks, but also looks at confined,
% immobile and directed simulations. Threshold values are chosen to
% minimize the error rate for adjacent distributions (i.e.
% confined/immobile, free/confined and directed/free). Thus the error rate
% for free is no longer fixed as in the old scheme, but depends on track
% length.
% Fully implemented for 2D data. 
% For 1D data, new scheme is currently used to distinguish between
% immobile, confined and free, but not yet implemented for directed. 
% Not implemented at all for 3D data yet. In the works.
%
%Khuloud Jaqaman, March 2008
%                 Feb 2015 : Expanded to include immobile.
% Tony Vega,      July 2016: New option to use thresholds based on
%                            distributions of all diffusion types

%% input

if nargin < 1
    disp('--trackMSSAnalysis: Please input tracks to analyze.');
    return
end

if nargin < 2 || isempty(probDim)
    probDim = 2;
end

if nargin < 3 || isempty(momentOrders)
    momentOrders = 0 : 6;
end
numOrders = length(momentOrders);

if nargin < 4 || isempty(alphaMSS)
    alphaMSS = 0.1;
end

%get number of tracks and frames
[numTracks,numFramesMovie] = size(tracks);
numFramesMovie = numFramesMovie / 8;

%find indices of tracks that are >= 20 frames long - do not attempt
%to calculate moments for shorter tracks
criteria.lifeTime.min = 20;
indx4diff = chooseTracks(tracks,criteria);
clear criteria

%% alpha-value for classification

%determine threshold based on alpha-value and dimensionality

if alphaMSS < 0 %NEW SCHEME

    switch probDim
        case 2
            load('newDiffTypeThreshold.mat');
            mssThreshImm = Mdl.mssThreshImm;
            mssThreshNeg = Mdl.mssThreshNeg;
            mssThreshPos = Mdl.mssThreshPos;
            mssThreshImm = [mssThreshImm(1:min(500,numFramesMovie)); mssThreshImm(500)*ones(max(0,numFramesMovie-500),1)];
            mssThreshNeg = [mssThreshNeg(1:min(500,numFramesMovie)); mssThreshNeg(500)*ones(max(0,numFramesMovie-500),1)];
            mssThreshPos = [mssThreshPos(1:min(500,numFramesMovie)); mssThreshPos(500)*ones(max(0,numFramesMovie-500),1)];

        case 1
            load('newDiffTypeThreshold1D.mat');
            mssThreshImm = Mdl1D.mssThreshImm;
            mssThreshNeg = Mdl1D.mssThreshNeg;
            mssThreshPos = Mdl1D.mssThreshPos;
            mssThreshImm = [mssThreshImm(1:min(500,numFramesMovie)); mssThreshImm(500)*ones(max(0,numFramesMovie-500),1)];
            mssThreshNeg = [mssThreshNeg(1:min(500,numFramesMovie)); mssThreshNeg(500)*ones(max(0,numFramesMovie-500),1)];
            mssThreshPos = [mssThreshPos(1:min(500,numFramesMovie)); mssThreshPos(500)*ones(max(0,numFramesMovie-500),1)];

        case 3
            disp('--trackMSSAnalysis: Optimal thresholds not implemented yet for 3D data.')
            return
    end
else %OLD SCHEME
    switch probDim
        case 1
            switch alphaMSS
                case 0.2 %10th percentile and 90th percentile
                    [mssThreshNeg,mssThreshPos] = threshMSS1D_p20(numFramesMovie);
                case 0.1 %5th percentile and 95th percentile
                    [mssThreshNeg,mssThreshPos] = threshMSS1D_p10(numFramesMovie);
                case 0.05 %2.5th percentile and 97.5th percentile
                    [mssThreshNeg,mssThreshPos] = threshMSS1D_p05(numFramesMovie);
                case 0.01 %0.5th percentile and 99.5th percentile
                    [mssThreshNeg,mssThreshPos] = threshMSS1D_p01(numFramesMovie);
            end
        case 2
            switch alphaMSS
                case 0.2 %10th percentile and 90th percentile
                    [mssThreshNeg,mssThreshPos] = threshMSS2D_p20(numFramesMovie);
                case 0.1 %5th percentile and 95th percentile
                    [mssThreshNeg,mssThreshPos] = threshMSS2D_p10(numFramesMovie);
                case 0.05 %2.5th percentile and 97.5th percentile
                    [mssThreshNeg,mssThreshPos] = threshMSS2D_p05(numFramesMovie);
                case 0.01 %0.5th percentile and 99.5th percentile
                    [mssThreshNeg,mssThreshPos] = threshMSS2D_p01(numFramesMovie);
            end
        case 3
            switch alphaMSS
                case 0.2 %10th percentile and 90th percentile
                    [mssThreshNeg,mssThreshPos] = threshMSS3D_p20(numFramesMovie);
                case 0.1 %5th percentile and 95th percentile
                    [mssThreshNeg,mssThreshPos] = threshMSS3D_p10(numFramesMovie);
                case 0.05 %2.5th percentile and 97.5th percentile
                    [mssThreshNeg,mssThreshPos] = threshMSS3D_p05(numFramesMovie);
                case 0.01 %0.5th percentile and 99.5th percentile
                    [mssThreshNeg,mssThreshPos] = threshMSS3D_p01(numFramesMovie);
            end
    end
    mssThreshImm = mssThreshNeg / 4;
end

%% memory for trajectory classification

%classification means ...
%0 = immobile
%1 = confined Brownian
%2 = pure Brownian
%3 = drift/directed
%NaN = unclassified

trackClass = NaN(numTracks,1);
mssSlope = NaN(numTracks,1);
genDiffCoef = NaN(numTracks,numOrders);
scalingPower = NaN(numTracks,numOrders);
normDiffCoef = NaN(numTracks,1);

%% moments and their scaling with time

for iTrack = indx4diff'

    %get track start and end time
    trackSEL = getTrackSEL(tracks(iTrack,:));
    startTime = trackSEL(1);
    endTime = trackSEL(2);
    numTimePoints = trackSEL(3);

    %extract track's coordinates and their standard deviations
    coordinates = [tracks(iTrack,1:8:end)' tracks(iTrack,2:8:end)' tracks(iTrack,3:8:end)'];
    coordinates = coordinates(startTime:endTime,:);
    standardDevs = [tracks(iTrack,5:8:end)' tracks(iTrack,6:8:end)' tracks(iTrack,7:8:end)'];
    standardDevs = standardDevs(startTime:endTime,:);

    %define maximum time lag for moment calculation
    maxLag = min(30,floor(numTimePoints/4));

    %calculate track moments
    trackMomentsT = calcTrackMoments(coordinates,standardDevs,momentOrders,maxLag);
    trackMoments = [trackMomentsT.momentValues];
    trackMoments = trackMoments(:,1:2:end);
    
    %estimate the moment scaling spectrum (MSS),
    %i.e. the scaling power for all moments
    scalingPowerT = NaN(1,numOrders);
    genDiffCoefT = NaN(1,numOrders);
    for iOrder = 1 : length(momentOrders)

        %caculate ln(lag) and ln(moment)
        lnTime = log((1:maxLag)');
        lnMoment = log(trackMoments(1:maxLag,iOrder));

        
        %remove any NaNs     
        indxGood = find(~isnan(lnMoment));
        lnTime = lnTime(indxGood);
        lnMoment = lnMoment(indxGood);
        
        %if there are moments to fit ...
        if length(lnMoment) > 1

            %fit a straight line in the plot of lnMoment vs. lnTime
            [slParam,sFull] = polyfit(lnTime,lnMoment,1);

                  scalingPowerT(iOrder) = slParam(1);
                  genDiffCoefT(iOrder) = exp(slParam(2)) / 2 / probDim;  

            %if this is the 2nd moment, calculate the "normal" diffusion
            %coefficient
            if momentOrders(iOrder)==2
                options = optimset('Display','off','Jacobian','on');
                lnSlope = lsqcurvefit(@strLineFun2,1,lnTime(1:min(5,...
                    length(lnTime))),lnMoment(1:min(5,length(lnMoment))),...
                    [],[],options);
                normDiffCoefT = exp(lnSlope) / 2 / probDim;                
            end
            
        end

    end
    
    
    %keep only non-NaN scaling powers
    indxGood = find(~isnan(scalingPowerT));
    momentOrders4fit = momentOrders(indxGood);
    scalingPowerT = scalingPowerT(indxGood);
    genDiffCoefT = genDiffCoefT(indxGood);

    %if there are non-NaN scaling powers
    if ~isempty(scalingPowerT)

        %fit a straight line to the MSS
        slParam = polyfit(momentOrders4fit,scalingPowerT,1);
%         slMoment2 = polyfit((1:maxLag)',trackMoments(1:maxLag,3),1);
%         plot((1:maxLag)',trackMoments(1:maxLag,3),'b');
        %get the slope of the line
        mssSlopeT = slParam(1);

        %classify track as ...
        %0 = immobile, if mSS slope <= mssThreshImm
        %1 = confined Brownian, if mssThreshImm < MSS slope < mssThreshNeg
        %2 = pure Brownian, if mssThreshNeg <= MSS slope <= mssThreshPos
        %3 = directed, if MSS slope > mssThreshPos
        
        if ~isnan(mssSlopeT)
            
                if mssSlopeT <= mssThreshImm(numTimePoints)
                    trackClass(iTrack) = 0;
                elseif mssSlopeT < mssThreshNeg(numTimePoints) && mssSlopeT > mssThreshImm(numTimePoints)
                    trackClass(iTrack) = 1;
                elseif mssSlopeT > mssThreshPos(numTimePoints)
                    trackClass(iTrack) = 3;
                else
                    trackClass(iTrack) = 2;
                end
        end
        
        %save additional output information
        mssSlope(iTrack) = mssSlopeT;
        genDiffCoef(iTrack,:) = genDiffCoefT;
        scalingPower(iTrack,:) = scalingPowerT;
        normDiffCoef(iTrack) = normDiffCoefT;
        
    end

end

%% subfunction 1


function [y,d] = strLineFun2(logSlope,x)

y = logSlope + x;
d = ones(size(x));


%% thresholds

% function mssThreshImm = threshMSSImm2D_p20(nTP)
% 
% %2D, alpha = 0.2 - UPDATE
% 
% %approximate threshold curve via the equation y = ax^b+c fitted once to
% %trajectory lengths of 20-130 and once to trajectory lengths of 130-500
% %take the max at any trajectory length, this gives best fit to data
% %this yields an inflection point around 120-125.
% %parameters obtained via curve fitting toolbox from 24000 simulations
% a1 = 3.543;
% b1 = -1.326;
% c1 = 0.003187;
% a2 = 0.09935;
% b2 = -0.5089;
% c2 = 0.0006031;
% 
% mssThreshImm = NaN(nTP,1);
% mssThreshImm(20:end) = max([a1*(20:nTP)'.^b1 + c1,a2*(20:nTP)'.^b2 + c2],[],2);
% if nTP > 500
%     mssThreshImm(500:end) = mssThreshImm(500);
% end

%%%%%%%%%%%%%%%%%%%%%%%

% function mssThreshImm = threshMSSImm2D_p10(nTP)
% 
% %2D, alpha = 0.1 - UPDATE
% 
% %approximate threshold curve via the equation y = ax^b+c fitted once to
% %trajectory lengths of 20-130 and once to trajectory lengths of 130-500
% %take the max at any trajectory length, this gives best fit to data
% %this yields an inflection point around 120-125.
% %parameters obtained via curve fitting toolbox from 24000 simulations
% a1 = 5.069;
% b1 = -1.284;
% c1 = 0.004329;
% a2 = 0.1554;
% b2 = -0.4845;
% c2 = -0.0002119;
% 
% mssThreshImm = NaN(nTP,1);
% mssThreshImm(20:end) = max([a1*(20:nTP)'.^b1 + c1,a2*(20:nTP)'.^b2 + c2],[],2);
% if nTP > 500
%     mssThreshImm(500:end) = mssThreshImm(500);
% end


%%%%%%%%%%%%%%%%%%%%%%%

% function mssThreshImm = threshMSSImm2D_p05(nTP)
% 
% %2D, alpha = 0.05 - UPDATE
% 
% %approximate threshold curve via the equation y = ax^b+c fitted once to
% %trajectory lengths of 20-130 and once to trajectory lengths of 130-500
% %take the max at any trajectory length, this gives best fit to data
% %this yields an inflection point around 120-125.
% %parameters obtained via curve fitting toolbox from 24000 simulations
% a1 = 6.082;
% b1 = -1.251;
% c1 = 0.005012;
% a2 = 0.2191;
% b2 = -0.4921;
% c2 = -0.0006046;
% 
% mssThreshImm = NaN(nTP,1);
% mssThreshImm(20:end) = max([a1*(20:nTP)'.^b1 + c1,a2*(20:nTP)'.^b2 + c2],[],2);
% if nTP > 500
%     mssThreshImm(500:end) = mssThreshImm(500);
% end

%%%%%%%%%%%%%%%%%%%%%%%

function mssThreshImm = threshMSSImm2D_p01(nTP)

%2D, alpha = 0.01

%approximate threshold curve via the equation y = ax^b+c fitted once to
%trajectory lengths of 20-130 and once to trajectory lengths of 130-500
%take the max at any trajectory length, this gives best fit to data
%this yields an inflection point around 120-125.
%parameters obtained via curve fitting toolbox from 24000 simulations
a1 = 9.252;
b1 = -1.263;
c1 = 0.00962;
% a2 = 0.2969;
% b2 = -0.4492;
% c2 = -0.004033;

mssThreshImm = NaN(nTP,1);
mssThreshImm(20:end) = a1*(20:nTP)'.^b1 + c1;
if nTP > 500
    mssThreshImm(500:end) = mssThreshImm(500);
end

%%%%%%%%%%%%%%%%%%%%%%%

function mssThreshImm = threshMSSImm2D_p001(nTP)

%2D, alpha = 0.001

%approximate threshold curve via the equation y = ax^b+c fitted once to
%trajectory lengths of 20-130 and once to trajectory lengths of 130-500
%take the max at any trajectory length, this gives best fit to data
%this yields an inflection point around 120-125.
%parameters obtained via curve fitting toolbox from 24000 simulations
a1 = 11.65;
b1 = -1.239;
c1 = 0.01412;
% a2 = 0.2969;
% b2 = -0.4492;
% c2 = -0.004033;

mssThreshImm = NaN(nTP,1);
mssThreshImm(20:end) = a1*(20:nTP)'.^b1 + c1;
if nTP > 500
    mssThreshImm(500:end) = mssThreshImm(500);
end

%%%%%%%%%%%%%%%%%%%%%%%

function mssThreshImm = threshMSSImm2D_p0001(nTP)

%2D, alpha = 0.0001

%approximate threshold curve via the equation y = ax^b+c 
%parameters obtained via curve fitting toolbox from 50000 simulations
a1 = 13.9;
b1 = -1.229;
c1 = 0.01902;
% a2 = -0.005641;
% b2 = 0.4097;
% c2 = 0.09628;

mssThreshImm = NaN(nTP,1);
mssThreshImm(20:end) = a1*(20:nTP)'.^b1 + c1;
if nTP > 500
    mssThreshImm(500:end) = mssThreshImm(500);
end

%%%%%%%%%%%%%%%%%%%%%%%

function mssThreshImm = threshMSSImm2D_p0(nTP)

%2D, alpha = 0

%approximate threshold curve via the equation y = ax^b+c 
%parameters obtained via curve fitting toolbox from 50000 simulations
a1 = 16.96;
b1 = -1.273;
c1 = 0.03433;
% a2 = -0.005641;
% b2 = 0.4097;
% c2 = 0.09628;

mssThreshImm = NaN(nTP,1);
mssThreshImm(20:end) = a1*(20:nTP)'.^b1 + c1;
if nTP > 500
    mssThreshImm(500:end) = mssThreshImm(500);
end

%%%%%%%%%%%%%%%%%%%%%%%

function [mssThreshNeg,mssThreshPos] = threshMSS1D_p20(nTP)

%1D, alpha = 0.2

%threshold curve parameters
turnPointsM = [20 60 200 500];
turnPointsP = [20 60 100 500];
slopeM = [0.002434794446981 0.000542795021531 0.000165329365168 0];
slopeP = [-0.001541239220102 -0.00036860800289 -0.00002638219148 0];
interseptM = [0.130233710469096 0.243753675996077 0.319246807268674 0.40191148985285];
interseptP = [0.672499246818346 0.602141373785626 0.567918792644692 0.554727696904506];

%threshold curve evaluation
[mssThreshNeg,mssThreshPos] = getThreshCurve(turnPointsM,turnPointsP,...
    slopeM,slopeP,interseptM,interseptP,nTP);

%%%%%%%%%%%%%%%%%%%%%%%

function [mssThreshNeg,mssThreshPos] = threshMSS1D_p10(nTP)

%1D, alpha = 0.1

%threshold curve parameters
turnPointsM = [20 50 200 500];
turnPointsP = [20 60 150 500];
slopeM = [0.00403865723972 0.00073054296121 0.000167789786604 0];
slopeP = [-0.001863373027824 -0.000207607837099 -0.000068295975124 0];
interseptM = [0.018637170190375 0.184042884115918 0.296593519037011 0.380488412339042];
interseptP = [0.727469044719516 0.628123133276013 0.607226353979748 0.5730783664177];

%threshold curve evaluation
[mssThreshNeg,mssThreshPos] = getThreshCurve(turnPointsM,turnPointsP,...
    slopeM,slopeP,interseptM,interseptP,nTP);

%%%%%%%%%%%%%%%%%%%%%%%

function [mssThreshNeg,mssThreshPos] = threshMSS1D_p05(nTP)

%1D, alpha = 0.05

%threshold curve parameters
turnPointsM = [20 60 200 500];
turnPointsP = [20 50 200 500];
slopeM = [0.004221335524384 0.000794743062117 0.000192611669839 0];
slopeP = [-0.002478698219027 -0.000203422663306 -0.000101280672122 0];
interseptM = [-0.059561585354062 0.146033962381945 0.266460240837543 0.362766075757204];
interseptP = [0.772420366381382 0.658656588595335 0.638228190358568 0.587587854297634];

%threshold curve evaluation
[mssThreshNeg,mssThreshPos] = getThreshCurve(turnPointsM,turnPointsP,...
    slopeM,slopeP,interseptM,interseptP,nTP);

%%%%%%%%%%%%%%%%%%%%%%%

function [mssThreshNeg,mssThreshPos] = threshMSS1D_p01(nTP)

%1D, alpha = 0.01

%threshold curve parameters
turnPointsM = [20 50 150 500];
turnPointsP = [20 45 100 500];
slopeM = [0.00748519566052008 0.00132812830129705 0.000291477165467116 0];
slopeP = [-0.00263804313196125 -0.000652822904393821 -0.000123920022338378 0];
interseptM = [-0.27508360480994 0.0327697631512118 0.188267433525701 0.334006016259259];
interseptP = [0.839817228313045 0.750482318072511 0.697592029866967 0.635632018697778];

%threshold curve evaluation
[mssThreshNeg,mssThreshPos] = getThreshCurve(turnPointsM,turnPointsP,...
    slopeM,slopeP,interseptM,interseptP,nTP);

%%%%%%%%%%%%%%%%%%%%%%%

function [mssThreshNeg,mssThreshPos] = threshMSS2D_p20(nTP)

%2D, alpha = 0.2

%threshold curve parameters
turnPointsM = [20 50 200 500];
turnPointsP = [20 60 150 500];
slopeM = [0.002365010936411 0.000327173962804 0.000137961052668 0];
slopeP = [-0.001158501020215 -0.000170757301414 -0.000033415212947 0];
interseptM = [0.227688919307209 0.329580767987585 0.367423350014689 0.436403876348735];
interseptP = [0.634413262884249 0.575148639756167 0.554547326486174 0.537839720012491];

%threshold curve evaluation
[mssThreshNeg,mssThreshPos] = getThreshCurve(turnPointsM,turnPointsP,...
    slopeM,slopeP,interseptM,interseptP,nTP);

%%%%%%%%%%%%%%%%%%%%%%%

function [mssThreshNeg,mssThreshPos] = threshMSS2D_p10(nTP)

%2D, alpha = 0.1

%NEW

mssThreshNeg = [...
    NaN
    NaN
    NaN
    NaN
    NaN
    NaN
    NaN
    NaN
    NaN
    NaN
    NaN
    NaN
    NaN
    NaN
    NaN
    NaN
    NaN
    NaN
    NaN
    0.22745412174663
    0.233161658076958
    0.238712335849316
    0.243942595943949
    0.248765797042666
    0.25319606019234
    0.257222712435148
    0.260820097125929
    0.264022922051307
    0.266958704365908
    0.269672315280133
    0.272160518632548
    0.274452617591953
    0.276641195820616
    0.278727477064239
    0.280671923536532
    0.28247532444392
    0.284208829102268
    0.285860941991394
    0.287378057585524
    0.288753002802575
    0.290071283927206
    0.291367428204686
    0.292623598509624
    0.293840775024627
    0.295072879740121
    0.296295727221108
    0.297444504094809
    0.298509395809849
    0.299561725921043
    0.30060498942953
    0.30160296929942
    0.302539697275081
    0.30345230271871
    0.304325700625445
    0.305133593757577
    0.305876861532026
    0.306607773063942
    0.307338820167741
    0.308044876026691
    0.308719468362311
    0.30939888182498
    0.310072992600017
    0.310709988314711
    0.311326320570336
    0.311994633765683
    0.312731698246595
    0.313505449420241
    0.314295810475204
    0.315118913163757
    0.315947218279532
    0.316733511082753
    0.317463410977594
    0.31817306411387
    0.318867078363294
    0.319524117937474
    0.320140509673117
    0.320749108687508
    0.321348933053907
    0.321914030863819
    0.32243530687815
    0.322933580726169
    0.323389130856063
    0.323767172308512
    0.324067984556657
    0.324336673775683
    0.324587993066945
    0.324811724125205
    0.32501270487864
    0.325222413638182
    0.325440555615016
    0.32565239061206
    0.325870736964501
    0.32613917517856
    0.326447877924341
    0.3267545200733
    0.327038798351875
    0.32732180403296
    0.327604253128306
    0.32787086458923
    0.328128885523064
    0.328422364979933
    0.328760995825149
    0.329126521994235
    0.329507449821629
    0.329917763166535
    0.330338236348803
    0.330735410709573
    0.331102609831489
    0.331463675658305
    0.331815107768013
    0.332139086364455
    0.332433929405062
    0.332723852877751
    0.333002266913764
    0.333239771054926
    0.333424148275459
    0.333583456328585
    0.33373100341678
    0.333871173017673
    0.334035432605662
    0.33428127941858
    0.334632903253784
    0.335080076975016
    0.33560316901398
    0.336184256673208
    0.336808209446308
    0.337459442116977
    0.338122660077784
    0.338786872271129
    0.339444818853822
    0.340090606417223
    0.340719439276555
    0.341330072164072
    0.341925717534233
    0.342511442719755
    0.343091088011455
    0.343672288496437
    0.344262857441147
    0.344866889856187
    0.345479535539688
    0.346091444559155
    0.346690433756783
    0.347271075584896
    0.34784415398524
    0.348421271526866
    0.349011618630643
    0.349616811684455
    0.350235527269212
    0.350864336548045
    0.351493821882612
    0.352114709737342
    0.352725274065241
    0.353325695040802
    0.353916152838517
    0.354496827632879
    0.355067899598381
    0.355629548909516
    0.356181955740775
    0.356725300266653
    0.357259762661642
    0.357785523100233
    0.358302761756921
    0.358811658806198
    0.359312394422556
    0.359805148780489
    0.360290102054488
    0.360767434419048
    0.361237326048659
    0.361699957117816
    0.362155507801011
    0.362604158272736
    0.363046088707484
    0.363481479279749
    0.363910510164022
    0.364333361534797
    0.364750213566565
    0.365161246433821
    0.365566640311056
    0.365966575372764
    0.366361231793436
    0.366750789747566
    0.367135429409647
    0.367515330954171
    0.367890674555631
    0.368261640388519
    0.368628408627329
    0.368991159446552
    0.369350073020683
    0.369705329524212
    0.370057109131635
    0.370405592017441
    0.370750958356126
    0.371093388322181
    0.371433062090099
    0.371770159834373
    0.372104861729495
    0.372437347949958
    0.372767798670255
    0.373096394064879
    0.373423314308323
    0.373748707425412
    0.374072592842313
    0.374394957835524
    0.374715789681543
    0.37503507565687
    0.375352803038005
    0.375668959101444
    0.375983531123689
    0.376296506381236
    0.376607872150587
    0.376917615708238
    0.37722572433069
    0.377532185294442
    0.377836985875991
    0.378140113351838
    0.37844155499848
    0.378741298092418
    0.37903932991015
    0.379335637728174
    0.379630208822991
    0.379923030471098
    0.380214089948994
    0.38050337453318
    0.380790871500152
    0.381076568126412
    0.381360451688456
    0.381642509462785
    0.381922728725897
    0.382201096754292
    0.382477600824468
    0.382752228212923
    0.383024966196158
    0.38329580205067
    0.38356472305296
    0.383831716479525
    0.384096769606865
    0.384359869711479
    0.384621004069865
    0.384880159958523
    0.385137324653951
    0.385392485432649
    0.385645629571115
    0.385896744345848
    0.386145817033347
    0.386392834910112
    0.386637785252641
    0.386880655337432
    0.387121432440986
    0.3873601038398
    0.387596656810374
    0.387831087681762
    0.388063428993232
    0.38829372233661
    0.388522009303719
    0.388748331486384
    0.38897273047643
    0.38919524786568
    0.38941592524596
    0.389634804209093
    0.389851926346903
    0.390067333251215
    0.390281066513854
    0.390493167726643
    0.390703678481408
    0.390912640369971
    0.391120094984158
    0.391326083915793
    0.391530648756701
    0.391733831098705
    0.391935672533629
    0.3921362146533
    0.392335499049539
    0.392533567314173
    0.392730461039024
    0.392926221815918
    0.39312089123668
    0.393314510893132
    0.3935071223771
    0.393698767280407
    0.393889487194879
    0.394079323712339
    0.394268318424612
    0.394456512923522
    0.394643948800894
    0.394830667648551
    0.395016711058318
    0.39520212062202
    0.395386937931481
    0.395571204578524
    0.395754962154976
    0.395938252252658
    0.396121116463397
    0.396303596379016
    0.39648573359134
    0.396667569692193
    0.3968491462734
    0.397030504926783
    0.397211687244169
    0.397392734817382
    0.397573689238245
    0.397754579905936
    0.397935387449051
    0.398116080303536
    0.39829662690534
    0.398476995690409
    0.398657155094692
    0.398837073554137
    0.399016719504692
    0.399196061382304
    0.39937506762292
    0.39955370666249
    0.39973194693696
    0.399909756882278
    0.400087104934393
    0.400263959529251
    0.400440289102801
    0.400616062090991
    0.400791246929767
    0.400965812055079
    0.401139725902874
    0.401312956909099
    0.401485473509703
    0.401657244140633
    0.401828237237836
    0.401998421237262
    0.402167764574857
    0.402336235686569
    0.402503803008346
    0.402670434976136
    0.402836100025887
    0.403000766593546
    0.403164403115061
    0.403326978026381
    0.403488459763451
    0.403648816762222
    0.40380801745864
    0.403966030288653
    0.404122823688209
    0.404278366093255
    0.40443262593974
    0.404585571663611
    0.404737171700816
    0.404887394487303
    0.405036208459019
    0.405183582051913
    0.405329483701932
    0.405473881845024
    0.405616744917136
    0.405758041354217
    0.405897739592214
    0.406035820637939
    0.406172315781653
    0.406307268884482
    0.406440723807553
    0.40657272441199
    0.406703314558921
    0.40683253810947
    0.406960438924763
    0.407087060865925
    0.407212447794084
    0.407336643570364
    0.407459692055891
    0.407581637111791
    0.40770252259919
    0.407822392379213
    0.407941290312986
    0.408059260261635
    0.408176346086286
    0.408292591648064
    0.408408040808095
    0.408522737427504
    0.408636725367419
    0.408750048488963
    0.408862750653263
    0.408974875721446
    0.409086467554636
    0.409197570013958
    0.40930822696054
    0.409418482255507
    0.409528379759984
    0.409637963335097
    0.409747276841972
    0.409856364141735
    0.409965269095511
    0.410074035564427
    0.410182707409607
    0.410291328492178
    0.410399942673265
    0.410508593813994
    0.410617325775491
    0.410726182418881
    0.410835207605291
    0.410944445195846
    0.411053939051672
    0.411163733033894
    0.411273871003638
    0.411384396822031
    0.411495354350197
    0.411606787449263
    0.411718739980353
    0.411831243447033
    0.411944279922618
    0.41205781912286
    0.412171830763515
    0.412286284560336
    0.412401150229077
    0.41251639748549
    0.412631996045331
    0.412747915624353
    0.412864125938309
    0.412980596702954
    0.41309729763404
    0.413214198447322
    0.413331268858554
    0.413448478583489
    0.413565797337881
    0.413683194837483
    0.41380064079805
    0.413918104935335
    0.414035556965092
    0.414152966603074
    0.414270303565036
    0.414387537566731
    0.414504638323913
    0.414621575552335
    0.414738318967752
    0.414854838285916
    0.414971103222582
    0.415087083493504
    0.415202748814435
    0.415318068901129
    0.415433013469339
    0.41554755223482
    0.415661654913325
    0.415775291220608
    0.415888430872423
    0.416001043584523
    0.416113099072661
    0.416224567052593
    0.416335417240071
    0.41644561935085
    0.416555143100682
    0.416663958205323
    0.416772034380524
    0.416879341342041
    0.416985848805627
    0.417091526487036
    0.41719634410202
    0.417300271366335
    0.417403277995734
    0.417505341923076
    0.417606473949638
    0.417706693093804
    0.417806018373958
    0.417904468808482
    0.41800206341576
    0.418098821214175
    0.418194761222111
    0.41828990245795
    0.418384263940076
    0.418477864686872
    0.418570723716721
    0.418662860048007
    0.418754292699113
    0.418845040688422
    0.418935123034317
    0.419024558755181
    0.419113366869399
    0.419201566395352
    0.419289176351425
    0.419376215756001
    0.419462703627462
    0.419548658984192
    0.419634100844575
    0.419719048226993
    0.41980352014983
    0.419887535631469
    0.419971113690293
    0.420054273344686
    0.42013703361303
    0.42021941351371
    0.420301432065108
    0.420383108285607
    0.420464461193592
    0.420545509807444
    0.420626273145547
    0.420706770226285
    0.420787020068041
    0.420867041689198
    0.420946854108139
    0.421026476343248
    0.421105927412907
    0.421185226335501
    0.421264392129411
    0.421343443813023
    0.421422400404718
    0.42150128092288
    0.421580104385892
    0.421658889812138
    0.42173765622];

mssThreshPos = [...
    NaN
    NaN
    NaN
    NaN
    NaN
    NaN
    NaN
    NaN
    NaN
    NaN
    NaN
    NaN
    NaN
    NaN
    NaN
    NaN
    NaN
    NaN
    NaN
    0.654314617903462
    0.651121703618836
    0.647869171910077
    0.644520555133229
    0.641137482542311
    0.637903321457796
    0.63492029309509
    0.632207407161199
    0.629758811559197
    0.627586008463287
    0.625643972537407
    0.62385997115495
    0.622186987322888
    0.62061323579898
    0.619105337803229
    0.617629145177402
    0.616182400201377
    0.614796525678107
    0.61347629909111
    0.612215523478362
    0.611019193702705
    0.609905622420528
    0.608874132727919
    0.607903724208895
    0.606975549542505
    0.606083268235531
    0.60521550511657
    0.604362118530077
    0.603528316670036
    0.602730247051137
    0.601975132813481
    0.601262510443939
    0.600594673453327
    0.599978970639928
    0.599409048563272
    0.598874413679956
    0.598367490414222
    0.597883655563545
    0.597417910899598
    0.596964849071257
    0.596519436650525
    0.596083431316087
    0.59565268357715
    0.595229222868723
    0.594820562846994
    0.594437303746982
    0.594079665805683
    0.59373928477206
    0.593413843784139
    0.59311040674089
    0.592830286005747
    0.592568178223267
    0.592309735010737
    0.592049584671471
    0.591790698081427
    0.591537469821343
    0.591288939181546
    0.591041605238755
    0.590790505128832
    0.590526977772355
    0.590248990681051
    0.589962657416809
    0.589668751396622
    0.589361684433013
    0.589045030039499
    0.588735585436852
    0.588449920113551
    0.588198347328549
    0.587982671329895
    0.587794838407755
    0.587622139378069
    0.587454418510344
    0.587288748340044
    0.587127936425178
    0.586978202975761
    0.586849116596437
    0.586747276738412
    0.586670983610322
    0.586611463626982
    0.586554850240366
    0.586486621815529
    0.586395315721713
    0.586276014920567
    0.586130100920972
    0.585961396464937
    0.585771514557314
    0.585569696748179
    0.585370077291245
    0.585183527014906
    0.585017817944233
    0.584873594723511
    0.584745878540664
    0.584626492095028
    0.584507760057278
    0.584382802217256
    0.584246470584476
    0.584094844369256
    0.583931155748866
    0.583763127784926
    0.58360097165287
    0.583458975575805
    0.583347272688156
    0.583265070129244
    0.583199753243692
    0.583134817848274
    0.583060248646393
    0.582971390924761
    0.582865257600286
    0.582739519936508
    0.582592287198727
    0.582426755743775
    0.582248928173202
    0.582066830245188
    0.581888987689774
    0.581720541653177
    0.581563871053679
    0.581418683198207
    0.581285068998552
    0.581160731795124
    0.581041256783231
    0.58092047367492
    0.580794812277068
    0.580664881909549
    0.580532664465509
    0.580398490404841
    0.580262229510879
    0.580123980772444
    0.579985663304701
    0.579849611235834
    0.579716586842009
    0.579586073057527
    0.579456652329342
    0.579327983972931
    0.579200063914547
    0.57907288808044
    0.578946452396861
    0.578820752790062
    0.578695785186293
    0.578571545511806
    0.578448029692851
    0.57832523365568
    0.578203153326545
    0.578081784631695
    0.577961123497383
    0.577841165849859
    0.577721907615374
    0.57760334472018
    0.577485473090528
    0.577368288652668
    0.577251787332853
    0.577135965057332
    0.577020817752357
    0.576906341344179
    0.57679253175905
    0.57667938492322
    0.576566896762941
    0.576455063204463
    0.576343880174039
    0.576233343597918
    0.576123449402351
    0.576014193513592
    0.575905571857889
    0.575797580361494
    0.575690214950659
    0.575583471551635
    0.575477346090672
    0.575371834494022
    0.575266932687935
    0.575162636598664
    0.575058942152459
    0.574955845275571
    0.574853341894251
    0.574751427934751
    0.574650099323321
    0.574549351986213
    0.574449181849677
    0.574349584839966
    0.574250556883329
    0.574152093906019
    0.574054191834285
    0.57395684659438
    0.57386005432578
    0.573763812020863
    0.573668116885235
    0.573572966124498
    0.573478356944257
    0.573384286550116
    0.573290752147679
    0.573197750942552
    0.573105280140336
    0.573013336946638
    0.57292191856706
    0.572831022207208
    0.572740645072685
    0.572650784369095
    0.572561437302043
    0.572472601077132
    0.572384272899968
    0.572296449976153
    0.572209129511293
    0.572122308710991
    0.572035984780852
    0.571950154926479
    0.571864816353477
    0.57177996626745
    0.571695601874002
    0.571611720378738
    0.571528318987261
    0.571445394905175
    0.571362945338085
    0.571280967491596
    0.57119945857131
    0.571118415782832
    0.571037836331767
    0.570957717423718
    0.57087805626429
    0.570798850059087
    0.570720096013713
    0.570641791333772
    0.570563933224868
    0.570486518892606
    0.570409545542589
    0.570333010380422
    0.570256910611709
    0.570181243442054
    0.570106006077062
    0.570031195722335
    0.56995680958348
    0.569882844866099
    0.569809298775796
    0.569736168518177
    0.569663450408612
    0.569591137201538
    0.569519220761158
    0.569447692951678
    0.569376545637301
    0.569305770682231
    0.569235359950672
    0.569165305306828
    0.569095598614902
    0.5690262317391
    0.568957196543624
    0.568888484892679
    0.568820088650469
    0.568751999681198
    0.568684209849069
    0.568616711018287
    0.568549495053056
    0.568482553817579
    0.56841587917606
    0.568349462992705
    0.568283297131715
    0.568217373457296
    0.568151683833652
    0.568086220124986
    0.568020974195502
    0.567955937909404
    0.567891103130897
    0.567826461724184
    0.56776200555347
    0.567697726482957
    0.567633616376851
    0.567569667099354
    0.567505870514672
    0.567442218487008
    0.567378702880566
    0.567315315559549
    0.567252048388163
    0.567188893230611
    0.567125841951096
    0.567062886413823
    0.567000018482997
    0.566937230022819
    0.566874512897496
    0.56681185897123
    0.566749260108226
    0.566686708172688
    0.566624195028819
    0.566561712540824
    0.566499252572906
    0.566436806989269
    0.566374369471378
    0.56631194096973
    0.566249524252085
    0.566187122086201
    0.566124737239837
    0.566062372480752
    0.566000030576704
    0.565937714295451
    0.565875426404753
    0.565813169672368
    0.565750946866054
    0.56568876075357
    0.565626614102675
    0.565564509681128
    0.565502450256686
    0.565440438597109
    0.565378477470155
    0.565316569643583
    0.565254717885151
    0.565192924962618
    0.565131193643743
    0.565069526696284
    0.565007926887999
    0.564946396986648
    0.564884939759989
    0.564823557975781
    0.564762254401782
    0.564701031805751
    0.564639892955446
    0.564578840618626
    0.56451787756305
    0.564457006556476
    0.564396230366663
    0.56433555176137
    0.564274973508354
    0.564214498375376
    0.564154129130192
    0.564093868540563
    0.564033719374245
    0.563973684398999
    0.563913766382583
    0.563853968092755
    0.563794292297274
    0.563734741763898
    0.563675319260387
    0.563616027554498
    0.563556869413991
    0.563497847606624
    0.563438964900155
    0.563380224062343
    0.563321628193631
    0.563263181725198
    0.563204889420906
    0.563146756044618
    0.563088786360197
    0.563030985131506
    0.562973357122407
    0.562915907096764
    0.562858639818438
    0.562801560051294
    0.562744672559193
    0.562687982105998
    0.562631493455572
    0.562575211371779
    0.56251914061848
    0.562463285959539
    0.562407652158818
    0.56235224398018
    0.562297066187488
    0.562242123544604
    0.562187420815392
    0.562132962763714
    0.562078754153433
    0.562024799748412
    0.561971104312514
    0.561917672609601
    0.561864509403535
    0.561811619458181
    0.561759007537401
    0.561706678405057
    0.561654636825012
    0.561602887561129
    0.561551435377271
    0.561500285037301
    0.561449441305081
    0.561398908944474
    0.561348692719344
    0.561298797393551
    0.561249227730961
    0.561199988495435
    0.561151084450835
    0.561102520361026
    0.561054300989869
    0.561006431101228
    0.560958915458965
    0.560911758826943
    0.560864965969024
    0.560818541649073
    0.56077249063095
    0.56072681767852
    0.560681525359265
    0.560636607455146
    0.560592055551748
    0.560547861234651
    0.56050401608944
    0.560460511701695
    0.560417339656999
    0.560374491540935
    0.560331958939085
    0.560289733437032
    0.560247806620357
    0.560206170074644
    0.560164815385474
    0.560123734138431
    0.560082917919096
    0.560042358313052
    0.560002046905881
    0.559961975283166
    0.559922135030489
    0.559882517733433
    0.559843114977579
    0.559803918348511
    0.55976491943181
    0.559726109813059
    0.559687481077841
    0.559649024811738
    0.559610732600332
    0.559572596029206
    0.559534606683942
    0.559496756150122
    0.559459036013329
    0.559421437859146
    0.559383953273154
    0.559346573840936
    0.559309291148075
    0.559272096780153
    0.559234982322751
    0.559197939361454
    0.559160959481843
    0.5591240342695
    0.559087155310008
    0.559050314188949
    0.559013502491906
    0.558976711804461
    0.558939933712196
    0.558903159800694
    0.558866381655538
    0.558829590862309
    0.55879277900659
    0.558755937673964
    0.55871905997872
    0.558682145149977
    0.558645193945562
    0.558608207123302
    0.558571185441022
    0.55853412965655
    0.558497040527712
    0.558459918812333
    0.558422765268242
    0.558385580653265
    0.558348365725227
    0.558311121241955
    0.558273847961277
    0.558236546641018
    0.558199218039005
    0.558161862913064
    0.558124482021022
    0.558087076120705
    0.558049645969941
    0.558012192326555
    0.557974715948373
    0.557937217593223
    0.557899698018931
    0.557862157983324
    0.557824598244227
    0.557787019559468
    0.557749422686873
    0.557711808384268
    0.55767417740948
    0.557636530520336
    0.557598868474662
    0.557561192030284
    0.557523501945028
    0.557485798976723
    0.557448083883193
    0.557410357422266
    0.557372620351768
    0.557334873429525
    0.557297117413364
    0.557259353061112
    0.557221581130595
    0.557183802379639
    0.557146017566071
    0.557108227447717
    0.557070432782405
    0.55703263432796
    0.556994832842209
    0.556957029082978
    0.556919223808095
    0.556881417775385];

mssThreshNeg = [mssThreshNeg(1:min(500,nTP)); mssThreshNeg(500)*ones(max(0,nTP-500),1)];
mssThreshPos = [mssThreshPos(1:min(500,nTP)); mssThreshPos(500)*ones(max(0,nTP-500),1)];

%%%%%%%%%%%%%%%%%%%%%%%

function [mssThreshNeg,mssThreshPos] = threshMSS2D_p05(nTP)

%2D, alpha = 0.05

%NEW

mssThreshNeg = [...
    NaN
    NaN
    NaN
    NaN
    NaN
    NaN
    NaN
    NaN
    NaN
    NaN
    NaN
    NaN
    NaN
    NaN
    NaN
    NaN
    NaN
    NaN
    NaN
    0.181502318197934
    0.18838154839312
    0.195093150236193
    0.201442121195694
    0.207306808998423
    0.212671284114308
    0.217532557256877
    0.221922548689565
    0.225932950092577
    0.229699916111634
    0.233254869995021
    0.236552670861111
    0.239577858589862
    0.242411883979074
    0.245071420325144
    0.24752759823054
    0.249788728574198
    0.251918804873189
    0.253897354787951
    0.255649897546368
    0.257179156137925
    0.258614575528493
    0.260031836782165
    0.261452589867952
    0.262895955219871
    0.264411374671013
    0.265954267057598
    0.267422783984008
    0.268774834751154
    0.27004432706212
    0.271212742725752
    0.27222691255208
    0.273083282167612
    0.273856763073577
    0.274571478179446
    0.275235026799138
    0.275888246341044
    0.276625829445396
    0.277455885316465
    0.278327934150515
    0.279214852070928
    0.280141531850763
    0.281070372256589
    0.28194568792854
    0.282754176677657
    0.283552543232945
    0.284351822876984
    0.285134951932275
    0.285910522477779
    0.286733405071101
    0.287603074387924
    0.28847846377164
    0.289344518753971
    0.290221884752678
    0.291086096888505
    0.291899333693971
    0.292654765103181
    0.293395567856581
    0.294121038994455
    0.294796425477663
    0.295404928327597
    0.295973031216056
    0.296495243035509
    0.296962554136571
    0.297386091414761
    0.297806251800542
    0.298221785133029
    0.298609627952852
    0.298968848762741
    0.299324796946564
    0.299656009252775
    0.299920405055219
    0.300113936894473
    0.300289897045359
    0.300465642712047
    0.300626347554218
    0.300772734121281
    0.300946839689796
    0.301161105705468
    0.301407626432376
    0.301699293564754
    0.302090771300807
    0.302597620306922
    0.303188096932729
    0.303828706002585
    0.304503647117815
    0.305156782347571
    0.305724734586455
    0.306185632342498
    0.30657911598172
    0.306929385133769
    0.307234447011637
    0.30750320909756
    0.307774567678586
    0.308045026162675
    0.308293662450211
    0.308519494296455
    0.308747270361001
    0.308971869842679
    0.309177841480646
    0.309397032188751
    0.309716334669726
    0.310182234779652
    0.310790790403923
    0.311508527402952
    0.312288439399305
    0.313084740117098
    0.313868425757027
    0.314632283201258
    0.31537604800628
    0.316106106881946
    0.316835567719803
    0.317576740361338
    0.31832627614131
    0.319068027983544
    0.319790445331619
    0.320491848838877
    0.321175898051449
    0.321847797499791
    0.322511776960654
    0.323167535303506
    0.32381425690187
    0.32444734520928
    0.3250569837267
    0.325635209870627
    0.326184914065035
    0.326713069256028
    0.327222064465545
    0.32771389695669
    0.328194976703951
    0.328674256987075
    0.329158154752876
    0.329647656229768
    0.330142502064348
    0.33064243290321
    0.331147189392948
    0.331656512180156
    0.332170141911429
    0.332687819233362
    0.333209284792549
    0.333734279235584
    0.334262543209063
    0.334793817359579
    0.335327842333727
    0.335864358778102
    0.336403107339297
    0.336943828663909
    0.33748626339853
    0.338030152189756
    0.338575235684181
    0.3391212545284
    0.339667949369007
    0.340215060852596
    0.340762329625763
    0.341309496335101
    0.341856301627205
    0.34240248614867
    0.34294779054609
    0.34349195546606
    0.344034721555173
    0.344575829460026
    0.345115019827211
    0.345652033303324
    0.346186610534959
    0.346718492168711
    0.347247418851174
    0.347773131228943
    0.348295369948611
    0.348813875656775
    0.349328389000027
    0.349838650624963
    0.350344401178178
    0.350845381306264
    0.351341331655818
    0.351831992873434
    0.352317105605706
    0.352796410499229
    0.353269648200596
    0.353736559356404
    0.354196884613245
    0.354650364617715
    0.355096803733521
    0.355536261192821
    0.355968859944884
    0.356394722938982
    0.356813973124385
    0.357226733450362
    0.357633126866185
    0.358033276321124
    0.358427304764448
    0.35881533514543
    0.359197490413338
    0.359573893517444
    0.359944667407017
    0.360309935031328
    0.360669819339648
    0.361024443281247
    0.361373929805394
    0.361718401861362
    0.362057982398419
    0.362392794365836
    0.362722960712884
    0.363048604388834
    0.363369848342954
    0.363686815524517
    0.363999628882791
    0.364308411367048
    0.364613285926558
    0.364914375510592
    0.365211803068419
    0.36550569154931
    0.365796163902535
    0.366083343077365
    0.366367352023071
    0.366648313688922
    0.366926351024189
    0.367201586978142
    0.367474144500052
    0.367744146539188
    0.368011716044823
    0.368276975966225
    0.368540049252665
    0.368801058853414
    0.369060127717742
    0.369317378794919
    0.369572935034215
    0.369826919384902
    0.370079454796249
    0.370330664217527
    0.370580670598006
    0.370829596886957
    0.371077545990059
    0.371324540638632
    0.371570583520403
    0.371815677323101
    0.372059824734454
    0.37230302844219
    0.372545291134039
    0.372786615497727
    0.373027004220985
    0.373266459991538
    0.373504985497117
    0.37374258342545
    0.373979256464264
    0.374215007301289
    0.374449838624251
    0.374683753120881
    0.374916753478905
    0.375148842386053
    0.375380022530053
    0.375610296598632
    0.37583966727952
    0.376068137260445
    0.376295709229134
    0.376522385873317
    0.376748169880722
    0.376973063939076
    0.377197070736108
    0.377420192959548
    0.377642433297122
    0.377863794436559
    0.378084279065588
    0.378303889871937
    0.378522629543333
    0.378740500767507
    0.378957506232185
    0.379173648625096
    0.379388930633969
    0.379603354946531
    0.379816924250512
    0.380029641233639
    0.380241508583641
    0.380452528988246
    0.380662705135182
    0.380872039712178
    0.381080535406962
    0.381288194907262
    0.381495020900807
    0.381701016075325
    0.381906183118544
    0.382110524718193
    0.382314041511462
    0.382516725933393
    0.38271856836849
    0.382919559201255
    0.383119688816193
    0.383318947597807
    0.3835173259306
    0.383714814199078
    0.383911402787743
    0.384107082081098
    0.384301842463648
    0.384495674319897
    0.384688568034347
    0.384880513991502
    0.385071502575867
    0.385261524171945
    0.385450569164239
    0.385638627937254
    0.385825690875492
    0.386011748363457
    0.386196790785654
    0.386380808526586
    0.386563791970756
    0.386745731502668
    0.386926617506826
    0.387106440367734
    0.387285190469894
    0.387462858197811
    0.387639433935989
    0.387814908068931
    0.38798927098114
    0.388162513057121
    0.388334624681377
    0.388505596238412
    0.388675418112729
    0.388844080688832
    0.389011574351225
    0.389177889484411
    0.389343016472895
    0.389506945701178
    0.389669667553767
    0.389831172415163
    0.389991450669871
    0.390150492702394
    0.390308288897236
    0.390464829638901
    0.390620105311892
    0.390774106300713
    0.390926822989868
    0.39107824576386
    0.391228371833239
    0.391377225712739
    0.391524838743139
    0.391671242265219
    0.391816467619758
    0.391960546147538
    0.392103509189336
    0.392245388085933
    0.392386214178109
    0.392526018806644
    0.392664833312316
    0.392802689035906
    0.392939617318193
    0.393075649499958
    0.393210816921979
    0.393345150925037
    0.393478682849911
    0.393611444037382
    0.393743465828227
    0.393874779563229
    0.394005416583165
    0.394135408228816
    0.394264785840962
    0.394393580760382
    0.394521824327855
    0.394649547884163
    0.394776782770084
    0.394903560326398
    0.395029911893884
    0.395155868813323
    0.395281462425494
    0.395406724071177
    0.395531685091152
    0.395656376826198
    0.395780830617095
    0.395905077804623
    0.396029149729561
    0.396153077732689
    0.396276893154787
    0.396400627336635
    0.396524311619012
    0.396647977342697
    0.396771655848472
    0.396895378477115
    0.397019176569406
    0.397143081466124
    0.397267124508051
    0.397391337035964
    0.397515750390644
    0.397640395912871
    0.397765295916105
    0.397890436604526
    0.398015795154997
    0.39814134874438
    0.398267074549535
    0.398392949747324
    0.39851895151461
    0.398645057028253
    0.398771243465115
    0.398897488002058
    0.399023767815944
    0.399150060083633
    0.399276341981989
    0.399402590687871
    0.399528783378143
    0.399654897229665
    0.3997809094193
    0.399906797123908
    0.400032537520351
    0.400158107785491
    0.40028348509619
    0.40040864662931
    0.400533569561711
    0.400658231070255
    0.400782608331805
    0.400906678523221
    0.401030418821365
    0.4011538064031
    0.401276818445286
    0.401399432124785
    0.401521624618459
    0.401643373103169
    0.401764654755777
    0.401885446753144
    0.402005726272133
    0.402125470489605
    0.402244656582421
    0.402363261727443
    0.402481263101533
    0.402598637881552
    0.402715363244362
    0.402831416366824
    0.402946774425801
    0.403061414598153
    0.403175314060743
    0.403288449990431
    0.40340079956408
    0.403512339958552
    0.403623048350707
    0.403732901917407
    0.403841884564866
    0.403950007116704
    0.404057287125893
    0.404163742145404
    0.40426938972821
    0.404374247427282
    0.404478332795592
    0.404581663386113
    0.404684256751815
    0.404786130445671
    0.404887302020653
    0.404987789029732
    0.40508760902588
    0.405186779562069
    0.405285318191272
    0.40538324246646
    0.405480569940604
    0.405577318166677
    0.40567350469765
    0.405769147086496
    0.405864262886186
    0.405958869649692
    0.406052984929986
    0.40614662628004
    0.406239811252825
    0.406332557401314
    0.406424882278478
    0.40651680343729
    0.40660833843072
    0.406699504811742
    0.406790320133326
    0.406880801948445
    0.406970967810071
    0.407060835271175
    0.407150421884729
    0.407239745203705
    0.407328822781075
    0.407417672169811
    0.407506310922885
    0.407594756593268
    0.407683026733933
    0.407771138897851
    0.407859110637994
    0.407946959507334
    0.408034703058843
    0.408122358845493
    0.408209944420255
    0.408297477336101
    0.408384975146004
    0.408472455402935];

mssThreshPos = [...
    NaN
    NaN
    NaN
    NaN
    NaN
    NaN
    NaN
    NaN
    NaN
    NaN
    NaN
    NaN
    NaN
    NaN
    NaN
    NaN
    NaN
    NaN
    NaN
    0.688317643178081
    0.684726982642884
    0.68109891780178
    0.677425966056172
    0.673772923654433
    0.670288565640082
    0.667028329047343
    0.663971485209905
    0.661098796538868
    0.658423236864972
    0.655915783776832
    0.653542912686538
    0.651304724479748
    0.649254881220427
    0.647402722076526
    0.64571892369454
    0.644173861928544
    0.642761044119332
    0.641442764872315
    0.640171172719201
    0.638926424664184
    0.637724234901367
    0.636561063652285
    0.635423767720294
    0.634316560623178
    0.633269160418724
    0.632282509308338
    0.631332603316206
    0.630403937638469
    0.629494519244589
    0.628598325433644
    0.627717703921463
    0.626877206866974
    0.626112399589143
    0.625433382597822
    0.624830074186726
    0.624272663307602
    0.623732990124725
    0.623188654303903
    0.622621882518039
    0.622018409224427
    0.621382154346502
    0.620734159692872
    0.620101273033383
    0.619513606228945
    0.619001390533913
    0.618570172161691
    0.618207471742387
    0.61789543013713
    0.617620008083951
    0.617368819561361
    0.617134673177448
    0.616912162432114
    0.61669548554837
    0.616482633313128
    0.61627408525523
    0.616069495218099
    0.61585972662781
    0.615645374006389
    0.615428467168467
    0.615204001060492
    0.614950113717697
    0.614648007493866
    0.614283787241613
    0.613854907780228
    0.613376822565723
    0.61288079011633
    0.612404917991244
    0.611979710511832
    0.611612841377266
    0.611296763087735
    0.611014387509404
    0.610750623745419
    0.610498967002745
    0.610260704299795
    0.610041184941574
    0.609841195252448
    0.609652559162008
    0.60946936691501
    0.609291220413442
    0.609125401270486
    0.60898090580048
    0.60886056117208
    0.608754670781189
    0.608646662655495
    0.608519212313296
    0.608359580758184
    0.608164909156993
    0.607938762999321
    0.607689273624184
    0.607429855604158
    0.607169871033368
    0.606918580257458
    0.606682763857289
    0.606465265393582
    0.606264217674601
    0.606072732783467
    0.605883775313072
    0.605691920712663
    0.605495671515588
    0.605301720282764
    0.605118361905657
    0.604953104735767
    0.604811833591361
    0.604699271334609
    0.604610098710317
    0.604531919619297
    0.604449657153941
    0.604349371760523
    0.604224907001748
    0.604072736775639
    0.603891192281183
    0.603684263288198
    0.603460619613105
    0.603231458616728
    0.60300509886874
    0.60279081690687
    0.60259239013901
    0.602412497491601
    0.602254141198669
    0.602116774525276
    0.601994020237082
    0.601873106297988
    0.601744314891167
    0.601600986274804
    0.601438352273111
    0.601256835807114
    0.6010618214076
    0.600857272726772
    0.600648222598096
    0.600441127912726
    0.600237228641353
    0.600036435488577
    0.599838696104781
    0.599643958140351
    0.59945216924567
    0.599263277071123
    0.599077229267095
    0.598893973483968
    0.598713457372128
    0.598535628581959
    0.598360434763844
    0.598187823568169
    0.598017742645318
    0.597850139645674
    0.597684962219622
    0.597522158017546
    0.597361674689831
    0.597203459886861
    0.59704746125902
    0.596893626456692
    0.596741903130262
    0.596592238930113
    0.59644458150663
    0.596298878510197
    0.5961550775912
    0.59601312640002
    0.595872972587044
    0.595734563802655
    0.595597847697237
    0.595462771921176
    0.595329284124854
    0.595197331958656
    0.595066863072967
    0.594937825118171
    0.594810165744652
    0.594683832602793
    0.594558773342981
    0.594434935615598
    0.594312267071029
    0.594190715359658
    0.59407022813187
    0.593950753038048
    0.593832237728578
    0.593714629853843
    0.593597877064227
    0.593481927010115
    0.59336672734189
    0.593252225709939
    0.593138369764643
    0.593025107156388
    0.592912392027208
    0.592800204485731
    0.592688531132239
    0.592577358567011
    0.592466673390325
    0.592356462202461
    0.5922467116037
    0.592137408194319
    0.5920285385746
    0.591920089344821
    0.591812047105262
    0.591704398456202
    0.591597129997921
    0.591490228330699
    0.591383680054814
    0.591277471770546
    0.591171590078175
    0.591066021577981
    0.590960752870242
    0.590855770555238
    0.590751061233249
    0.590646611504555
    0.590542407969434
    0.590438437228166
    0.590334685881031
    0.590231140528308
    0.590127787770276
    0.590024614207216
    0.589921606439406
    0.589818751067127
    0.589716034690657
    0.589613443910276
    0.589510965326264
    0.5894085855389
    0.589306291148464
    0.589204068755234
    0.589101904959491
    0.588999786361514
    0.588897699561582
    0.588795631159976
    0.588693567756974
    0.588591495952856
    0.588489402347901
    0.588387273542389
    0.5882850961366
    0.588182856730813
    0.588080541925307
    0.587978138320362
    0.587875632516257
    0.587773011113272
    0.587670265223054
    0.587567404002717
    0.587464441120744
    0.587361390245618
    0.58725826504582
    0.587155079189832
    0.587051846346137
    0.586948580183218
    0.586845294369556
    0.586742002573634
    0.586638718463934
    0.586535455708938
    0.586432227977129
    0.586329048936988
    0.586225932256999
    0.586122891605642
    0.586019940651402
    0.585917093062759
    0.585814362508197
    0.585711762656197
    0.585609307175242
    0.585507009733813
    0.585404884000394
    0.585302943643466
    0.585201202331512
    0.585099673733014
    0.584998371516454
    0.584897309350315
    0.584796500903079
    0.584695959843227
    0.584595699839243
    0.584495734559609
    0.584396077672806
    0.584296742847318
    0.584197743751626
    0.584099094054212
    0.58400080742356
    0.583902897528151
    0.583805378036467
    0.583708262616991
    0.583611564938205
    0.583515298668591
    0.583419477476631
    0.583324115030809
    0.583229224999606
    0.583134821051503
    0.583040916854985
    0.582947526078532
    0.582854662390628
    0.582762339459753
    0.582670569344341
    0.582579357662615
    0.582488708422753
    0.582398625632927
    0.582309113301314
    0.582220175436088
    0.582131816045425
    0.582044039137498
    0.581956848720483
    0.581870248802555
    0.581784243391889
    0.581698836496659
    0.581614032125041
    0.58152983428521
    0.58144624698534
    0.581363274233606
    0.581280920038183
    0.581199188407247
    0.581118083348971
    0.581037608871532
    0.580957768983103
    0.58087856769186
    0.580800009005978
    0.580722096933631
    0.580644835482995
    0.580568228662244
    0.580492280479554
    0.580416994943099
    0.580342376061053
    0.580268427841593
    0.580195154292893
    0.580122559423128
    0.580050647240473
    0.579979421753102
    0.579908886969191
    0.579839046896914
    0.579769905544447
    0.579701466919964
    0.57963373503164
    0.57956671388765
    0.57950040749617
    0.579434819865373
    0.579369955003435
    0.579305816918531
    0.579242409618836
    0.579179737112525
    0.579117803407771
    0.579056612512752
    0.57899616843564
    0.578936475184612
    0.57887753226413
    0.578819321163809
    0.57876181886955
    0.578705002367256
    0.57864884864283
    0.578593334682173
    0.578538437471188
    0.578484133995778
    0.578430401241845
    0.57837721619529
    0.578324555842016
    0.578272397167927
    0.578220717158923
    0.578169492800907
    0.578118701079782
    0.57806831898145
    0.578018323491813
    0.577968691596773
    0.577919400282233
    0.577870426534095
    0.577821747338262
    0.577773339680635
    0.577725180547118
    0.577677246923611
    0.577629515796019
    0.577581964150242
    0.577534568972184
    0.577487307247746
    0.577440155962831
    0.577393092103341
    0.577346092655179
    0.577299134604247
    0.577252194936447
    0.577205250637681
    0.577158278693852
    0.577111256090863
    0.577064159814614
    0.57701696685101
    0.576969654185951
    0.576922198805341
    0.576874577695081
    0.576826767841075
    0.576778746229223
    0.576730489845429
    0.576681975675596
    0.576633180705624
    0.576584081921417
    0.576534656308877
    0.576484880853906
    0.576434732542406
    0.576384193652116
    0.576333267628115
    0.57628196320732
    0.576230289126646
    0.57617825412301
    0.576125866933326
    0.576073136294511
    0.576020070943481
    0.575966679617151
    0.575912971052437
    0.575858953986256
    0.575804637155522
    0.575750029297153
    0.575695139148063
    0.575639975445168
    0.575584546925385
    0.575528862325629
    0.575472930382816
    0.575416759833862
    0.575360359415682
    0.575303737865193
    0.57524690391931
    0.57518986631495
    0.575132633789027
    0.575075215078458
    0.575017618920159
    0.574959854051046
    0.574901929208033
    0.574843853128038
    0.574785634547976
    0.574727282204763
    0.574668804835314
    0.574610211176546
    0.574551509965375
    0.574492709938715
    0.574433819833484
    0.574374848386596
    0.574315804334968
    0.574256696415516
    0.574197533365154
    0.574138323920801
    0.57407907681937
    0.574019800797778
    0.57396050459294
    0.573901196941774
    0.573841886581193
    0.573782582248115
    0.573723292679455
    0.573664026612129
    0.573604792783052
    0.573545598336405
    0.573486444045419
    0.573427329090593
    0.573368252652422
    0.573309213911404
    0.573250212048034
    0.573191246242811
    0.573132315676229
    0.573073419528787
    0.57301455698098
    0.572955727213306
    0.572896929406261
    0.572838162740342
    0.572779426396045
    0.572720719553867
    0.572662041394306
    0.572603391097857
    0.572544767845017
    0.572486170816283
    0.572427599192152
    0.57236905215312
    0.572310528879684
    0.572252028552341
    0.572193550351587
    0.572135093457919
    0.572076657051834
    0.572018240313829
    0.5719598424244
    0.571901462564043
    0.571843099913256
    0.571784753652535
    0.571726422962378
    0.571668107023279
    0.571609805015737
    0.571551516120248
    0.571493239517309
    0.571434974387416
    0.571376719911066
    0.571318475268756
    0.571260239640982
    0.571202012208242
    0.571143792151031
    0.571085578649847
    0.571027370885185
    0.570969168037544
    0.570910969287419
    0.570852773815308
    0.570794580801706
    0.570736389427111
    0.57067819887202];

mssThreshNeg = [mssThreshNeg(1:min(500,nTP)); mssThreshNeg(500)*ones(max(0,nTP-500),1)];
mssThreshPos = [mssThreshPos(1:min(500,nTP)); mssThreshPos(500)*ones(max(0,nTP-500),1)];

%%%%%%%%%%%%%%%%%%%%%%%

function [mssThreshNeg,mssThreshPos] = threshMSS2D_p01(nTP)

%2D, alpha = 0.01

mssThreshNeg = [...
    NaN
    NaN
    NaN
    NaN
    NaN
    NaN
    NaN
    NaN
    NaN
    NaN
    NaN
    NaN
    NaN
    NaN
    NaN
    NaN
    NaN
    NaN
    NaN
    0.0915486471971317
    0.100671616999023
    0.10955875894705
    0.118011316449105
    0.125887144108242
    0.133131289897661
    0.139737005238879
    0.145705995027657
    0.151112359984369
    0.156116370113472
    0.160762844472145
    0.165040671537703
    0.168986102045849
    0.172745327012372
    0.176396211132299
    0.179890461183287
    0.183179581516945
    0.186324128001645
    0.189312219341217
    0.192057512725438
    0.194543351409659
    0.19686401647993
    0.199027837443714
    0.201007741228426
    0.202815435121145
    0.204569899703765
    0.206285392283955
    0.207922811901324
    0.209488237078343
    0.21106571111798
    0.212656039210603
    0.214201355442208
    0.215664673494099
    0.217100369103699
    0.218518129108796
    0.219842456946586
    0.221042500405771
    0.222194807337113
    0.223323108866386
    0.224392374718687
    0.225378598561683
    0.226312088482021
    0.227196258649204
    0.228021290134248
    0.228795724532316
    0.22956637008591
    0.230333334413044
    0.231051575162275
    0.231718271951513
    0.232411606349154
    0.233167094538452
    0.233969907270292
    0.23479588018697
    0.235654650167993
    0.23650797621557
    0.237283968701925
    0.237939777825461
    0.238514955122998
    0.239019646431625
    0.239455360292773
    0.239859569531789
    0.240339076079567
    0.240962143803522
    0.241736325218023
    0.242646005137539
    0.243673798661359
    0.244705490801636
    0.245582657313175
    0.246225636624915
    0.246662840085691
    0.246937802051194
    0.247087468664853
    0.247185626269538
    0.247345981475725
    0.247611175374529
    0.247967147763268
    0.248420470864584
    0.249039545486146
    0.249836171851885
    0.250761993534409
    0.251753953378734
    0.252764538789389
    0.253707819505046
    0.254510239535313
    0.255160322959598
    0.255686609596412
    0.256097900337804
    0.25638543292716
    0.256570478785386
    0.256721788766842
    0.256879441756437
    0.257065639645503
    0.257316094144063
    0.257670554948797
    0.258118921083491
    0.258619293978429
    0.25913597887017
    0.25966798909132
    0.260202425460088
    0.260726204395079
    0.261244478486087
    0.261789413423779
    0.262370256151468
    0.26298020787078
    0.263605592409436
    0.264243344259435
    0.264898937926137
    0.265578041020037
    0.266287786325582
    0.267045020853011
    0.267867082118982
    0.268758506631976
    0.269717180148353
    0.270718313478107
    0.271719370284138
    0.272685180432856
    0.273603896855532
    0.274477698940946
    0.275311298454889
    0.276115784214232
    0.276894716268988
    0.277653187680488
    0.278404128778922
    0.279161090188406
    0.279924497201986
    0.280700062659582
    0.281506239802355
    0.282349766064853
    0.283227440466056
    0.284124261058993
    0.285024372076602
    0.285912714842098
    0.286786328624336
    0.28764547314468
    0.288490408124494
    0.289321393285141
    0.290138688347987
    0.290942553034394
    0.291733247065727
    0.292511030163349
    0.293276162048625
    0.294028902442917
    0.294769511067591
    0.29549824764401
    0.296215371893538
    0.296921143537539
    0.297615822297377
    0.298299667894415
    0.298972940050017
    0.299635898485549
    0.300288802922372
    0.300931913081852
    0.301565488685352
    0.302189789454236
    0.302805075109868
    0.303411605373612
    0.304009639966831
    0.304599438610891
    0.305181261027153
    0.305755366936984
    0.306322016061745
    0.306881468122802
    0.307433982841518
    0.307979819939257
    0.308519239137383
    0.30905250015726
    0.309579862720252
    0.310101586547722
    0.310617931361035
    0.311129156881555
    0.311635522830644
    0.312137288929669
    0.312634714899991
    0.313128060462975
    0.313617585339985
    0.314103549252385
    0.314586211921539
    0.315065833068811
    0.315542672415564
    0.316016989683162
    0.31648904459297
    0.316959049690909
    0.317427028821132
    0.317892958652347
    0.318356815853266
    0.318818577092599
    0.319278219039055
    0.319735718361345
    0.320191051728178
    0.320644195808265
    0.321095127270316
    0.321543822783041
    0.321990259015149
    0.322434412635351
    0.322876260312358
    0.323315778714878
    0.323752944511623
    0.324187734371301
    0.324620124962624
    0.325050092954301
    0.325477615015042
    0.325902667813558
    0.326325228018558
    0.326745272298753
    0.327162777322852
    0.327577719759565
    0.327990076277603
    0.328399823545676
    0.328806938232494
    0.329211397006766
    0.329613176537203
    0.330012253492516
    0.330408604541413
    0.330802206352604
    0.331193035594802
    0.331581068936714
    0.331966283047051
    0.332348654594523
    0.332728160247841
    0.333104776675714
    0.333478480546853
    0.333849248529966
    0.334217057293766
    0.334581883506961
    0.334943703838261
    0.335302494956377
    0.335658233530019
    0.336010896227896
    0.336360459718719
    0.336706900671199
    0.337050195754043
    0.337390333782881
    0.337727352161005
    0.338061300438624
    0.33839222816595
    0.338720184893191
    0.339045220170558
    0.339367383548262
    0.339686724576511
    0.340003292805516
    0.340317137785487
    0.340628309066634
    0.340936856199167
    0.341242828733296
    0.34154627621923
    0.341847248207181
    0.342145794247358
    0.34244196388997
    0.342735806685229
    0.343027372183343
    0.343316709934523
    0.343603869488979
    0.343888900396921
    0.344171852208559
    0.344452774474103
    0.344731716743762
    0.345008728567748
    0.345283859496269
    0.345557159079537
    0.34582867686776
    0.346098462411149
    0.346366565259914
    0.346633034964265
    0.346897921074411
    0.347161273140564
    0.347423140712932
    0.347683573341726
    0.347942620577156
    0.348200331969432
    0.348456757068764
    0.348711945425362
    0.348965946589435
    0.349218810111194
    0.34947058554085
    0.349721322428611
    0.349971070324687
    0.35021987877929
    0.350467797342629
    0.350714875564913
    0.350961162996353
    0.351206709187159
    0.351451552286378
    0.351695684838409
    0.351939087986488
    0.352181742873849
    0.352423630643729
    0.352664732439364
    0.35290502940399
    0.353144502680842
    0.353383133413156
    0.353620902744167
    0.353857791817112
    0.354093781775227
    0.354328853761746
    0.354562988919907
    0.354796168392944
    0.355028373324093
    0.355259584856591
    0.355489784133673
    0.355718952298574
    0.355947070494531
    0.356174119864779
    0.356400081552555
    0.356624936701093
    0.35684866645363
    0.357071251953401
    0.357292674343643
    0.35751291476759
    0.357731954368479
    0.357949774289546
    0.358166355674025
    0.358381679665154
    0.358595727406168
    0.358808480040302
    0.359019918710793
    0.359230024560876
    0.359438778733787
    0.359646162372761
    0.359852156621035
    0.360056742621844
    0.360259901518424
    0.360461614454011
    0.360661862571841
    0.360860627015149
    0.361057888927171
    0.361253629451142
    0.3614478297303
    0.361640470907879
    0.361831534127115
    0.362021000531244
    0.362208851263502
    0.362395077878643
    0.362579713577497
    0.362762801972413
    0.362944386675739
    0.363124511299823
    0.363303219457014
    0.363480554759661
    0.363656560820112
    0.363831281250715
    0.364004759663819
    0.364177039671773
    0.364348164886925
    0.364518178921623
    0.364687125388216
    0.364855047899052
    0.36502199006648
    0.365187995502849
    0.365353107820507
    0.365517370631802
    0.365680827549083
    0.365843522184698
    0.366005498150996
    0.366166799060326
    0.366327468525035
    0.366487550157473
    0.366647087569988
    0.366806124374928
    0.366964704184642
    0.367122870611478
    0.367280667267785
    0.367438137765911
    0.367595325718205
    0.367752274737016
    0.367909028434691
    0.368065630423579
    0.36822212431603
    0.36837855372439
    0.368534962261009
    0.368691393538236
    0.368847891168418
    0.369004498763905
    0.369161259937044
    0.369318218300185
    0.369475417465675
    0.369632901045863
    0.369790712653099
    0.369948895899729
    0.370107494398103
    0.370266551760569
    0.370426111599476
    0.370586204357528
    0.370746807798854
    0.370907886517936
    0.37106940510926
    0.371231328167309
    0.371393620286568
    0.37155624606152
    0.37171917008665
    0.371882356956442
    0.372045771265379
    0.372209377607945
    0.372373140578626
    0.372537024771905
    0.372700994782265
    0.372865015204191
    0.373029050632168
    0.373193065660678
    0.373357024884207
    0.373520892897237
    0.373684634294254
    0.373848213669741
    0.374011595618183
    0.374174744734063
    0.374337625611866
    0.374500202846074
    0.374662441031174
    0.374824304761648
    0.374985758631981
    0.375146767236657
    0.375307295170159
    0.375467307026972
    0.375626767401581
    0.375785640888468
    0.375943892082118
    0.376101485577015
    0.376258385967643
    0.376414557848487
    0.37656996581403
    0.376724574458756
    0.376878348377149
    0.377031252163694
    0.377183250412875
    0.377334307719175
    0.377484388677078
    0.377633457881069
    0.377781479925632
    0.37792841940525
    0.378074240914409
    0.378218909047591
    0.378362388399281
    0.378504653545506
    0.378645718988463
    0.378785609211892
    0.378924348699534
    0.379061961935128
    0.379198473402415
    0.379333907585135
    0.379468288967028
    0.379601642031835
    0.379733991263295
    0.379865361145148
    0.379995776161135
    0.380125260794996
    0.380253839530471
    0.3803815368513
    0.380508377241224
    0.380634385183982
    0.380759585163315
    0.380884001662963
    0.381007659166666
    0.381130582158164
    0.381252795121197
    0.381374322539506
    0.381495188896831
    0.381615418676911
    0.381735036363487
    0.3818540664403
    0.381972533391089
    0.382090461699594
    0.382207875849556
    0.382324800324715
    0.38244125960881
    0.382557278185583
    0.382672880538773
    0.38278809115212
    0.382902934509366
    0.383017435094248
    0.383131617390509
    0.383245505881888
    0.383359125052125
    0.38347249938496
    0.383585653364134
    0.383698611473387
    0.383811398196458
    0.383924038017089
    0.384036555419019
    0.384148974885988
    0.384261320901737
    0.384373617950005
    0.384485890514533];

mssThreshPos = [...
    NaN
    NaN
    NaN
    NaN
    NaN
    NaN
    NaN
    NaN
    NaN
    NaN
    NaN
    NaN
    NaN
    NaN
    NaN
    NaN
    NaN
    NaN
    NaN
    0.750265572468563
    0.746217200245581
    0.742117813924218
    0.737971145291753
    0.733868862858166
    0.729990679197887
    0.72642797020255
    0.723151820426775
    0.720100058054654
    0.717243686891916
    0.714539753569261
    0.711942629982225
    0.709437924910794
    0.707063448960662
    0.704827892167264
    0.702742930416691
    0.700809934456745
    0.699002230389337
    0.697251475783319
    0.695506848752802
    0.693755881259079
    0.692005833213999
    0.690250199400663
    0.688488829435655
    0.686749928971535
    0.685094841765023
    0.683543128639221
    0.682076119929991
    0.68068969839937
    0.67939232180094
    0.678169225133095
    0.676997114155101
    0.675870542406679
    0.674798084022674
    0.673787244132564
    0.672840773668072
    0.6719629173696
    0.6711712197284
    0.670481788364631
    0.669888635448306
    0.669352281682291
    0.668821849278578
    0.668272966170585
    0.667694158042765
    0.667080362161112
    0.666443037485745
    0.665806294582891
    0.665191658158698
    0.664605261628074
    0.664038809345026
    0.663475049814771
    0.662891295378162
    0.66227266327622
    0.661602462278384
    0.660878944212457
    0.660119830056248
    0.659361808983804
    0.658654097553155
    0.658023632321308
    0.657479430340614
    0.657022613472294
    0.656639767181658
    0.656311041066085
    0.656003958270266
    0.655679899899229
    0.655305809459256
    0.654863493266331
    0.654356042971839
    0.653810104486025
    0.653268568000375
    0.652770343902119
    0.652330127806941
    0.651938949223061
    0.651572547174639
    0.651209378257629
    0.650831713846563
    0.650431564273119
    0.650015229474996
    0.649592760916608
    0.649176997882014
    0.648775290746419
    0.648385441168244
    0.647996575812798
    0.647606452025784
    0.647227027965201
    0.646880436291042
    0.646585088989155
    0.646345926231941
    0.646159401222635
    0.646013206093368
    0.645893667166536
    0.645787332011113
    0.64568377692405
    0.645585544725221
    0.645512402866599
    0.645481670944532
    0.645497610094702
    0.645551740336929
    0.645642775059394
    0.645756876297718
    0.645865950633395
    0.645926981080937
    0.64591322017193
    0.645833324540451
    0.6457212646252
    0.64560709332276
    0.645503369754937
    0.645403119626725
    0.645289488284613
    0.645151492243763
    0.644985882199931
    0.64479473590996
    0.64458363298583
    0.644346693193633
    0.644070425128451
    0.643748116315153
    0.643382120158113
    0.642982469896317
    0.642566714319327
    0.642152722007253
    0.641756605306366
    0.641392947963751
    0.641076536621273
    0.640818315156262
    0.64061838927512
    0.640466558943228
    0.640345736458033
    0.640235104749114
    0.640121271598142
    0.63999725911874
    0.639860286735131
    0.639714002494725
    0.639559716993427
    0.639397711794269
    0.63922826846028
    0.63905166855449
    0.638868193639929
    0.638678125279626
    0.638481745036612
    0.638279334473917
    0.63807117515457
    0.637857548641601
    0.63763873649804
    0.637415020286916
    0.637186681571261
    0.636954001914103
    0.636717262878473
    0.636476746027399
    0.636232732923914
    0.635985505131045
    0.635735344211823
    0.635482531729277
    0.635227349246438
    0.634970078326336
    0.634711000532
    0.63445039742646
    0.634188550572747
    0.633925741533889
    0.633662251872917
    0.63339836315286
    0.633134356936749
    0.632870514787613
    0.632607118268483
    0.632344448942388
    0.632082788372357
    0.631822418121421
    0.63156361975261
    0.631306674828953
    0.631051864913481
    0.630799471569223
    0.630549776359209
    0.630303060846469
    0.630059606594032
    0.62981969516493
    0.629583608122191
    0.629351627028845
    0.629124033447922
    0.628901108942453
    0.628683135075466
    0.628470393409992
    0.628263165509061
    0.628061661345803
    0.627865804533751
    0.62767544709654
    0.627490441057802
    0.627310638441172
    0.627135891270282
    0.626966051568769
    0.626800971360264
    0.626640502668401
    0.626484497516815
    0.626332807929139
    0.626185285929007
    0.626041783540053
    0.62590215278591
    0.625766245690212
    0.625633914276593
    0.625505010568687
    0.625379386590128
    0.625256894364548
    0.625137385915583
    0.625020713266865
    0.624906728442029
    0.624795283464708
    0.624686230358536
    0.624579421147147
    0.624474707854175
    0.624371942503253
    0.624270977118015
    0.624171663722095
    0.624073854339126
    0.623977400992743
    0.623882155706579
    0.623787970504268
    0.623694697409444
    0.62360218844574
    0.62351029563679
    0.623418871006228
    0.623327766577688
    0.623236834374803
    0.623145926421208
    0.623054894740535
    0.622963591356419
    0.622871868292494
    0.622779577572393
    0.622686571219751
    0.6225927012582
    0.622497819711374
    0.622401778602908
    0.622304429956435
    0.622205625795589
    0.622105251936466
    0.622003329365013
    0.62189991285964
    0.621795057198756
    0.621688817160771
    0.621581247524094
    0.621472403067135
    0.621362338568304
    0.621251108806009
    0.621138768558661
    0.621025372604668
    0.620910975722441
    0.62079563269039
    0.620679398286922
    0.620562327290449
    0.620444474479379
    0.620325894632123
    0.620206642527089
    0.620086772942687
    0.619966340657327
    0.619845400449419
    0.619724007097371
    0.619602215379593
    0.619480080074496
    0.619357655960488
    0.619234997815978
    0.619112160419378
    0.618989198549095
    0.618866166983539
    0.618743120501121
    0.61862011388025
    0.618497201899334
    0.618374439336784
    0.61825188097101
    0.61812958158042
    0.618007595943424
    0.617885978838433
    0.617764785043854
    0.617644069338099
    0.617523886499575
    0.617404291306694
    0.617285338537864
    0.617167082971495
    0.617049579385997
    0.616932882559779
    0.61681704727125
    0.616702128298821
    0.6165881804209
    0.616475258415897
    0.616363417062222
    0.616252695681619
    0.616143071769166
    0.616034507363278
    0.61592696450237
    0.615820405224854
    0.615714791569145
    0.615610085573656
    0.615506249276802
    0.615403244716996
    0.615301033932652
    0.615199578962184
    0.615098841844006
    0.614998784616532
    0.614899369318175
    0.61480055798735
    0.61470231266247
    0.61460459538195
    0.614507368184203
    0.614410593107643
    0.614314232190684
    0.614218247471739
    0.614122600989223
    0.61402725478155
    0.613932170887133
    0.613837311344387
    0.613742638191724
    0.61364811346756
    0.613553699210308
    0.613459357458381
    0.613365050250194
    0.613270739624161
    0.613176387618695
    0.613081956272211
    0.612987407623121
    0.612892703709841
    0.612797806570784
    0.612702678244363
    0.612607280768993
    0.612511576183088
    0.612415526525061
    0.612319093833326
    0.612222240146298
    0.61212492750239
    0.612027117940016
    0.611928773497589
    0.611829856213525
    0.611730328126235
    0.611630151274136
    0.611529287695639
    0.61142769942916
    0.611325359209217
    0.611222282554745
    0.611118495680787
    0.611014024802383
    0.610908896134574
    0.610803135892401
    0.610696770290905
    0.610589825545126
    0.610482327870107
    0.610374303480889
    0.610265778592511
    0.610156779420015
    0.610047332178442
    0.609937463082833
    0.60982719834823
    0.609716564189672
    0.609605586822201
    0.609494292460858
    0.609382707320685
    0.609270857616721
    0.609158769564008
    0.609046469377587
    0.6089339832725
    0.608821337463786
    0.608708558166487
    0.608595671595644
    0.608482703966298
    0.60836968149349
    0.608256630392261
    0.608143576877652
    0.608030547164704
    0.607917567468457
    0.607804664003954
    0.607691862986234
    0.60757919063034
    0.607466673151311
    0.607354336764189
    0.607242207684016
    0.607130312125831
    0.607018676304676
    0.606907326435591
    0.606796288733619
    0.606685589413799
    0.606575254691174
    0.606465310780783
    0.606355783897668
    0.60624670025687
    0.60613808607343
    0.606029967562388
    0.605922370938787
    0.605815318742873
    0.605708818815719
    0.605602875323605
    0.60549749243281
    0.605392674309615
    0.605288425120299
    0.605184749031141
    0.605081650208421
    0.604979132818419
    0.604877201027414
    0.604775859001685
    0.604675110907513
    0.604574960911177
    0.604475413178957
    0.604376471877132
    0.604278141171981
    0.604180425229785
    0.604083328216823
    0.603986854299374
    0.603891007643718
    0.603795792416136
    0.603701212782905
    0.603607272910306
    0.603513976964619
    0.603421329112123
    0.603329333519098
    0.603237994351823
    0.603147315776578
    0.603057301959643
    0.602967957067296
    0.602879285265819
    0.602791290721489
    0.602703977600587
    0.602617350069393
    0.602531412294186
    0.602446168441246
    0.602361622676852
    0.602277779167284
    0.602194642078821
    0.602112215577743
    0.60203050383033
    0.601949511002861
    0.601869241261616
    0.601789698772874
    0.601710887702916
    0.60163281221802
    0.601555476484466
    0.601478884668534
    0.601403040936504
    0.601327949454654
    0.601253611173497
    0.601180014180471
    0.601107143347248
    0.601034983545497
    0.600963519646891
    0.600892736523098
    0.600822619045791
    0.600753152086641
    0.600684320517317
    0.600616109209491
    0.600548503034833
    0.600481486865014
    0.600415045571706
    0.600349164026578
    0.600283827101303
    0.600219019667549
    0.600154726596989
    0.600090932761293
    0.600027623032131
    0.599964782281175
    0.599902395380095
    0.599840447200563
    0.599778922614249
    0.599717806492823
    0.599657083707957
    0.599596739131321
    0.599536757634586
    0.599477124089423
    0.599417823367502
    0.599358840340496
    0.599300159880073
    0.599241766857905
    0.599183646145664
    0.599125782615018
    0.59906816113764
    0.599010766585201
    0.59895358382937
    0.598896597741819
    0.598839793194219
    0.59878315505824
    0.598726668205553
    0.598670317507829
    0.598614087836739
    0.598557964063953
    0.598501931061143
    0.598445973699978
    0.59839007685213
    0.59833422538927
    0.598278404183069
    0.598222598105197];

mssThreshNeg = [mssThreshNeg(1:min(500,nTP)); mssThreshNeg(500)*ones(max(0,nTP-500),1)];
mssThreshPos = [mssThreshPos(1:min(500,nTP)); mssThreshPos(500)*ones(max(0,nTP-500),1)];

%%%%%%%%%%%%%%%%%%%%%%%

function [mssThreshNeg,mssThreshPos] = threshMSS3D_p20(nTP)

%3D, alpha = 0.2

%threshold curve parameters
turnPointsM = [20 60 200 500];
turnPointsP = [20 60 100 500];
slopeM = [0.001359327197525 0.000280732353905 0.000108640351462 0];
slopeP = [-0.000746906632969 -0.000348357342053 -0.000037074457479 0];
interseptM = [0.299405762325546 0.364121452942741 0.39853985343145 0.452860029162242];
interseptP = [0.609183824504497 0.585270867049586 0.554142578592106 0.535605349852791];

%threshold curve evaluation
[mssThreshNeg,mssThreshPos] = getThreshCurve(turnPointsM,turnPointsP,...
    slopeM,slopeP,interseptM,interseptP,nTP);

%%%%%%%%%%%%%%%%%%%%%%%

function [mssThreshNeg,mssThreshPos] = threshMSS3D_p10(nTP)

%3D, alpha = 0.1

%threshold curve parameters
turnPointsM = [20 60 200 500];
turnPointsP = [20 50 100 500];
slopeM = [0.001587591010958 0.000334860802492 0.000129492806027 0];
slopeP = [-0.001033885285809 -0.000409579718459 -0.000060967307118 0];
interseptM = [0.258642983727984 0.333806796235989 0.374880395529001 0.439626798542322];
interseptP = [0.644618927677761 0.613403649310276 0.57854240817611 0.548058754617302];

%threshold curve evaluation
[mssThreshNeg,mssThreshPos] = getThreshCurve(turnPointsM,turnPointsP,...
    slopeM,slopeP,interseptM,interseptP,nTP);

%%%%%%%%%%%%%%%%%%%%%%%

function [mssThreshNeg,mssThreshPos] = threshMSS3D_p05(nTP)

%3D, alpha = 0.05

%threshold curve parameters
turnPointsM = [20 60 200 500];
turnPointsP = [20 50 80 500];
slopeM = [0.001897773535643 0.000390310143108 0.00014504631569 0];
slopeP = [-0.001281266825658 -0.000815336722161 -0.000079326637385 0];
interseptM = [0.217433363211081 0.307881166763189 0.356933932246803 0.429457090091877];
interseptP = [0.683369058696754 0.660072553521921 0.601191746739831 0.56152842804728];

% threshold curve evaluation

%threshold curve evaluation
[mssThreshNeg,mssThreshPos] = getThreshCurve(turnPointsM,turnPointsP,...
    slopeM,slopeP,interseptM,interseptP,nTP);

function [mssThreshNeg,mssThreshPos] = threshMSS3D_p01(nTP)

%3D, alpha = 0.01

%threshold curve parameters
turnPointsM = [20 50 200 500];
turnPointsP = [20 90 200 500];
slopeM = [0.00330864748993062 0.000607571513122873 0.000183704012624133 0];
slopeP = [-0.0011001388082068 -0.000401685416953439 -7.31736162754304e-05 0];
interseptM = [0.100955106735165 0.236008905575552 0.3207824056753 0.412634411987367];
interseptP = [0.747201850658008 0.684341045445206 0.618638685309604 0.582051877171889];

%threshold curve evaluation
[mssThreshNeg,mssThreshPos] = getThreshCurve(turnPointsM,turnPointsP,...
    slopeM,slopeP,interseptM,interseptP,nTP);

%%%%%%%%%%%%%%%%%%%%%%%

%% threshold curve evaluation subfunction

function [mssThreshNeg,mssThreshPos] = getThreshCurve(turnPointsM,...
    turnPointsP,slopeM,slopeP,interseptM,interseptP,nTP)

%confined diffusion threshold
mssThreshNeg = NaN(1,turnPointsM(1)-1);
for i = 1 : length(turnPointsM)-1
    x = turnPointsM(i) : turnPointsM(i+1)-1;
    mssThreshNeg = [mssThreshNeg slopeM(i)*x+interseptM(i)]; %#ok<AGROW>
end
x = turnPointsM(end) : nTP;
mssThreshNeg = [mssThreshNeg slopeM(end)*x+interseptM(end)];

%directed diffusion threshold
mssThreshPos = NaN(1,turnPointsP(1)-1);
for i = 1 : length(turnPointsP)-1
    x = turnPointsP(i) : turnPointsP(i+1)-1;
    mssThreshPos = [mssThreshPos slopeP(i)*x+interseptP(i)]; %#ok<AGROW>
end
x = turnPointsP(end) : nTP;
mssThreshPos = [mssThreshPos slopeP(end)*x+interseptP(end)];

