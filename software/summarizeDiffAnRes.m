function [probMotionType,motionChar,errFlag] = summarizeDiffAnRes(tracks,...
    minTrackLen,probDim,diffAnalysisRes,extractType)
%SUMMARIZEDIFFANRES calculates motion type probabilities and motion characteristics from diffusion analysis
%
%SYNOPSIS [probMotionType,motionChar,errFlag] = summarizeDiffAnRes(tracks,...
%    minTrackLen,probDim,diffAnalysisRes,extractType)
%
%INPUT  tracks     : Output of trackCloseGapsKalman.
%       minTrackLen: Minimum length of a track to be used in getting
%                    motion type statistics.
%                    Optional. Default: 5.
%       probDim    : Dimensionality - 2 for 2D, 3 for 3D.
%                    Optional. Default: 2.
%       diffAnalysisRes: Diffusion analysis results (output of
%                    trackDiffAnalysis1). Optional. If not input, it will
%                    be calculated.
%       extractType: 1 - Analyze every track segment separately.
%                    2 - Extract from each compound track the longest
%                        trajectory to use in analysis - NOT IMPLEMENTED
%                        YET.
%                    Must use same extractType as in trackDiffusionAnalysis1.
%                    Variable irrelevant if tracks are input as a matrix.
%                    Optional. Default: 1.
%
%OUTPUT probMotionType: 11-by-3 array. Rows refer to:
%                    (1) Linear & 1D immobile,
%                    (2) Linear & 1D confined,
%                    (3) Linear & 1D Brownian,
%                    (4) Linear & 1D directed,
%                    (5) Linear & diffusion undetermined,
%                    (6) Not linear/unclassified & 2D immobile,
%                    (7) Not linear/unclassied & 2D confined,
%                    (8) Not linear/unclassified & 2D Brownian,
%                    (9) Not linear/unclassified & 2D directed,
%                   (10) Not linear & diffusion undetermined, and
%                   (11) Completely undetermined.
%                    1st column: Each motion type's probability.
%                    2nd column: Probability of a motion type within its
%                    category, i.e. 1-5 are in the linear category, while
%                    6-10 or 6-9&11 are in the not-linear category (6-10 are
%                    used when asymmetry is checked for, 6-9&11 are used
%                    when asymmetry is not checked for).
%                    3rd column: Number of features per frame in each
%                    motion category, same rows as first column.
%       motionChar : Structure summarizing motion characteristics.
%                    Contains the fields "linear" and "notLinear". 
%                       The field "linear" contains the sub-field "all"
%                    (i.e. all trajectories classified as linear are lumped
%                    together).
%                       The field "notLinear" contains the sub-fields
%                    "immobile","confined", "brownian" and "directed",
%                    showing the characteristics of each of these motion
%                    categories.
%                       For each of these sub-fields, there are two
%                    sub-sub-fields: 
%           .distribution: Nx4 array. Each row belongs to a separate
%                          trajectory in the category. The columns
%                          store the diffusion coefficient, the confinement
%                          radius, the other confinement dimension in the
%                          case of linear trajectories, and trajectory
%                          lifetime.
%           .meanStd     : 2x4 array. Columns are same as "distribution".
%                          1st row is the mean, 2nd row is the standard
%                          deviation.
%           .median      : 1x4 array showing the median of the columns
%                          in "distribution".
%
%Khuloud Jaqaman, March 2010
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

%% output

probMotionType = [];
motionChar = [];
errFlag = 0;

%% input

if nargin < 1 || isempty(tracks)
    disp('summarizeDiffAnRes: Missing input argument!');
    return
end

if nargin < 2 || isempty(minTrackLen)
    minTrackLen = 5;
end

if nargin < 3 || isempty(probDim)
    probDim = 2;
end

if nargin < 4 || isempty(diffAnalysisRes)
    [diffAnalysisRes,errFlag] = trackDiffusionAnalysis1(tracks,1,probDim,...
        1,[0.05 0.05],0);
    if errFlag
        return
    end
end

if nargin < 5 || isempty(extractType)
    extractType = 1;
else
    if ~any(extractType == [1 2])
        disp('--trackDiffusionAnalysis1: Variable extractType should be 1 or 2.');
        errFlag = 1;
    end
end

%% preamble

%keep only tracks with length >= minTrackLen
criteria.lifeTime.min = minTrackLen;
indx = chooseTracks(tracks,criteria);
clear criteria
if ~isempty(indx)
    tracks = tracks(indx);
    diffAnalysisRes = diffAnalysisRes(indx);
else
    disp('Ignoring minTrackLen because imposing it leaves no tracks to analyze.')
end

%get number of tracks and number of frames
numTracks = length(tracks);
seqOfEvents = vertcat(tracks.seqOfEvents);
numFrames = max(seqOfEvents(:,1));

%put tracks in matrix format
if extractType == 1
    [tracksMat,tracksIndxMat,trackStartRow] = convStruct2MatIgnoreMS(tracks);
end

%get number of track segments
numTrackSegments = size(tracksMat,1);

%get track lengths
trackSegmentLft = getTrackSEL(tracksMat);
trackSegmentLft = trackSegmentLft(:,3);

%% features

%get average number of features per frame
numFeatTot = length(find(tracksIndxMat(:)));
aveFeatPerFrame = numFeatTot / numFrames;

%% motion type probabilities

%get track segment classification from diffusion analysis
trackSegmentClass = vertcat(diffAnalysisRes.classification);

%get indices of the different motion types
indxLinImm      = find( trackSegmentClass(:,1) == 1   & trackSegmentClass(:,3) == 0   );
indxLinConf     = find( trackSegmentClass(:,1) == 1   & trackSegmentClass(:,3) == 1   );
indxLinBrown    = find( trackSegmentClass(:,1) == 1   & trackSegmentClass(:,3) == 2   );
indxLinDir      = find( trackSegmentClass(:,1) == 1   & trackSegmentClass(:,3) == 3   );
indxLinUndet    = find( trackSegmentClass(:,1) == 1   & isnan(trackSegmentClass(:,3)) );
indxNonlinImm   = find( trackSegmentClass(:,1) ~= 1   & trackSegmentClass(:,2) == 0   );
indxNonlinConf  = find( trackSegmentClass(:,1) ~= 1   & trackSegmentClass(:,2) == 1   );
indxNonlinBrown = find( trackSegmentClass(:,1) ~= 1   & trackSegmentClass(:,2) == 2   );
indxNonlinDir   = find( trackSegmentClass(:,1) ~= 1   & trackSegmentClass(:,2) == 3   );
indxNonlinUndet = find( trackSegmentClass(:,1) == 0   & isnan(trackSegmentClass(:,2)) );
indxUndetUndet  = find( isnan(trackSegmentClass(:,1)) & isnan(trackSegmentClass(:,2)) );

%calculate number of track segments per motion type
numSegmentsType = [length(indxLinImm) length(indxLinConf) ...
    length(indxLinBrown) length(indxLinDir) length(indxLinUndet) ...
    length(indxNonlinImm) length(indxNonlinConf) length(indxNonlinBrown) ...
    length(indxNonlinDir) length(indxNonlinUndet) length(indxUndetUndet)]';

%calculate fraction of track segments falling in each motion type
fracSegmentsType = numSegmentsType / numTrackSegments;

%calculate number of features in each motion type
numFeatType = [length(find(tracksIndxMat(indxLinImm,:))) ...
    length(find(tracksIndxMat(indxLinConf,:))) ...
    length(find(tracksIndxMat(indxLinBrown,:))) ...
    length(find(tracksIndxMat(indxLinDir,:))) ...
    length(find(tracksIndxMat(indxLinUndet,:))) ...
    length(find(tracksIndxMat(indxNonlinImm,:))) ...
    length(find(tracksIndxMat(indxNonlinConf,:))) ...
    length(find(tracksIndxMat(indxNonlinBrown,:))) ...
    length(find(tracksIndxMat(indxNonlinDir,:))) ...
    length(find(tracksIndxMat(indxNonlinUndet,:))) ...
    length(find(tracksIndxMat(indxUndetUndet,:)))]';

%get fraction of features undergoing each motion type - this is the
%probability of a feature undergoing a certain motion type
probFeatType = numFeatType / numFeatTot;

%within the overall linear and non-linear motion types, calculate the
%probabilities of the different sub-types
probFeatLinSubType = probFeatType(1:5) / sum(probFeatType(1:5));
if all(isnan(trackSegmentClass(:,1)))
    probFeatNonlinSubType = probFeatType([6:9 11]) / sum(probFeatType([6:9 11]));
    probFeatSubType = [probFeatLinSubType; probFeatNonlinSubType(1:4); NaN; ...
        probFeatNonlinSubType(5)];
else
    probFeatNonlinSubType = probFeatType(6:10) / sum(probFeatType(6:10));
    probFeatSubType = [probFeatLinSubType; probFeatNonlinSubType; NaN];
end

%combine the motion type and sub-type probabilities
probMotionType = [probFeatType probFeatSubType numFeatType/numFrames];

%% motion characteristics

%extract diffusion coefficients and confinement radii
diffCoefAll = catStruct(1,'diffAnalysisRes.fullDim.genDiffCoef(:,3)');
confRadAll = catStruct(1,'diffAnalysisRes.confRadInfo.confRadius');

%distribute motion characteristics based on motion types

%linear, all together
motionCharTmp = [diffCoefAll([indxLinConf;indxLinBrown;indxLinDir]) ...
    confRadAll([indxLinConf;indxLinBrown;indxLinDir],:) ...
    trackSegmentLft([indxLinConf;indxLinBrown;indxLinDir])];
motionChar.linear.all.distribution = motionCharTmp;
motionChar.linear.all.meanStd = [nanmean(motionCharTmp,1); nanstd(motionCharTmp,[],1)];
motionChar.linear.all.median = nanmedian(motionCharTmp,1);

%non-linear & immobile
motionCharTmp = [diffCoefAll(indxNonlinImm) confRadAll(indxNonlinImm,:) ...
    trackSegmentLft(indxNonlinImm)];
motionChar.notLinear.immobile.distribution = motionCharTmp;
motionChar.notLinear.immobile.meanStd = [nanmean(motionCharTmp,1); nanstd(motionCharTmp,[],1)];
motionChar.notLinear.immobile.median = nanmedian(motionCharTmp,1);

%non-linear & confined
motionCharTmp = [diffCoefAll(indxNonlinConf) confRadAll(indxNonlinConf,:) ...
    trackSegmentLft(indxNonlinConf)];
motionChar.notLinear.confined.distribution = motionCharTmp;
motionChar.notLinear.confined.meanStd = [nanmean(motionCharTmp,1); nanstd(motionCharTmp,[],1)];
motionChar.notLinear.confined.median = nanmedian(motionCharTmp,1);

%non-linear & Brownian
motionCharTmp = [diffCoefAll(indxNonlinBrown) confRadAll(indxNonlinBrown,:) ...
    trackSegmentLft(indxNonlinBrown)];
motionChar.notLinear.brownian.distribution = motionCharTmp;
motionChar.notLinear.brownian.meanStd = [nanmean(motionCharTmp,1); nanstd(motionCharTmp,[],1)];
motionChar.notLinear.brownian.median = nanmedian(motionCharTmp,1);

%non-linear & directed
motionCharTmp = [diffCoefAll(indxNonlinDir) confRadAll(indxNonlinDir,:) ...
    trackSegmentLft(indxNonlinDir)];
motionChar.notLinear.directed.distribution = motionCharTmp;
motionChar.notLinear.directed.meanStd = [nanmean(motionCharTmp,1); nanstd(motionCharTmp,[],1)];
motionChar.notLinear.directed.median = nanmedian(motionCharTmp,1);

%% ~~~ the end ~~~

