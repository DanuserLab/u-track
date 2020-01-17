function [diffAnalysisRes,errFlag] = trackDiffusionAnalysis1(tracks,...
    extractType,probDim,checkAsym,alphaValues,plotRes,confRadMin)
%TRACKDIFFUSIONANALYSIS performs diffusion analysis, checking first for asymmetric tracks
%
%SYNOPSIS [diffAnalysisRes,errFlag] = trackDiffusionAnalysis1(tracks,...
%    extractType,probDim,checkAsym,alphaValues,plotRes,confRadMin)
%
%INPUT  tracks      : -- EITHER --
%                     Output of trackWithGapClosing (matrix),
%                     -- OR --
%                     Output of trackCloseGapsKalman (structure, possibly
%                     with merges/splits.
%       extractType : 1 - Analyze every track segment separately.
%                     2 - NOTTHING IMPLEMENTED AT THE MOMENT.
%                     Variable irrelevant if tracks are input as a matrix.
%                     Optional. Default: 1.
%       probDim     : Problem dimensionality. Optional. Default: 2.
%       checkAsym   : 1 to check for asymmetric tracks and to analyze their
%                     diffusion after dimensionality reduction, 0
%                     otherwise. Optional. Default: 0.
%       alphaValues : Row vector with 2 entries. 
%                     *** First entry is the alpha-value for MSS analysis.
%                     Can take the values 0.2, 0.1, 0.05 or 0.01; or the
%                     negative of these values for NEW CLASSIFICATION SCHEME.
%                     See trackMSSAnalysis for most up-to-date values
%                     allowed and explanation of NEW CLASSIFICATION SCHEME.
%                     *** Second entry is the alpha-value for asymmetry
%                     determination. Can take the values 0.2, 0.1, 0.05 or
%                     0.01. See help of asymDeterm2D3D for most up-to-date
%                     values allowed.
%                     *** Optional. Default: [0.05 0.1]. If only one value is
%                     entered, it is taken as the alpha-value for MSS
%                     analysis, while the alpha-value for asymmetry
%                     analysis is given the default value.
%
%       plotRes     : 1 to plot results, 0 otherwise. Optional. Default: 0.
%                     Results can be plotted only if problem is 2D.
%                     See plotTracksDiffAnalysis2D for color code.
%       confRadMin  : 1 to calculate the confinement radius of confined
%                     particles using the minimum positional standard
%                     deviation; OR
%                     0 to calculate it using the mean positional standard
%                     deviation; OR
%                     2 to approximate the confinement area by a rectangle
%                     and calculate the length of both edges.
%                     Optional. Default: 0.
%
%OUTPUT diffAnalysisRes : Structure array with the following fields per
%                         track:
%           .classification: Number of segment x 3 matrix. 
%                           *Column 1: Classification based on asymmetry.
%                            1 = asymmetric, 0 = not asymmetric.
%                           *Column 2: Classification based on moment
%                            scaling spectrum analysis applied to the
%                            tracks using their full dimensionality.
%                            0 = immobile, 1 = confined Brownian, 
%                            2 = pure Brownian, 3 = directed motion.
%                           *Column 3: Classification of motion along
%                            the preferred direction for linear tracks,
%                            also based on moment scaling spectrum analysis.
%                            0 = Immobile, 1 = confined Brownian,
%                            2 = pure Brownian, 3 = directed motion.
%           .fullDim       : MSS analysis results for full dimensionality.
%                            Structure with fields:
%               .mssSlope    : Slope of the line representing moment 
%                              scaling power vs. moment order.
%               .genDiffCoef : Generalized diffusion coefficient for each
%                              order employed. The "normal" (MSD) diffusion
%                              coefficient is the 3rd entry.
%               .scalingPower: The moment scaling power for each order
%                              employed. The scaling power of the MSD is
%                              the 3rd entry.
%           .oneDim        : MSS analysis results for reduced dimensionality.
%                            Structure with same fields as fullDim.
%           .confRadius    : For particles undergoing confined motion, the
%                            confinement radius assuming circular
%                            confinement, or the 2 confinement dimensions
%                            assuming rectangular confinement.
%                            For particles undergoing linear motion, the 2
%                            confinement dimensions parallel and
%                            perpendicular to the direction of motion.
%           .summary       : Summary of diffusion characteristics for all
%                            tracks combined, as output by
%                            summarizeDiffAnRes. Field filled only in first
%                            array entry, with all others empty, since
%                            there is no need for repetition.
%       errFlag         : 0 if executed normally, 1 otherwise.
%
%REMARKS 
%(1)While tracks do not have to be linear in order to be asymmetric,
%the last analysis step assumes that tracks are linear.
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

%Khuloud Jaqaman, March 2008

%% Output

diffAnalysisRes = [];
errFlag = 0;

%% Input

%check whether tracks were input
if nargin < 1
    disp('--trackDiffusionAnalysis1: Please input at least the tracks to be analyzed!');
    errFlag = 1;
    return
end

if nargin < 2 || isempty(extractType)
    extractType = 1;
else
    if ~any(extractType == [1 2])
        disp('--trackDiffusionAnalysis1: Variable extractType should be 1 or 2.');
        errFlag = 1;
    end
end

if nargin < 3 || isempty(probDim)
    probDim = 2;
end

if nargin < 4 || isempty(checkAsym)
    checkAsym = 0;
end

if nargin < 5 || isempty(alphaValues)
    alphaValues = [0.05 0.1];
elseif length(alphaValues) == 1
    alphaValues = [alphaValues 0.1];
end

if nargin < 6 || isempty(plotRes)
    plotRes = 0;
elseif plotRes == 1 && probDim ~= 2
    disp('--trackDiffusionAnalysis1: Cannot plot tracks if problem is not 2D!');
    plotRes = 0;
end

if nargin < 7 || isempty(confRadMin)
    confRadMin = 0;
end

if errFlag
    disp('--trackDiffusionAnalysis1: Please fix input variables');
    return
end

%% track extraction for analysis

%store input tracks in a new variable
tracksInput = tracks;

%extract segments for analysis if tracks were input as a structure that
%might contain merges and splits
%the point is to reduce compound tracks that contain merges and splits into
%simple separate tracks
%thus this step is not necessary if the tracks were input as a matrix,
%which by definition does not contain unresolved compound tracks.
if isstruct(tracks)

    %get number of input tracks from structure
    numInputTracks = length(tracksInput);

    clear tracks

    switch extractType

        case 1 %retrieve every track segment separately

            [tracks,dummy,compTrackStartRow,numSegments] = ...
                convStruct2MatIgnoreMS(tracksInput,1);
            
        case 2 %some other type

            disp('Sorry - not implemented yet!')
            errFlag = 1;
            return

    end

else

    %get number of input tracks from matrix
    numInputTracks = size(tracksInput,1);

    %indicate rows where tracks start (trivial in this case)
    compTrackStartRow = (1 : numInputTracks)';

    %indicate number of segments in each track (1 for all tracks)
    numSegments = ones(numInputTracks,1);

end

%get number of track segments to be analyzed
numTrackSegments = size(tracks,1);

%% moment scaling spectrum analysis on full-dimensionality data

%this analysis is based on Ewers et al (PNAS 2005) & Ferrari et al (Physica
%D 2001)
%it is applied to the tracks using the real problem dimensionality

%call the moment scaling spectrum analysis code
momentOrders = 0 : 6;
[trackClassMSS,mssSlope,genDiffCoef,scalingPower,normDiffCoef] = ...
    trackMSSAnalysis(tracks,probDim,momentOrders,alphaValues(1));

%% track classification based on asymmetry

%this classification scheme is taken from Huet et al (BJ 2006)
%it classifies tracks as asymmetric or not, based on the scatter of
%positions along them

%reserve memory for results
trackClassAsym = NaN(numTrackSegments,1);
asymParam = NaN(numTrackSegments,1);
indxAsym = [];

if checkAsym

    %assign alpha-value for threshold determination
    alphaAsym = alphaValues(2);

    %find indices of tracks whose length >= 5 frames
    criteria.lifeTime.min = 5;
    indx4asymClass = chooseTracks(tracks,criteria);
    clear criteria

    %go over all of these tracks
    if ~isempty(indx4asymClass')
        for iTrack = indx4asymClass'

            %get the particle positions along the track
            coordX = tracks(iTrack,1:8:end)';
            coordY = tracks(iTrack,2:8:end)';
            coordZ = tracks(iTrack,3:8:end)';
            coordXYZ = [coordX coordY coordZ];

            %determine whether the track is sufficiently asymmetric
            [asymParamT,asymFlag] = asymDeterm2D3D(coordXYZ(:,1:probDim),alphaAsym);

            %classify track as ...
            %1 = linear, if the asymmetry parameter is larger than the threshold
            %0 = not linear, if the asymmetry parameter is smaller than the
            %threshold
            %otherwise, keep track classification as undetermined
            trackClassAsym(iTrack,:) = asymFlag;

            %also save asymmetry parameter
            asymParam(iTrack,:) = asymParamT;

        end
    end
    
    %find indices of all tracks classified as asymmetric
    indxAsym = find(trackClassAsym(:,1) == 1);
    
end

%% moment scaling spectrum analysis on reduced-dimensionality data

%this analysis is also based on Ewers et al (PNAS 2005) & Ferrari et al
%(Physica D 2001)
%it is applied to the tracks using their reduced dimensionality
%it assumes that tracks classified as asymmetric are in fact linear
%thus, it finds their preferred direction of motion, and then performs the
%analysis on the steps taken along that preferred direction

%reserve memory for results
trackClass1D = NaN(numTrackSegments,1);
mssSlope1D = NaN(numTrackSegments,1);
genDiffCoef1D = NaN(numTrackSegments,length(momentOrders));
scalingPower1D = NaN(numTrackSegments,length(momentOrders));
normDiffCoef1D = NaN(numTrackSegments,1);

if checkAsym && ~isempty(indxAsym)

    %project positions of linear tracks onto direction of motion
    numCol = size(tracks,2) / 8;
    tracksAsym = NaN(length(indxAsym),numCol*8);
    iAsym = 0;
    for iTrack = indxAsym'

        %get the positions in this track and their standard deviations
        trackCoordX = tracks(iTrack,1:8:end)';
        deltaCoordX = tracks(iTrack,5:8:end)';
        trackCoordY = tracks(iTrack,2:8:end)';
        deltaCoordY = tracks(iTrack,6:8:end)';
        trackCoordZ = tracks(iTrack,3:8:end)';
        deltaCoordZ = tracks(iTrack,7:8:end)';
        trackCoord = [trackCoordX trackCoordY trackCoordZ];
        deltaCoord = [deltaCoordX deltaCoordY deltaCoordZ];
        trackCoord = trackCoord(:,1:probDim);
        deltaCoord = deltaCoord(:,1:probDim);

        %project onto direction of motion
        [posAlongDir,deltaPosAlongDir] = projectCoordOntoDir(trackCoord,...
            deltaCoord,[],[]);

        %construct matrix of linear tracks with projected positions
        trackCoord2 = [posAlongDir zeros(numCol,3) deltaPosAlongDir zeros(numCol,3)]';
        trackCoord2 = trackCoord2(:)';
        iAsym = iAsym + 1;
        tracksAsym(iAsym,:) = trackCoord2;

    end

    %call the moment scaling spectrum analysis code
    momentOrders = 0 : 6;
    [trackClassT,mssSlopeT,genDiffCoefT,scalingPowerT,normDiffCoefT] = ...
        trackMSSAnalysis(tracksAsym,1,momentOrders,alphaValues(1));

    %since not all track segments are linear, put analysis results in their
    %proper place among all track segment
    trackClass1D(indxAsym) = trackClassT;
    mssSlope1D(indxAsym) = mssSlopeT;
    genDiffCoef1D(indxAsym,:) = genDiffCoefT;
    scalingPower1D(indxAsym,:) = scalingPowerT;
    normDiffCoef1D(indxAsym) = normDiffCoefT;
    
end

%% confinement radius estimation

%find all tracks classified as confined or immobile
indxConf = find( trackClassMSS <= 1 );

%remove from the list those tracks classified as linear
indxConf = setdiff(indxConf,indxAsym);

%reserve memory
confRadius = NaN(numTrackSegments,2);
trackCenter = NaN(numTrackSegments,probDim);
prefDir = NaN(numTrackSegments,probDim);

%estimate the confinement radius of confined tracks
if ~isempty(indxConf)
    for iTrack = indxConf'
        %     for iTrack = 1 : numTrackSegments

        %get track coordinates
        xCoord = (tracks(iTrack,1:8:end))';
        yCoord = (tracks(iTrack,2:8:end))';
        zCoord = (tracks(iTrack,3:8:end))';
        xyzCoord = [xCoord yCoord zCoord];
        
        %calculate the variance-covariance matrix of this track's positions
        xyzCov = nancov(xyzCoord(:,1:probDim));

        %find the eignevalues and eigenvectors of the variance-covariance
        %matrix
        eigenVal = eig(xyzCov);

        %calculate the track's confinement radius
        if confRadMin == 1
            confRadius(iTrack,1) = sqrt( min(eigenVal) * (probDim + 2) );
        else
            confRadius(iTrack,1) = sqrt( mean(eigenVal) * (probDim + 2) );
        end
        

    end
end

%add confined tracks to asymmetric tracks if requested
if confRadMin == 2
    indxAsym = [indxAsym; indxConf];
end

%estimate the confinement radii of asymmetric tracks
if ~isempty(indxAsym)
    for iTrack = indxAsym'

        %get track coordinates
        xCoord = (tracks(iTrack,1:8:end))';
        yCoord = (tracks(iTrack,2:8:end))';
        zCoord = (tracks(iTrack,3:8:end))';
        xyzCoord = [xCoord yCoord zCoord];

        %find the eignevalues of the variance-covariance matrix of this track's
        %positions
        [eigenVec,eigenVal] = eig(nancov(xyzCoord(:,1:probDim)));
        eigenVal = diag(eigenVal);

        %calculate the confinement radius along the preferred direction of
        %motion
        confRadius(iTrack,2) = sqrt( max(eigenVal) * (3) );

        %calculate the confinement radius perpendicular to the preferred
        %direction of motion
        confRadius(iTrack,1) = sqrt( mean(eigenVal(eigenVal~=max(eigenVal))) * (probDim + 1) );

        %calculate the track's center
        trackCenter(iTrack,:) = nanmean(xyzCoord(:,1:probDim));

        %store the preferred direction of motion
        prefDir(iTrack,:) = eigenVec(:,eigenVal==max(eigenVal))';

    end
end

%% save results in output structure

%reserve memory
diffAnalysisRes = repmat(struct('classification',[],'fullDim',[],'oneDim',...
    [],'confRadInfo',[],'summary',[]),numInputTracks,1);

%go over all input tracks
for iTrack = 1 : numInputTracks
    
    %store classification
    %column 1: asymmetry, column 2: MSS on full D, column 3: MSS on 1D
    diffAnalysisRes(iTrack).classification = [trackClassAsym(...
        compTrackStartRow(iTrack):compTrackStartRow(iTrack)+...
        numSegments(iTrack)-1) trackClassMSS(compTrackStartRow(iTrack):...
        compTrackStartRow(iTrack)+numSegments(iTrack)-1) ...
        trackClass1D(compTrackStartRow(iTrack):compTrackStartRow(iTrack)...
        +numSegments(iTrack)-1)];
    
    %store parameters of full D classification
    fullDim.mssSlope = mssSlope(compTrackStartRow(iTrack):...
        compTrackStartRow(iTrack)+numSegments(iTrack)-1);
    fullDim.genDiffCoef = genDiffCoef(compTrackStartRow(iTrack):...
        compTrackStartRow(iTrack)+numSegments(iTrack)-1,:);
    fullDim.scalingPower = scalingPower(compTrackStartRow(iTrack):...
        compTrackStartRow(iTrack)+numSegments(iTrack)-1,:);
    fullDim.normDiffCoef = normDiffCoef(compTrackStartRow(iTrack):...
        compTrackStartRow(iTrack)+numSegments(iTrack)-1);
    diffAnalysisRes(iTrack).fullDim = fullDim;

    %store parameters of 1D classification
    oneDim.mssSlope = mssSlope1D(compTrackStartRow(iTrack):...
        compTrackStartRow(iTrack)+numSegments(iTrack)-1);
    oneDim.genDiffCoef = genDiffCoef1D(compTrackStartRow(iTrack):...
        compTrackStartRow(iTrack)+numSegments(iTrack)-1,:);
    oneDim.scalingPower = scalingPower1D(compTrackStartRow(iTrack):...
        compTrackStartRow(iTrack)+numSegments(iTrack)-1,:);
    oneDim.normDiffCoef = normDiffCoef1D(compTrackStartRow(iTrack):...
        compTrackStartRow(iTrack)+numSegments(iTrack)-1);
    diffAnalysisRes(iTrack).oneDim = oneDim;
    
    %store confinement radius information
    confRadInfo.confRadius = confRadius(compTrackStartRow(iTrack):...
        compTrackStartRow(iTrack)+numSegments(iTrack)-1,:);
    confRadInfo.trackCenter = trackCenter(compTrackStartRow(iTrack):...
        compTrackStartRow(iTrack)+numSegments(iTrack)-1,:);
    confRadInfo.prefDir = prefDir(compTrackStartRow(iTrack):...
        compTrackStartRow(iTrack)+numSegments(iTrack)-1,:);
    if diffAnalysisRes(iTrack).classification(3) ==3
        xCoord = (tracks(iTrack,1:8:end))';
        confRadInfo.driftSpeed = (2*confRadius(compTrackStartRow(iTrack)+numSegments(iTrack)-1,2))/length(xCoord);
    else
        confRadInfo.driftSpeed = NaN;
    end
    diffAnalysisRes(iTrack).confRadInfo = confRadInfo;

end

%call code to summarize diffusion analysis results
%store in .summary field of first track - rest stay empty
if isstruct(tracksInput)
    minTrackLen = 5;
    [probMotionType,motionChar,errFlag] = summarizeDiffAnRes(tracksInput,minTrackLen,probDim,diffAnalysisRes,extractType);
    diffAnalysisRes(1).summary.probMotionType = probMotionType;
    diffAnalysisRes(1).summary.motionChar = motionChar;
end

%% plotting

%plot results if requested
if plotRes
    plotTracksDiffAnalysis2D(tracksInput,diffAnalysisRes,[],1);
end

%% ~~~ the end ~~~
