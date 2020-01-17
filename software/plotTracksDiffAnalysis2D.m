function plotTracksDiffAnalysis2D(trackedFeatureInfo,diffAnalysisRes,timeRange,...
    newFigure,image,showConf,simplifyLin,offset,hideUnclass)
%PLOTTRACKSDIFFANALYSIS2D plots tracks in 2D highlighting the different diffusion types
%
%SYNOPSIS plotTracksDiffAnalysis2D(trackedFeatureInfo,diffAnalysisRes,timeRange,...
%    newFigure,image,showConf,simplifyLin,offset)
%
%INPUT  trackedFeatureInfo: -- EITHER -- 
%                           Output of trackWithGapClosing:
%                           Matrix indicating the positions and amplitudes 
%                           of the tracked features to be plotted. Number 
%                           of rows = number of tracks, while number of 
%                           columns = 8*number of time points. Each row 
%                           consists of 
%                           [x1 y1 z1 a1 dx1 dy1 dz1 da1 x2 y2 z2 a2 dx2 dy2 dz2 da2 ...]
%                           in image coordinate system (coordinates in
%                           pixels). NaN is used to indicate time points 
%                           where the track does not exist.
%                           -- OR -- 
%                           Output of trackCloseGapsKalman:
%                           Structure array with number of entries equal to
%                           the number of tracks (or compound tracks when
%                           merging/splitting are considered). Contains the
%                           fields:
%           .tracksCoordAmpCG: The positions and amplitudes of the tracked
%                              features, after gap closing. Number of rows
%                              = number of track segments in compound
%                              track. Number of columns = 8 * number of 
%                              frames the compound track spans. Each row
%                              consists of 
%                              [x1 y1 z1 a1 dx1 dy1 dz1 da1 x2 y2 z2 a2 dx2 dy2 dz2 da2 ...]
%                              NaN indicates frames where track segments do
%                              not exist.
%           .seqOfEvents     : Matrix with number of rows equal to number
%                              of events happening in a track and 4
%                              columns:
%                              1st: Frame where event happens;
%                              2nd: 1 - start of track, 2 - end of track;
%                              3rd: Index of track segment that ends or starts;
%                              4th: NaN - start is a birth and end is a death,
%                                   number - start is due to a split, end
%                                   is due to a merge, number is the index
%                                   of track segment for the merge/split.
%       diffAnalysisRes   : -- EITHER --
%                           Structure of diffusion analysis results as
%                           output by the code "trackDiffusionAnalysis1".
%                           -- OR --
%                           Structure of diffusion mode analysis results as
%                           output by the code "classifyTrackDiffMode".
%       timeRange         : 2-element row vector indicating time range to plot. 
%                           Optional. Default: whole movie.
%       newFigure         : 1 if plot should be made in a new figure
%                           window, 0 otherwise (in which case it will be
%                           plotted in an existing figure window).
%                           Optional. Default: 1.
%       image             : An image that the tracks will be overlaid on if
%                           newFigure=1. It will be ignored if newFigure=0.
%                           Optional. Default: no image.
%       showConf          : 1 to show confinement radii, 0 otherwise.
%                           Optional. Default: 0.
%       simplifyLin       : 1 to color all linear groups as one category;
%                           2 to include "random+superdiffusive" also in
%                           the overall linear category;
%                           0 to keep each group separate.
%                           Optional. Default: 0.
%       offset            : [dx,dy] that is to be added to the coordinates.
%                           Optional. Default: [0,0]
%
%       hideUnclass       : 1 to not display any tracks shorter than 20
%                           frames
%
%OUTPUT The plot.
%       Color coding:
%
%   IF OUTPUT OF trackDiffusionAnalysis1
%
%       linear & 1D immobile -> gray
%       linear & 1D confined diffusion -> orange
%       linear & 1D normal diffusion -> red
%       linear & 1D super diffusion -> green
%       linear & too short to analyze 1D diffusion -> yellow
%       random/unclassified & 2D immobile -> brown
%       random/unclassified & 2D confined diffusion -> blue
%       random/unclassified & 2D normal diffusion -> cyan
%       random/unclassified & 2D super diffusion -> magenta
%       random & too short to analyze 2D diffusion -> purple
%       too short for any analysis -> black, or white if overlaying on an
%       image
%
%       If simplifyLin = 1, all linear groups will be colored red
%       If simplifyLin = 2, all linear groups + random&super-diffusive will
%       be colored red
%
%   IF OUTPUT OF classifyTrackDiffMode
%
%       Mode 1 -> blue
%       Mode 2 -> cyan
%       Mode 3 -> red
%       Mode 4 -> green
%       Mode 5 -> magenta
%       too short to classify -> black, or white if overlaying on an image
%
%       If there are more than 5 modes, for now the code will complain and
%       exit
%
%Khuloud Jaqaman, March 2008
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

%% Input

%check whether correct number of input arguments was used
if nargin < 2
    disp('--plotTracksDiffAnalysis2D Incorrect number of input arguments!');
    return
end

%get number of tracks and number of time points
if isstruct(trackedFeatureInfo) %if tracks are in structure format
    numTracks = length(trackedFeatureInfo);
    tmp = vertcat(trackedFeatureInfo.seqOfEvents);
    numTimePoints = max(tmp(:,1));
    clear tmp
else %if tracks are in matrix format
    [numTracks,numTimePoints] = size(trackedFeatureInfo);
    numTimePoints = numTimePoints/8;
end

errFlag = 0;

%check whether a time range for plotting was input
if nargin < 3 || isempty(timeRange)
    timeRange = [1 numTimePoints];
else
    if timeRange(1) < 1 || timeRange(2) > numTimePoints
        disp('--plotTracksDiffAnalysis2D: Wrong time range for plotting!');
        errFlag = 1;
    end
end

%check whether newFigure was input
if nargin < 4 || isempty(newFigure)
    newFigure = 1;
else
    if newFigure ~= 0 && newFigure ~= 1
        disp('--plotTracksDiffAnalysis2D: newFigure should be 0 or 1!');
        errFlag = 1;
    end
end

%check whether user supplied an image
if nargin < 5 || isempty(image)
    image = [];
end

%check whether to plot confinement radii
if nargin < 6 || isempty(showConf)
    showConf = 0;
end

%check whether to simplify color-coding of linear tracks
if nargin < 7 || isempty(simplifyLin)
    simplifyLin = 0;
end

%check whether there is an offset
if nargin < 8 || isempty(offset)
    offset = [0 0];
end

%check whether there to include tracks below 20 frames
if nargin < 9 
    hideUnclass = 0;
end

%exit if there are problems in input variables
if errFlag
    disp('--plotTracksDiffAnalysis2D: Please fix input data!');
    return
end

%if input diffusion analysis is in fact diffusion mode analysis, disable
%showConf and simplifyLin
if ~isfield(diffAnalysisRes,'classification')
    showConf = 0;
    simplifyLin = 0;
end

%% Pre-processing

if isstruct(trackedFeatureInfo) %if tracks are input in structure format

    %store the input structure as a variable with a different name
    inputStructure = trackedFeatureInfo;
    clear trackedFeatureInfo;
    
    %get number of segments making each track
    numSegments = zeros(numTracks,1);
    for i = 1 : numTracks
        numSegments(i) = size(inputStructure(i).tracksCoordAmpCG,1);
    end

    %if all tracks have only one segment ...
    if max(numSegments) == 1

        %indicate that there are no compound tracks with merging and splitting branches
        mergeSplit = 0;

        %locate the row of the first track of each compound track in the
        %big matrix of all tracks (to be constructed in the next step) 
        %in this case of course every compound track is simply one track
        %without branches
        trackStartRow = (1:numTracks)';

        %store tracks in a matrix
        trackedFeatureInfo = NaN*ones(numTracks,8*numTimePoints);
        for i = 1 : numTracks
            startTime = inputStructure(i).seqOfEvents(1,1);
            endTime   = inputStructure(i).seqOfEvents(end,1);
            trackedFeatureInfo(i,8*(startTime-1)+1:8*endTime) = inputStructure(i).tracksCoordAmpCG;
        end
        
    else %if some tracks have merging/splitting branches
        
        %indicate that in the variable mergeSplit
        mergeSplit = 1;
        
        %locate the row of the first track of each compound track in the
        %big matrix of all tracks (to be constructed in the next step)
        trackStartRow = ones(numTracks,1);
        for iTrack = 2 : numTracks
            trackStartRow(iTrack) = trackStartRow(iTrack-1) + numSegments(iTrack-1);            
        end
        
        %put all tracks together in a matrix
        trackedFeatureInfo = NaN*ones(trackStartRow(end)+numSegments(end)-1,8*numTimePoints);
        for i = 1 : numTracks
            startTime = inputStructure(i).seqOfEvents(1,1);
            endTime   = inputStructure(i).seqOfEvents(end,1);
            trackedFeatureInfo(trackStartRow(i):trackStartRow(i)+...
                numSegments(i)-1,8*(startTime-1)+1:8*endTime) = ...
                inputStructure(i).tracksCoordAmpCG;
        end
        
    end    
    
else %if tracks are not input in structure format

    %indicate that there are no compound tracks with merging and splitting branches
    mergeSplit = 0;
    
    %locate the row of the first track of each compound track in the
    %big matrix of all tracks
    %in this case of course every compound track is simply one track
    %without branches
    trackStartRow = (1:numTracks)';

end

%get the x,y-coordinates of features in all tracks
tracksX = trackedFeatureInfo(:,1:8:end)' + offset(1);
tracksY = trackedFeatureInfo(:,2:8:end)' + offset(2);

%find x-coordinate limits
minXCoord = min(floor(min(tracksX(:))),0);
maxXCoord =  ceil(max(tracksX(:)));

%find y-coordinate limits
minYCoord = min(floor(min(tracksY(:))),0);
maxYCoord =  ceil(max(tracksY(:)));

%get number of track segment to be plotted
numTrackSegments = size(tracksX,2);

%% track segment color based on track type

if isfield(diffAnalysisRes,'classification') %output of trackDiffusionAnalysis1
    
    %get track segment types from diffusion analysis
    trackSegmentType = vertcat(diffAnalysisRes.classification);
    if hideUnclass==1
        checkClass = ~isnan(trackSegmentType(:,2));
        trackSegmentType = trackSegmentType(checkClass,:);
        numTrackSegments = size(trackSegmentType,1);
    end
    %color coding:
    %       linear & 1D immobile -> gray
    %       linear & 1D confined diffusion -> orange
    %       linear & 1D normal diffusion -> red
    %       linear & 1D super diffusion -> green
    %       linear & too short to analyze 1D diffusion -> yellow
    %       random/unclassified & 2D immobile -> brown
    %       random/unclassified & 2D confined diffusion -> blue
    %       random/unclassified & 2D normal diffusion -> cyan
    %       random/unclassified & 2D super diffusion -> magenta
    %       random & too short to analyze 2D diffusion -> purple
    %       too short for any analysis -> black, or white if overlaying on an
    %       image
    %
    %       If simplifyLin = 1, all linear groups will be colored red.
    %       If simplifyLin = 2, all linear groups + random&super-diffusive will
    %       be colored red.
    
    %default color is that of unclassified
    if isempty(image)
        trackSegmentColor = repmat([0 0 0],numTrackSegments,1);
    else
        trackSegmentColor = repmat([0.7 0.7 0.7],numTrackSegments,1);
    end
    %linear + 1D immobile
    indx = find(trackSegmentType(:,1) == 1 & trackSegmentType(:,3) == 0);
    switch simplifyLin
        case 0
            trackSegmentColor(indx,:) = repmat([0.3 0.3 0.3],length(indx),1);
        case {1,2}
            trackSegmentColor(indx,:) = repmat([1 0 0],length(indx),1);
    end
    %linear + 1D confined
    indx = find(trackSegmentType(:,1) == 1 & trackSegmentType(:,3) == 1);
    switch simplifyLin
        case 0
            trackSegmentColor(indx,:) = repmat([1 0.7 0],length(indx),1);
        case {1,2}
            trackSegmentColor(indx,:) = repmat([1 0 0],length(indx),1);
    end
    %linear + 1D normal
    indx = find(trackSegmentType(:,1) == 1 & trackSegmentType(:,3) == 2);
    trackSegmentColor(indx,:) = repmat([1 0 0],length(indx),1);
    indx = find(trackSegmentType(:,1) == 1 & trackSegmentType(:,3) == 3);
    switch simplifyLin
        case 0
            trackSegmentColor(indx,:) = repmat([0 1 0],length(indx),1);
        case {1,2}
            trackSegmentColor(indx,:) = repmat([1 0 0],length(indx),1);
    end
    %linear + 1D super
    indx = find(trackSegmentType(:,1) == 1 & isnan(trackSegmentType(:,3)));
    switch simplifyLin
        case 0
            trackSegmentColor(indx,:) = repmat([1 1 0],length(indx),1);
        case {1,2}
            trackSegmentColor(indx,:) = repmat([1 0 0],length(indx),1);
    end
    %random/unclassified + 2D immobile[0.5 0.3 0]
    indx = find(trackSegmentType(:,1) ~= 1 & trackSegmentType(:,2) == 0);
    trackSegmentColor(indx,:) = repmat([0 0 0],length(indx),1);
    %random/unclassified + 2D confined
    indx = find(trackSegmentType(:,1) ~= 1 & trackSegmentType(:,2) == 1);
    trackSegmentColor(indx,:) = repmat([0 0.5 0],length(indx),1);
    %random/unclassified + 2D normal
    indx = find(trackSegmentType(:,1) ~= 1 & trackSegmentType(:,2) == 2);
    trackSegmentColor(indx,:) = repmat([0 1 1],length(indx),1);
    %random/unclassified + 2D super
    indx = find(trackSegmentType(:,1) ~= 1 & trackSegmentType(:,2) == 3);
    switch simplifyLin
        case {0,1}
            trackSegmentColor(indx,:) = repmat([1 0 1],length(indx),1);
        case 2
            trackSegmentColor(indx,:) = repmat([1 0 0],length(indx),1);
    end
    %random + 2D unclassified
    indx = find(trackSegmentType(:,1) == 0 & isnan(trackSegmentType(:,2)));
    trackSegmentColor(indx,:) = repmat([0.6 0 1],length(indx),1);
    
else %output of classifyTrackDiffMode
    
    %get track segment diffusion modes
    trackSegmentType = vertcat(diffAnalysisRes.diffMode);
    
    %color coding:
    %       Mode 1 -> blue
    %       Mode 2 -> cyan
    %       Mode 3 -> red
    %       Mode 4 -> green
    %       Mode 5 -> magenta
    %       too short to classify -> black, or white if overlaying on an image
    
    %find number of diffusion modes
    numMode = max(trackSegmentType);
    
    %exit if there are more than 5 modes
    if numMode > 5
        disp('--plotTracksDiffAnalysis2D: Code cannot currently color-code more than 5 modes!');
        return
    end
    
    %the default, for unclassified tracks
    if isempty(image)
        trackSegmentColor = repmat([0 0 0],numTrackSegments,1);
    else
        trackSegmentColor = repmat([0.7 0.7 0.7],numTrackSegments,1);
    end
    
    %go over the modes
    indx = find(trackSegmentType==1);
    trackSegmentColor(indx,:) = repmat([0 0 1],length(indx),1);
    indx = find(trackSegmentType==2);
    trackSegmentColor(indx,:) = repmat([0 1 1],length(indx),1);
    indx = find(trackSegmentType==3);
    trackSegmentColor(indx,:) = repmat([1 0 0],length(indx),1);
    indx = find(trackSegmentType==4);
    trackSegmentColor(indx,:) = repmat([0 1 0],length(indx),1);
    indx = find(trackSegmentType==5);
    trackSegmentColor(indx,:) = repmat([1 0 1],length(indx),1);
    
end

%% confinement radius information

if showConf
    
    %get track segment center, confinement radii and preferred direction of
    %motion
    trackSegmentCenter = catStruct(1,'diffAnalysisRes.confRadInfo.trackCenter');
    trackSegmentConfRad = catStruct(1,'diffAnalysisRes.confRadInfo.confRadius');
    trackSegmentPrefDir = catStruct(1,'diffAnalysisRes.confRadInfo.prefDir');
    
    %determine indices of tracks with one confinement radius
    indxCircle = find( ~isnan(trackSegmentConfRad(:,1)) & isnan(trackSegmentConfRad(:,2)) );
    
    %determine indices of tracks with 2 confinement radii
    indxRectangle = find( ~isnan(trackSegmentConfRad(:,2)) );
    
end

%% Plotting

%if the user wants to plot in a new figure window
if newFigure

    %open new figure window
    figure

    if ~isempty(image) %if user supplied an image
        tmpImage = image(image~=0);
        minImage = min(tmpImage);
        maxImage = max(tmpImage);
        imshow(image,[minImage maxImage]); %plot the image
    else %if user did not supply an image
        imshow(ones(maxYCoord,maxXCoord),[]); %plot an empty image
    end

    %set figure axes limits
    axis([minXCoord maxXCoord minYCoord maxYCoord]);

    %show coordinates on axes
    ah = gca;
    set(ah,'visible','on');

    %label axes
    xlabel('x-coordinate (pixels)');
    ylabel('y-coordinate (pixels)');

end

%hold on figure
hold on

%extract the portion of tracksX and tracksY that is of interest
tracksXP = tracksX(timeRange(1):timeRange(2),:);
tracksYP = tracksY(timeRange(1):timeRange(2),:);

%plot tracks with their appropriate line color
%missing intervals are indicated by a dotted line
for i = 1 : numTrackSegments
    obsAvail = find(~isnan(tracksXP(:,i)));
    plot(tracksXP(obsAvail,i),tracksYP(obsAvail,i),'k:');
    plot(tracksXP(:,i),tracksYP(:,i),'Color',trackSegmentColor(i,:));
end

%show merges and splits
if mergeSplit

    %go over all tracks
    for iTrack = 1 : numTracks

        %parse sequence of events of this compound track and find merges and
        %splits
        seqOfEvents = inputStructure(iTrack).seqOfEvents;
        indxSplit = (find(seqOfEvents(:,2) == 1 & ~isnan(seqOfEvents(:,4)) ...
            & seqOfEvents(:,1) > timeRange(1) & seqOfEvents(:,1) <= timeRange(2)))';
        indxMerge = (find(seqOfEvents(:,2) == 2 & ~isnan(seqOfEvents(:,4)) ...
            & seqOfEvents(:,1) > timeRange(1) & seqOfEvents(:,1) <= timeRange(2)))';

        %go over all splits
        for iSplit = indxSplit

            %get time of splitting
            timeSplit = seqOfEvents(iSplit,1);

            %determine row where starting track is located
            rowS = trackStartRow(iTrack) + seqOfEvents(iSplit,3) - 1;

            %determine row where splitting track is located
            rowSp = trackStartRow(iTrack) + seqOfEvents(iSplit,4) - 1;

            %plot split as a dash-dotted line
            plot([tracksX(timeSplit,rowS) tracksX(timeSplit-1,rowSp)], ...
                [tracksY(timeSplit,rowS) tracksY(timeSplit-1,rowSp)],'k-.');

        end

        %go over all merges
        for iMerge = indxMerge

            %get time of merging
            timeMerge = seqOfEvents(iMerge,1);

            %determine row where ending track is located
            rowE = trackStartRow(iTrack) + seqOfEvents(iMerge,3) - 1;

            %determine row where merging track is located
            rowM = trackStartRow(iTrack) + seqOfEvents(iMerge,4) - 1;

            %plot merge as a dashed line
            plot([tracksX(timeMerge-1,rowE) tracksX(timeMerge,rowM)], ...
                [tracksY(timeMerge-1,rowE) tracksY(timeMerge,rowM)],'k--');

        end

    end %(for iTrack = 1 : numTracks)

end %(if mergeSplit)

%show confinement areas if requested
if showConf

    if isempty(image)
        confColor = [0 0 0];
    else
        confColor = [0.7 0.7 0.7];
    end
    
    %generate circle to plot
    theta = (0:pi/10:2*pi); %angle
    xy = [cos(theta') sin(theta')]; %x and y-coordinates
    
    %go over symmetric tracks
    for iTrack = indxCircle'
        
        %plot a circle of radius = confinement radius and centered at the
        %center of this track
        circleVal = xy .* trackSegmentConfRad(iTrack,1);
        plot(trackSegmentCenter(iTrack,1)+circleVal(:,1),...
            trackSegmentCenter(iTrack,2)+circleVal(:,2),'Color',trackSegmentColor(iTrack,:));
        
    end
    
    %go over linear tracks
    for iTrack = indxRectangle'
    
        %get the confinement axes
        axisPara = trackSegmentPrefDir(iTrack,:);
        axisPerp = [-axisPara(2) axisPara(1)] * trackSegmentConfRad(iTrack,1);
        axisPara = axisPara * trackSegmentConfRad(iTrack,2);
        
        %find the 4 corners of the confinement rectangle
        cornerCoord = [-axisPara - axisPerp; -axisPara + axisPerp; ...
            axisPara + axisPerp; axisPara - axisPerp; -axisPara - axisPerp] ...
            + repmat(trackSegmentCenter(iTrack,:),5,1);
        
        %plot the rectangle
        plot(cornerCoord(:,1),cornerCoord(:,2),'Color',confColor);

    end
    
end

%%%%% ~~ the end ~~ %%%%%

