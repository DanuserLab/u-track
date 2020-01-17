function plotTracks2D(trackedFeatureInfo, varargin)
%function plotTracks2D(trackedFeatureInfo,timeRange,colorTime,markerType,...
%    indicateSE,newFigure,image,flipXY,ask4sel,offset,minLength)
%PLOTTRACKS2D plots a group of tracks in 2D and allows user to click on them and extract track information
%
%SYNOPSIS plotTracks2D(trackedFeatureInfo,timeRange,colorTime,markerType,...
%    indicateSE,newFigure,image,flipXY,ask4sel,offset)
%
%INPUT  trackedFeatureInfo: -- EITHER --
%                           Output of trackWithGapClosing:
%                           Matrix indicating the positions and amplitudes
%                           of the tracked features to be plotted. Number
%                           of rows = number of tracks, while number of
%                           columns = 8*number of time points. Each row consists of
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
%       timeRange         : 2-element row vector indicating time range to plot.
%                           Optional. Default: whole movie.
%       colorTime         : String with the following options:
%                           -'1' if time is to be color-coded (green in the
%                           beginning, blue in the middle, red in the end).
%                           -'k', 'b', 'r', etc. if all tracks are in black,
%                              blue, red, etc.
%                           -'2' if tracks are colored by cycling through
%                              the plot's default color order.
%                           -'3' as 2, but using extendedColors (23
%                              different colors).
%                           Optional. Default: 'k'.
%       markerType        : String indicating marker type for plotting.
%                           Only used if colorTime is not '1'.
%                           Optional. Default: 'none'.
%       indicateSE        : 1 if track starts and ends are to be indicated
%                           with circles and squares, respectively; 0
%                           otherwise. Optional. Default: 1.
%       newFigure         : 1 if plot should be made in a new figure
%                           window, 0 otherwise (in which case it will be
%                           plotted in an existing figure window).
%                           newFigure can also be a handle to the axes in
%                           which to plot.
%                           Optional. Default: 1.
%       image             : An image that the tracks will be overlaid on if
%                           newFigure=1. It will be ignored if newFigure=0.
%                           Optional. Default: no image
%       flipXY            : 1 if x and y coord should be flipped for
%                           plotting. Optional. Default: 0.
%       ask4sel           : 1 if user should be asked to select tracks in
%                           plot in order to show track information.
%                           Optional. Default: 1.
%       offset            : [dx,dy] that is to be added to the coordinates.
%                           Optional. Default: [0,0]
%       minLength         : Minimum length of tracks to be ploted.
%                           Optional. Default: 1.
%
%OUTPUT The plot.
%
%Khuloud Jaqaman, August 2006
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

ip = inputParserRetrofit;
ip.addRequired('trackedFeatureInfo', ...
    @(s) isstruct(s) || isnumeric(s));
ip.addArgument('timeRange',[1 numTimePoints], ...
    @(t) all(size(t) == [1 2]) && ...
    t(1) >= 1 && t(2) <= numTimePoints );
ip.addArgument('colorTime','k', ...
    @isscalar);
ip.addArgument('markerType','none', ...
    @ischar);
ip.addArgument('indicateSE',1, ...
    @(x) ismember(x,[ 0 1 ]));
ip.addArgument('newFigure',1, ...
    @(x) ismember(x,[0 1]) || ...
    ishandle(x) && strncmp(get(x,'Type'),'axes',4));
ip.addArgument('image',[], ...
    @isnumeric);
ip.addArgument('flipXY',false, ...
    @(x) ismember(x,[0 1]) );
ip.addArgument('ask4sel',true, ...
    @(x) ismember(x,[0 1]));
ip.addArgument('offset',[ 0 0 ], ...
    @(x) isnumeric(x) && all(size(x) == [1 2]));
ip.addArgument('minLength',1, ...
    @(x) isnumeric(x) && isscalar(x));
ip.parse(trackedFeatureInfo,varargin{:});

trackedFeatureInfo = ip.Results.trackedFeatureInfo;
timeRange = ip.Results.timeRange;
colorTime = ip.Results.colorTime;
markerType = ip.Results.markerType;
indicateSE = ip.Results.indicateSE;
newFigure = ip.Results.newFigure;
image = ip.Results.image;
flipXY = ip.Results.flipXY;
ask4sel = ip.Results.ask4sel;
offset = ip.Results.offset;
minLength = ip.Results.minLength;

%keep only tracks with minimum requested length
if minLength > 1
    criteria.lifeTime.min = minLength;
    indx = chooseTracks(trackedFeatureInfo,criteria);
    trackedFeatureInfo = trackedFeatureInfo(indx,:);
end


% make sure colorTime 1,2 are strings
if isnumeric(colorTime) && isscalar(colorTime)
    colorTime = num2str(colorTime);
end

% check for axes handle
if ishandle(newFigure) && strcmp(get(newFigure,'Type'),'axes')
    axH = newFigure;
    newFigure = 0;
elseif newFigure == 0
    axH = gca;
end


%% Pre-processing

if isstruct(trackedFeatureInfo) %if tracks are input in structure format
    
    %store the input structure as a variable with a different name
    inputStructure = trackedFeatureInfo;
    clear trackedFeatureInfo;
    
    [trackedFeatureInfo,~,trackStartRow,numSegments] = convStruct2MatIgnoreMS(inputStructure);
    
    %indicate merge/split flag
    if max(numSegments) > 1
        mergeSplit = 1;
    else
        mergeSplit = 0;
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
if flipXY
    tracksY = trackedFeatureInfo(:,1:8:end)' + offset(1);
    tracksX = trackedFeatureInfo(:,2:8:end)' + offset(2);
else
    tracksX = trackedFeatureInfo(:,1:8:end)' + offset(1);
    tracksY = trackedFeatureInfo(:,2:8:end)' + offset(2);
end

% % add offset
% tracksX = tracksX ;
% tracksY = tracksY ;

%find x-coordinate limits
minXCoord = min(floor(min(tracksX(:))),0);
maxXCoord =  ceil(max(tracksX(:)));

%find y-coordinate limits
minYCoord = min(floor(min(tracksY(:))),0);
maxYCoord =  ceil(max(tracksY(:)));

%calculate the number of time points to be plotted
numTimePlot = timeRange(2) - timeRange(1) + 1;

%define colors to loop through in case colorTime = '2'
if isempty(image)
    colorLoop = [0 0 0; 1 0 0; 0 1 0; 0 0 1; 1 1 0; 1 0 1; 0 1 1]; %colors: k,r,g,b,y,m,c
else
    colorLoop = [0.7 0.7 0.7; 1 0 0; 0 1 0; 0 0 1; 1 1 0; 1 0 1; 0 1 1]; %colors: 'light pink',r,g,b,y,m,c
end

%% Plotting

%if the user wants to plot in a new figure window
if newFigure
    
    %open new figure window
    figure;
    
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
    axH = gca;
    set(axH,'visible','on');
    
    %label axes
    if flipXY
        xlabel('y-coordinate (pixels)');
        ylabel('x-coordinate (pixels)');
    else
        xlabel('x-coordinate (pixels)');
        ylabel('y-coordinate (pixels)');
    end
    
end

%hold on axes
set(axH,'NextPlot','add');

removeGapParams.Delimeter = Inf;
removeGapParams.RemoveOtherGaps = true;
removeGapParams.Color = 'k';
removeGapParams.LineStyle = ':';

%extract the portion of tracksX and tracksY that is of interest
tracksXP = tracksX(timeRange(1):timeRange(2),:);
tracksYP = tracksY(timeRange(1):timeRange(2),:);

axes(axH);

switch colorTime
    
    case '1' %if user wants to color-code time
        
        lineWithGaps(tracksXP,tracksYP,removeGapParams);
        
        %get the overall color per time interval
        colorOverTime = timeColormap(numTimePlot);
        
        %overlay tracks with color coding wherever a feature has been detected
        for i=1:numTimePlot-1
            validData=~all(isnan(tracksXP(i:i+1,:)),1);
            lineWithGaps(tracksXP(i:i+1,validData), ...
                tracksYP(i:i+1,validData), ...
                'Color',colorOverTime(i,:));
        end
        
    case '2' %no time color-coding, loop through series of colors to color tracks
        
        %plot tracks by looping through colors
        %missing intervals are indicated by a dotted line
        lineWithGaps(tracksXP,tracksYP,removeGapParams);
        for i=1:7
            lineWithGaps(tracksXP(:,i:7:end),tracksYP(:,i:7:end),'Color',colorLoop(i,:),'Marker',markerType);
        end
        
    case '3' % no time color-coding, use extendedColors
        
        %plot tracks by looping through colors
        %missing intervals are indicated by a dotted line
        lineWithGaps(tracksXP,tracksYP,removeGapParams);
        for i=1:23
            lineWithGaps(tracksXP(:,i:23:end),tracksYP(:,i:23:end),'Color',extendedColors(i),'Marker',markerType);
        end
        
    otherwise %no time color-coding, all tracks same color
        
        %plot tracks with the line color indicated
        %missing intervals are indicated by a dotted line
        lineWithGaps(tracksXP,tracksYP,removeGapParams);
        lineWithGaps(tracksXP,tracksYP,'Color',colorTime,'Marker',markerType);
        
end %(switch colorTime)

%show merges and splits
if mergeSplit
    
    split.X = [];
    split.Y = [];
    merge.X = [];
    merge.Y = [];
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
            
            split.X = [split.X tracksX(timeSplit,rowS) tracksX(timeSplit-1,rowSp) NaN];
            split.Y = [split.Y tracksY(timeSplit,rowS) tracksY(timeSplit-1,rowSp) NaN];
            
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
            merge.X = [merge.X tracksX(timeMerge-1,rowE) tracksX(timeMerge,rowM) NaN];
            merge.Y = [merge.Y tracksY(timeMerge-1,rowE) tracksY(timeMerge,rowM) NaN];
            
        end
        
    end %(for iTrack = 1 : numTracks)
    
    %plot split as a dash-dotted line
    plot(axH,split.X, split.Y,'k-.');
    %plot merge as a dashed line
    plot(axH,merge.X, merge.Y,'k--');
    
end %(if mergeSplit)

if indicateSE %if user wants to indicate starts and ends
    
    allStartInfo = [];
    allEndInfo = [];
    %if there are merges and splits
    if mergeSplit
        
        %go over all tracks
        for iTrack = 1 : numTracks
            
            %parse sequence of events of this compound track and find starts and
            %ends
            seqOfEvents = inputStructure(iTrack).seqOfEvents;
            indxStart = (find(seqOfEvents(:,2) == 1 & isnan(seqOfEvents(:,4)) ...
                & seqOfEvents(:,1) >= timeRange(1) & seqOfEvents(:,1) <= timeRange(2)))';
            indxEnd = (find(seqOfEvents(:,2) == 2 & isnan(seqOfEvents(:,4)) ...
                & seqOfEvents(:,1) >= timeRange(1) & seqOfEvents(:,1) <= timeRange(2)))';
            
            %get the information of the starts
            startInfo = zeros(length(indxStart),3);
            for i = 1 : length(indxStart)
                iStart = indxStart(i);
                
                %get start time
                timeStart = seqOfEvents(iStart,1);
                
                %determine row where starting track is located in big matrix
                %of tracks
                rowS = trackStartRow(iTrack) + seqOfEvents(iStart,3) - 1;
                
                %get coordinates at the start
                startInfo(i,:) = [tracksX(timeStart,rowS) tracksY(timeStart,rowS) timeStart];
                
            end
            
            %get the information of the ends
            endInfo = zeros(length(indxEnd),3);
            for i = 1 : length(indxEnd)
                iEnd = indxEnd(i);
                
                %get end time
                timeEnd = seqOfEvents(iEnd,1);
                
                %determine row where ending track is located in big matrix
                %of tracks
                rowE = trackStartRow(iTrack) + seqOfEvents(iEnd,3) - 1;
                
                %get coordinates at the end
                endInfo(i,:) = [tracksX(timeEnd,rowE) tracksY(timeEnd,rowE) timeEnd];
                
            end
            
            allStartInfo = [ allStartInfo ; startInfo ];
            allEndInfo = [ allEndInfo ; endInfo ];
            
        end %(for iTrack = 1 : numTracks)
        
    else %if there are no merges and splits
        
        %find the beginning and end of each track
        for i=numTracks:-1:1
            timePoint = find(~isnan(tracksX(:,i)));
            if(~isempty(timePoint))
                startInfo(i,:) = [tracksX(timePoint(1),i) ...
                    tracksY(timePoint(1),i) timePoint(1)];
                endInfo(i,:) = [tracksX(timePoint(end),i) ...
                    tracksY(timePoint(end),i) timePoint(end)];
            end
        end
        
        allStartInfo = [ allStartInfo ; startInfo(startInfo(:,3)>=timeRange(1) & startInfo(:,3)<=timeRange(2),:) ];
        allEndInfo = [ allEndInfo ; endInfo(endInfo(:,3)>=timeRange(1) & endInfo(:,3)<=timeRange(2),:) ];
        
    end %(if mergeSplit)
    
    switch colorTime
        case {'1','2','3'}
            scatterColor = 'k';
        otherwise
            scatterColor = colorTime;
    end
    scatter(axH,allStartInfo(:,1),allStartInfo(:,2),[],scatterColor,'o');
    scatter(axH,allEndInfo(:,1),allEndInfo(:,2),[],scatterColor,'square');
    
end %(if indicateSE)

if ask4sel
    hSelectButton = uicontrol('Style','pushbutton','String','Select Track','Callback',@selectPoint);
    pos = get(hSelectButton,'Position');
    pos(3) = 120;
    set(hSelectButton,'Position',pos);
end



    function selectPoint(~,~)
        disp('Select multiple points using the mouse. Double click on the last point to finish.');
        
        %let the user choose the points of interest
        [x,y] = getpts;
        
        %find the time points of the indicated points
        for npt=1:length(x)
            
            %find the distances between those points and the tracks
            distTrack2Point = (tracksXP-x(npt)).^2+(tracksYP-y(npt)).^2;
            
            %determine the minimum distance for each chosen point
            [frameChosen,rowChosen] = find(distTrack2Point==min(distTrack2Point(:)));
            
            %go over all chosen rows
            for j = 1 : length(rowChosen)
                
                %find the track corresponding to each minimum distance
                trackChosen = find(trackStartRow <= rowChosen(j),1,'last');
                segmentChosen = rowChosen(j) - trackStartRow(trackChosen) + 1;
                
                disp(['Track: ' num2str(trackChosen) ...
                    '   Segment: ' num2str(segmentChosen) ...
                    '   Frame: ' num2str(frameChosen(j)+timeRange(1)-1) ...
                    '   Coordinates: ' num2str(tracksXP(frameChosen(j),rowChosen(j))) ...
                    ' ' num2str(tracksYP(frameChosen(j),rowChosen(j)))  ]);
                
            end
            
        end
    end

end

%% %%%%% ~~ the end ~~ %%%%%

