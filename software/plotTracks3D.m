function axH = plotTracks3D(trackedFeatureInfo,timeRange,colorTime,markerType,...
    indicateSE, newFigure, figName)
%PLOTTRACKS3D plots a group of tracks in 3D with time coloring
%
%SYNOPSIS plotTracks3D(trackedFeatureInfo,timeRange,colorTime,markerType,...
%    indicateSE,newFigure)
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
%       timeRange         : 2-element row vector indicating time range to plot. 
%                           Optional. Default: whole movie.
%       colorTime         : String with the following options:
%                           -'1' if time is to be color-coded (green in the
%                           beginning, blue in the middle, red in the end).
%                           -'k', 'b', 'r', etc. if all tracks are in black,
%                           blue, red, etc.
%                           - 'R' color randomizer
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
%       figName           : Figure Name. 
%                           Default: "PlotTracks3D_wholeMovie"
%
%OUTPUT The plot 
% returns the figure handle h.
%
%Khuloud Jaqaman, July 2007
%edits Philippe Roudot, October 2014
% adapted for 3DUtrackPackage GUI simple figure display of output.
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%check whether correct number of input arguments was used
if nargin < 1
    disp('--plotTracks3D: Incorrect number of input arguments!');
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
if nargin < 2 || isempty(timeRange)
    timeRange = [1 numTimePoints];
else
    if timeRange(1) < 1 || timeRange(2) > numTimePoints
        disp('--plotTracks3D: Wrong time range for plotting!');
        errFlag = 1;
    end
end

%check whether colorTime was input
if nargin < 3 || isempty(colorTime)
    colorTime = 'k';
end

%check whether markerType was input
if nargin < 4 || isempty(markerType)
    markerType = 'none';
end

%check whether indicateSE was input
if nargin < 5 || isempty(indicateSE)
    indicateSE = 1;
else
    if indicateSE ~= 0 && indicateSE ~= 1
        disp('plotTracks3D: indicateSE should be 0 or 1!');
        errFlag = 1;
    end
end

%check whether newFigure was input
if nargin < 6 || isempty(newFigure)
    newFigure = 1;
else
    if newFigure ~= 0 && newFigure ~= 1 && ~(ishandle(newFigure) && strmatch(get(newFigure,'Type'),'axes'))
        disp('--plotTracks3D: newFigure should be 0 or 1 or an axes handle!');
        errFlag = 1;
    end
end



% check for axes handle
if ishandle(newFigure) && strcmp(get(newFigure,'Type'),'axes')
    axH = newFigure;
    newFigure = 0;
elseif newFigure == 0
    axH = gca;
end

%check whether newFigure was input
if nargin < 7 || isempty(figName)
    figName = 'PlotTracks3D';
end



%exit if there are problem in input variables
if errFlag
    disp('--plotTracks3D: Please fix input data!');
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Pre-processing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
    
    %indicate that each track consists of one segment
    numSegments = ones(numTracks,1);

    %locate the row of the first track of each compound track in the
    %big matrix of all tracks
    %in this case of course every compound track is simply one track
    %without branches
    trackStartRow = (1:numTracks)';

end

%get the coordinates of features in all tracks
tracksX = trackedFeatureInfo(:,1:8:end)';
tracksY = trackedFeatureInfo(:,2:8:end)';
tracksZ = trackedFeatureInfo(:,3:8:end)';

%calculate the number of time points to be plotted
numTimePlot = timeRange(2) - timeRange(1) + 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Plotting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%open a new figure window if requested
if newFigure
    hfig = figure;
    % read axes handle of new figure
    axH = gca;
end


%hold on figure
hold on

%extract the portion of tracksX, tracksY and tracksZ that is of interest
tracksXP = tracksX(timeRange(1):timeRange(2),:);
tracksYP = tracksY(timeRange(1):timeRange(2),:);
tracksZP = tracksZ(timeRange(1):timeRange(2),:);

switch colorTime
    case '1' %if user wants to color-code time

        %plot tracks ignoring missing points
        %gaps are depicted as a dotted black line
        for i = 1 : trackStartRow(end) + numSegments(end) - 1
            obsAvail = find(~isnan(tracksXP(:,i)));
            plot3(axH,tracksXP(obsAvail,i),tracksYP(obsAvail,i),tracksZP(obsAvail,i),'k:');
        end

        %get the fraction of each color in each time interval to be plotted
        numTimePlotOver2 = ceil((numTimePlot-1)/2); %needed to change blue color over time
        redVariation = (0:numTimePlot-2)'/(numTimePlot-2);
        greenVariation = (numTimePlot-2:-1:0)'/(numTimePlot-2);
        blueVariation = [(0:numTimePlotOver2-1)'/(numTimePlotOver2-1);...
            (numTimePlot-numTimePlotOver2-2:-1:0)'/(numTimePlot-numTimePlotOver2-1)];

        %get the overall color per time interval
        colorOverTime = [redVariation greenVariation blueVariation];

        %overlay tracks with color coding wherever a feature has been detected
        for i=1:numTimePlot-1
            plot3(axH,tracksXP(i:i+1,:),tracksYP(i:i+1,:),tracksZP(i:i+1,:),'color',colorOverTime(i,:));
        end
    case 'R'
        for i = 1 : trackStartRow(end) + numSegments(end) - 1
            %obsAvail = find(~isnan(tracksXP(:,i)));
            %plot3(axH,tracksXP(obsAvail,i),tracksYP(obsAvail,i),tracksZP(obsAvail,i),[colorTime ':']);
            
        end
        plot3(axH,tracksXP,tracksYP,tracksZP);

    otherwise
        disp('test')
        %plot tracks with the line color indicated, where missing intervals are
        %indicated by a dotted line
        for i = 1 : trackStartRow(end) + numSegments(end) - 1
            obsAvail = find(~isnan(tracksXP(:,i)));
            plot3(axH,tracksXP(obsAvail,i),tracksYP(obsAvail,i),tracksZP(obsAvail,i),[colorTime ':']);
            plot3(axH,tracksXP(:,i),tracksYP(:,i),tracksZP(:,i),colorTime,'marker',markerType);
        end

end %(if colorTime == '1' ... else ...)

%label axes
xlabel('x-coordinate (pixels)');
ylabel('y-coordinate (pixels)');
zlabel('z-coordinate (pixels)');

%show merges and splits
if mergeSplit

    %draw line in black if time color coding was used initially
    if colorTime == '1'
        colorTime = 'k';
    end

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
            plot3(axH,[tracksX(timeSplit,rowS) tracksX(timeSplit-1,rowSp)], ...
                [tracksY(timeSplit,rowS) tracksY(timeSplit-1,rowSp)],...
                [tracksZ(timeSplit,rowS) tracksZ(timeSplit-1,rowSp)],...
                [colorTime '-.']);

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
            plot3(axH,[tracksX(timeMerge-1,rowE) tracksX(timeMerge,rowM)], ...
                [tracksY(timeMerge-1,rowE) tracksY(timeMerge,rowM)],...
                [tracksZ(timeMerge-1,rowE) tracksZ(timeMerge,rowM)],...
                [colorTime '--']);

        end

    end %(for iTrack = 1 : numTracks)

end %(if mergeSplit)

if indicateSE %if user wants to indicate starts and ends

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
            startInfo = [];
            for i = 1 : length(indxStart)
                iStart = indxStart(i);

                %get start time
                timeStart = seqOfEvents(iStart,1);

                %determine row where starting track is located in big matrix
                %of tracks
                rowS = trackStartRow(iTrack) + seqOfEvents(iStart,3) - 1;

                %get coordinates at the start
                startInfo(i,:) = [tracksX(timeStart,rowS) tracksY(timeStart,rowS) ...
                    tracksZ(timeStart,rowS) timeStart];

            end

            %get the information of the ends
            endInfo = [];
            for i = 1 : length(indxEnd)
                iEnd = indxEnd(i);
                
                %get end time
                timeEnd = seqOfEvents(iEnd,1);

                %determine row where ending track is located in big matrix
                %of tracks
                rowE = trackStartRow(iTrack) + seqOfEvents(iEnd,3) - 1;

                %get coordinates at the end
                endInfo(i,:) = [tracksX(timeEnd,rowE) tracksY(timeEnd,rowE) ...
                    tracksZ(timeEnd,rowE) timeEnd];

            end

            %place circles at track starts and squares at track ends
            if ~isempty(startInfo)
                plot3(axH,startInfo(:,1),startInfo(:,2),startInfo(:,3),colorTime,...
                    'LineStyle','none','marker','o');
            end
            if ~isempty(endInfo)
                plot3(axH,endInfo(:,1),endInfo(:,2),endInfo(:,3),colorTime,...
                    'LineStyle','none','marker','square');
            end

        end %(for iTrack = 1 : numTracks)

    else %if there are no merges and splits

        %find the beginning and end of each track
        for i=numTracks:-1:1
            timePoint = find(~isnan(tracksX(:,i)));
            startInfo(i,:) = [tracksX(timePoint(1),i) ...
                tracksY(timePoint(1),i) tracksZ(timePoint(1),i) timePoint(1)];
            endInfo(i,:) = [tracksX(timePoint(end),i) ...
                tracksY(timePoint(end),i) tracksZ(timePoint(end),i) timePoint(end)];
        end

        %place circles at track starts and squares at track ends if they happen to
        %be in the plotting region of interest
        if colorTime == '1'
            indx = find(startInfo(:,4)>=timeRange(1) & startInfo(:,4)<=timeRange(2));
            plot3(axH,startInfo(indx,1),startInfo(indx,2),startInfo(indx,3),'k','LineStyle','none','marker','o');
            indx = find(endInfo(:,4)>=timeRange(1) & endInfo(:,4)<=timeRange(2));
            plot3(axH,endInfo(indx,1),endInfo(indx,2),endInfo(indx,3),'k','LineStyle','none','marker','square');
        else
            indx = find(startInfo(:,4)>=timeRange(1) & startInfo(:,4)<=timeRange(2));
            plot3(axH,startInfo(indx,1),startInfo(indx,2),startInfo(indx,3),colorTime,...
                'LineStyle','none','marker','o');
            indx = find(endInfo(:,4)>=timeRange(1) & endInfo(:,4)<=timeRange(2));
            plot3(axH,endInfo(indx,1),endInfo(indx,2),endInfo(indx,3),colorTime,...
                'LineStyle','none','marker','square');
        end

    end %(if mergeSplit)

end %(if indicateSE)
        

%%%%% ~~ the end ~~ %%%%%

% % %I guess this only works for 2D
% % 
% % %ask the user whether to click on figure and get frame information
% % userEntry = input('select points in figure? y/n ','s');
% % 
% % while strcmp(userEntry,'y')
% % 
% %     %let the user choose the points of interest
% %     [x,y] = getpts;
% % 
% %     %find the time points of the indicated points
% %     for i=1:length(x)
% %         
% %         %find the distances between those points and the tracks
% %         distTrack2Point = (tracksXP-x(i)).^2+(tracksYP-y(i)).^2;
% %         
% %         %determine the minimum distance for each chosen point
% %         [frameChosen,rowChosen] = find(distTrack2Point==min(distTrack2Point(:)));
% %         
% %         %go over all chosen rows
% %         for j = 1 : length(rowChosen)
% %             
% %             %find the track corresponding to each minimum distance
% %             trackChosen = find(trackStartRow <= rowChosen(j),1,'last');
% %             segmentChosen = rowChosen(j) - trackStartRow(trackChosen) + 1;
% %         
% %             disp(['Track: ' num2str(trackChosen) ...
% %                 '   Segment: ' num2str(segmentChosen) ...
% %                 '   Frame: ' num2str(frameChosen(j)+timeRange(1)-1) ...
% %                 '   Coordinates: ' num2str(tracksXP(frameChosen(j),rowChosen(j))) ...
% %                 ' ' num2str(tracksYP(frameChosen(j),rowChosen(j)))  ]);
% %             
% %         end
% %         
% %     end
% %         
% %     %ask the user again whether to click on figure and get frame information
% %     userEntry = input('select points again? y/n ','s');
% % 
% % end

