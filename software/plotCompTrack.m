function plotCompTrack(trackedFeatureInfo,varargin)
% function plotCompTrack(trackedFeatureInfo,plotX,plotY,plotA,inOneFigure,...
%     plotAS,timeStep,markMS)
%PLOTCOMPTRACK plots the x-coordinates, y-coordinates and/or intensities along a compound track, indicating merges, splits and gaps
%
%SYNOPSIS plotCompTrackAmp(trackedFeatureInfo,plotX,plotY,plotA,inOneFigure,...
%    plotAS,timeStep,markMS)
%
%INPUT  trackedFeatureInfo: Output of trackCloseGapsKalman for one track:
%                           Contains the fields:
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
%           .aggregationState: This field results from running the function
%                              aggregStateFromCompTracks. Only needed if
%                              aggregation state is to be plotted.
%       plotX             : 1 if x-coordinate is to be plotted, zero
%                           otherwise. Optional. Default: 1.
%       plotY             : 1 if y-coordinate is to be plotted, zero
%                           otherwise. Optional. Default: 1.
%       plotA             : 1 if amplitude is to be plotted, zero
%                           otherwise. Optional. Default: 1.
%       inOneFigure       : 1 if all plots appear in one figure window (one
%                           above the other), 0 if each figure is in its
%                           own window. Optional. Default: 1.
%       plotAS   : 1 to plot particle aggregation state (if
%                           supplied), 0 otherwise. Optional. Default: 1.
%       timeStep          : Time step for x-axis. Can be whatever units are
%                           relevant. Optional. Default: 1.
%                           NOTE: If timeStep is supplie, x-axis starts at
%                                 0. If not supplied, x-axis starts at 1.
%       markMS            : 1 to mark merges and splits, 0 otherwise.
%                           Optional. Default: 0.
%
%OUTPUT The plot(s).
%
%REMARKS Gaps are dotted black lines, splits are dash-dotted black lines
%        and merges are dashed lines. If requested, splits and merges are
%        also marked as black "+" and red "x" symbols, respectively.
%
%        KJ: Made a few changes on 16 Dec 2014. Not sure how these changes
%        will affect display in case of "double frequency" time. Will test
%        when case arises again.
%
%Khuloud Jaqaman, May 2007

%% Input

ip = inputParserRetrofit;
ip.addRequired('trackedFeatureInfo', ...
    @(s) isstruct(s) || isnumeric(s));
ip.addArgument('plotX',1, ...
    @(x) islogical(x) || ismember(x,[0 1]));
ip.addArgument('plotY',1, ...
    @(x) islogical(x) || ismember(x,[0 1]));
ip.addArgument('plotA',1, ...
    @(x) islogical(x) || ismember(x,[0 1]));
ip.addArgument('inOneFigure',1, ...
    @(x) islogical(x) || ismember(x,[0 1]));
ip.addArgument('plotAS',1, ...
    @(x) islogical(x) || ismember(x,[0 1]));
ip.addArgument('timeStep',[], ...
    @(x) isnumeric(x) || isempty(x));
ip.addArgument('markMS',0, ...
    @(x) islogical(x) || ismember(x,[0 1]));
ip.parse(trackedFeatureInfo, varargin{:});

assignFieldsHere(ip.Results);

if(isempty(ip.Results.timeStep))
    timeStep = 1;
    axisOffset = 0;
else
    axisOffset = 1;
end

%extract information from input
seqOfEvents = trackedFeatureInfo.seqOfEvents;
tracksCoordAmpCG = trackedFeatureInfo.tracksCoordAmpCG;
if isfield(trackedFeatureInfo,'aggregState')
    aggregState = trackedFeatureInfo.aggregState;
else
    plotAS = 0;
end
followerInfo = 0;
if isfield(trackedFeatureInfo,'followerCoordAmp')
    followerCoordAmp = trackedFeatureInfo.followerCoordAmp;
    followerInfo = 1;
end

%convert sparse to full if necessary
%code below is a copy-paste from trackCloseGapsKalmanSparse
if issparse(tracksCoordAmpCG)
    tracksCoordAmpCG = full(tracksCoordAmpCG);
    for iRow = 1 : size(tracksCoordAmpCG,1)
        colZero = find(tracksCoordAmpCG(iRow,:)==0);
        colZero = colZero(:)';
        xCoordCol = colZero - mod(colZero-1,8*ones(size(colZero)));
        colZero = colZero(tracksCoordAmpCG(iRow,xCoordCol)==0);
        tracksCoordAmpCG(iRow,colZero) = NaN;
    end
end

%determine whether track is sampled regularly or with doubled frequency
doubleFreq = any(mod(seqOfEvents(:,1)*2,2)==1);

%% Plotting

%find merge and split times if to be marked
if markMS
    splitTimes = seqOfEvents(seqOfEvents(:,2)==1&~isnan(seqOfEvents(:,4)),1);
    mergeTimes = seqOfEvents(seqOfEvents(:,2)==2&~isnan(seqOfEvents(:,4)),1);
end
    
%x-coordinates
if plotX

    %open new figure window and hold on to it
    figure, hold on

    %specify subplot
    if inOneFigure
        h = subplot(plotX+plotY+plotA+plotAS,1,1);
    else
        h = subplot(1,1,1);
    end
    hold on

    %extract x-coordinates from input
    xCoordSequence = tracksCoordAmpCG(:,1:8:end);
    if followerInfo
        xCoordFollower = followerCoordAmp(:,1:8:end);
    end

    %mark merges and splits
    if markMS
        yaxisMin = min(xCoordSequence(:));
        yaxisMax = max(xCoordSequence(:));
        plot((splitTimes-axisOffset)*timeStep,(yaxisMin-0.1*(yaxisMax-yaxisMin))*ones(size(splitTimes)),'k+');
        plot((mergeTimes-axisOffset)*timeStep,(yaxisMin-0.1*(yaxisMax-yaxisMin))*ones(size(mergeTimes)),'rx');
        legend({'Splits','Merges'})
    end
    
    %plot the x-coordinate, taking into account closed gaps, merges and
    %splits and sampling frequency doubling
    if followerInfo
        plotCompTrackCore(xCoordSequence,seqOfEvents,doubleFreq,timeStep,axisOffset,xCoordFollower,h);
    else
        plotCompTrackCore(xCoordSequence,seqOfEvents,doubleFreq,timeStep,axisOffset);
    end
    
    %put axes labels
    if ~inOneFigure
        if axisOffset
            xlabel('Time');
        else
            xlabel('Frame number');
        end
    end
    ylabel('X-coordinate (pixels)');

    %hold off of figure if each plot is in a separate figure window or if
    %this is the last plot
    if ~inOneFigure || (~plotY && ~plotA && ~plotAS)
        hold off
    end

end %(if plotX)

%y-coordinates
if plotY

    %open new figure window if needed and hold on to it
    if ~inOneFigure || ~plotX
        figure, hold on
    end

    %specify subplot
    if inOneFigure
        h = subplot(plotX+plotY+plotA+plotAS,1,plotX+1);
    else
        h = subplot(1,1,1);        
    end
    hold on

    %extract y-coordinates from input
    yCoordSequence = tracksCoordAmpCG(:,2:8:end);
    if followerInfo
        yCoordFollower = followerCoordAmp(:,2:8:end);
    end

    %mark merges and splits
    if markMS
        yaxisMin = min(yCoordSequence(:));
        yaxisMax = max(yCoordSequence(:));
        plot((splitTimes-axisOffset)*timeStep,(yaxisMin-0.1*(yaxisMax-yaxisMin))*ones(size(splitTimes)),'k+');
        plot((mergeTimes-axisOffset)*timeStep,(yaxisMin-0.1*(yaxisMax-yaxisMin))*ones(size(mergeTimes)),'rx');
        legend({'Splits','Merges'})
    end
    
    %plot the y-coordinate, taking into account closed gaps, merges and
    %splits and sampling frequency doubling
    if followerInfo
        plotCompTrackCore(yCoordSequence,seqOfEvents,doubleFreq,timeStep,axisOffset,yCoordFollower,h);
    else
        plotCompTrackCore(yCoordSequence,seqOfEvents,doubleFreq,timeStep,axisOffset);
    end

    %put axes labels
    if ~inOneFigure
        if axisOffset
            xlabel('Time');
        else
            xlabel('Frame number');
        end
    end
    ylabel('Y-coordinate (pixels)');

    %hold off of figure if each plot is in a separate figure window or if
    %this is the last plot
    if ~inOneFigure || (~plotA && ~plotAS)
        hold off
    end

end %(if plotY)

%amplitudes
if plotA

    %open new figure window if needed and hold on to it
    if ~inOneFigure || (~plotX && ~plotY)
        figure, hold on
    end

    %specify subplot
    if inOneFigure
        h = subplot(plotX+plotY+plotA+plotAS,1,plotX+plotY+1);
    else
        h = subplot(1,1,1);
    end
    hold on

    %extract amplitudes from input
    ampSequence = tracksCoordAmpCG(:,4:8:end);
    if followerInfo
        ampFollower = followerCoordAmp(:,4:8:end);
    end

    %mark merges and splits
    if markMS
        yaxisMin = min(ampSequence(:));
        yaxisMax = max(ampSequence(:));
        plot((splitTimes-axisOffset)*timeStep,(yaxisMin-0.1*(yaxisMax-yaxisMin))*ones(size(splitTimes)),'k+');
        plot((mergeTimes-axisOffset)*timeStep,(yaxisMin-0.1*(yaxisMax-yaxisMin))*ones(size(mergeTimes)),'rx');
        legend({'Splits','Merges'})
    end
    
    %plot the amplitudes, taking into account closed gaps, merges and
    %splits and sampling frequency doubling
    if followerInfo
        plotCompTrackCore(ampSequence,seqOfEvents,doubleFreq,timeStep,axisOffset,ampFollower,h);
    else
        plotCompTrackCore(ampSequence,seqOfEvents,doubleFreq,timeStep,axisOffset);
    end

    %put axes labels
    if ~inOneFigure
        if axisOffset
            xlabel('Time');
        else
            xlabel('Frame number');
        end
    end
    ylabel('Amplitude (a.u.)');

    %hold off of figure if each plot is in a separate figure window or if
    %this is the last plot
    if ~inOneFigure || (~plotAS)
        hold off
    end

end %(if plotA)

%aggregation state
if plotAS

    %open new figure window if needed and hold on to it
    if ~inOneFigure || (~plotX && ~plotY && ~plotA)
        figure, hold on
    end

    %if all figures will be in the same window, specify subplot
    if inOneFigure
        subplot(plotX+plotY+plotA+plotAS,1,plotX+plotY+plotA+1)
        hold on
    end

    %replace zeros with NaNs in aggregation state matrix
    aggregState(aggregState==0) = NaN;
    
    %mark merges and splits
    if markMS
        yaxisMin = min(aggregState(:));
        yaxisMax = max(aggregState(:));
        plot((splitTimes-axisOffset)*timeStep,(yaxisMin-0.1*(yaxisMax-yaxisMin))*ones(size(splitTimes)),'k+');
        plot((mergeTimes-axisOffset)*timeStep,(yaxisMin-0.1*(yaxisMax-yaxisMin))*ones(size(mergeTimes)),'rx');
        legend({'Splits','Merges'})
    end
    
    %plot the aggregation states, taking into account closed gaps, merges and
    %splits and sampling frequency doubling
    plotCompTrackCore(aggregState,seqOfEvents,doubleFreq,timeStep,axisOffset);

    %put axes labels
    if ~inOneFigure
        if axisOffset
            xlabel('Time');
        else
            xlabel('Frame number');
        end
    end
    ylabel('Aggregation state');

    %hold off of figure
    hold off
    
end %(if plotAS)

%put x-axis label at bottom in case of everything in one figure
if inOneFigure
    if axisOffset
        xlabel('Time');
    else
        xlabel('Frame number');
    end
end


%% Subfunction

function plotCompTrackCore(valuesMatrix,seqOfEvents,doubleFreq,timeStep,axisOffset,valuesFollower,h)

if nargin < 4 || isempty(timeStep)
    timeStep = 1;
end
if nargin < 5 || isempty(axisOffset)
    axisOffset = 0;
end
if nargin < 6
    valuesFollower = [];
    h = [];
end

%get first frame, last frame and number of frames
firstFrame = seqOfEvents(1,1);
lastFrame = seqOfEvents(end,1);

%get sampling frequency
samplingFreq = 1 / (1+doubleFreq);

%get number of segments making compound track
numSegments = size(valuesMatrix,1);

%plot values as dotted black lines, closing gaps
for i = 1 : numSegments
    indx = find(~isnan(valuesMatrix(i,:)));
    plot(((indx-1)*samplingFreq+firstFrame-axisOffset)*timeStep,valuesMatrix(i,indx),'k:');
end

%plot values in color, leaving gaps as blank (so that they appear as
%dotted lines in the final figure)
plot(((firstFrame:samplingFreq:lastFrame)'-axisOffset)*timeStep,valuesMatrix','marker','.');
if ~isempty(h)
    h.ColorOrderIndex=1;
end
if ~isempty(valuesFollower)
    plot(((firstFrame:samplingFreq:lastFrame)'-axisOffset)*timeStep,valuesFollower','marker','o');
end

%find merges and splits
indxSplit = (find(seqOfEvents(:,2) == 1 & ~isnan(seqOfEvents(:,4))))';
indxMerge = (find(seqOfEvents(:,2) == 2 & ~isnan(seqOfEvents(:,4))))';

%go over all splits
for iSplit = indxSplit
    
    %get time of splitting
    timeSplit = seqOfEvents(iSplit,1);
    
    %determine location in valuesMatrix
    splitLoc = (timeSplit - firstFrame) / samplingFreq + 1;
    
    %determine index of starting track
    rowS = seqOfEvents(iSplit,3);
    
    %determine index of splitting track
    rowSp = seqOfEvents(iSplit,4);
    
    %plot split as a black dash-dotted line
    plot(([timeSplit-samplingFreq timeSplit]-axisOffset)*timeStep,[valuesMatrix(rowSp,splitLoc-1) ...
        valuesMatrix(rowS,splitLoc)],'k-.')
    
end

%go over all merges
for iMerge = indxMerge
    
    %get time of merging
    timeMerge = seqOfEvents(iMerge,1);
    
    %determine location in valuesMatrix
    mergeLoc = (timeMerge - firstFrame) / samplingFreq + 1;
    
    %determine index of ending track
    rowE = seqOfEvents(iMerge,3);
    
    %determine index of merging track
    rowM = seqOfEvents(iMerge,4);
    
    %plot merge as a black dashed line
    plot(([timeMerge-samplingFreq timeMerge]-axisOffset)*timeStep,[valuesMatrix(rowE,mergeLoc-1) ...
        valuesMatrix(rowM,mergeLoc)],'k--')
    
end


%% ~~~ the end ~~~

