classdef TracksDisplay < MovieDataDisplay
    %Conrete class for displaying flow
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
    properties
        Linestyle='-';
        Linewidth=1;
        GapLinestyle='--';
        Color = 'r';
        MergeColor = 'y';
        MergeMarker = 's';
        SplitColor = 'w';
        SplitMarker = 's';
        useDragtail=true;
        dragtailLength=10;
        showLabel=false;
        markMergeSplit=false;
        ButtonDownFcn=[];
    end
    methods
        function obj=TracksDisplay(varargin)
            nVarargin = numel(varargin);
            if nVarargin > 1 && mod(nVarargin,2)==0
                for i=1 : 2 : nVarargin-1
                    obj.(varargin{i}) = varargin{i+1};
                end
            end
        end
        function h=initDraw(obj, tracks, tag, varargin)
            if(~ishandle(varargin{1}))
                hIn = [];
            else
                hIn = varargin{1};
                varargin = varargin(2:end);
            end
            if isempty(tracks)
                set(hIn,'Visible','off');
                h = -1;
                return;
            end
            % Get track length and filter valid tracks
            trackLengths = cellfun('prodofsize',{tracks.xCoord});
            validTracks = trackLengths>0;
            tracks = tracks(validTracks);
            trackLengths = trackLengths(validTracks);
            
            nTracks = numel(tracks);
            
            % Constraing the dragtail length between 2 and the maximum
            % track length
            if obj.useDragtail
                dLength = max(2,min(obj.dragtailLength,max(trackLengths)));
            else
                dLength = max(trackLengths);
            end
            
            % Concatenate data in a matrix of size dragtailLength x nTracks
            xData = NaN(dLength, nTracks);
            yData = NaN(dLength, nTracks);
            uTrackLengths = unique(trackLengths);
            for i = uTrackLengths(:)'
                selected = trackLengths == i;
                xTemp = vertcat(tracks(selected).xCoord)';
                yTemp = vertcat(tracks(selected).yCoord)';
                if(i < dLength)
                    xData(1:i,selected) = xTemp;
                    yData(1:i,selected) = yTemp;
                else
                    xData(:,selected) = xTemp(end-dLength+1:end,:);
                    yData(:,selected) = yTemp(end-dLength+1:end,:);
                end
            end
            displayLength = trackLengths;
            displayLength(trackLengths > dLength) = dLength;
            
            % Initialize matrix for gaps
            xGapData = NaN(size(xData));
            yGapData = NaN(size(xData));
            
            % Label gaps: series of NaNs not-connected to the border
            I = isnan(xData);
            I = [I; zeros(size(I))];
            I = reshape(I, size(I,1)/2, size(I,2)*2);
            I = [zeros(size(I,1), 1) I];
            I = imclearborder(I);
            
            cc = bwconncomp(I);
            gapLengths = cellfun('length',cc.PixelIdxList);

            % Obtain first and last coordinate indices
            iFirstLast = zeros(2,cc.NumObjects);
            uGapLengths = unique(gapLengths);
            for i=uGapLengths(:)'
                s = i == gapLengths;
                temp = [cc.PixelIdxList{s}];
                iFirstLast(:,s) = temp([1 i],:);
            end

            % Remove extra columns by adjusting index
            m = mod(iFirstLast,cc.ImageSize(1));
            iFirstLast = (iFirstLast + m - cc.ImageSize(1))/2;
            
            
            iFirst = iFirstLast(1,:) - 1;
            iLast = iFirstLast(2,:) + 1;
            xFirst = xData(iFirst);
            xLast = xData(iLast);
            yFirst = yData(iFirst);
            yLast = yData(iLast);
            
            % Fill gaps
            for gapLength = uGapLengths(:)'
                s = gapLengths == gapLength;
                px = [cc.PixelIdxList{s}];
                
                % Remove extra columns by adjusting index
                m = mod(px,cc.ImageSize(1));
                px = (px + m - cc.ImageSize(1))/2;
                
                a = (1:gapLength)'/(gapLength+1);
                xGapData(iFirst) = xFirst;
                xGapData(iLast) = xLast;
                yGapData(iFirst) = yFirst;
                yGapData(iLast) = yLast;
                xGapData(px) = ...
                    bsxfun(@plus,bsxfun(@times,a,xLast(s) - xFirst(s)),xFirst(s));
                yGapData(px) = ...
                    bsxfun(@plus,bsxfun(@times,a,yLast(s) - yFirst(s)),yFirst(s));
            end
            
            eventsExist = isfield(tracks,'splitEvents') || isfield(tracks,'mergeEvents');
            dragtailWindows = [trackLengths - displayLength + 1 ; trackLengths];
            
            %% Initialize matrix for split events
            xSplitData = NaN(dLength, nTracks);
            ySplitData = NaN(dLength, nTracks);
            
            if(eventsExist)
                hasSplitEvents = ~cellfun('isempty',{tracks.splitEvents})';
                eventTimes = [tracks(hasSplitEvents).splitEvents];
                eventTracks = zeros(size(eventTimes));
                eventTrackIdx = cumsum([1 cellfun('length',{tracks(hasSplitEvents).splitEvents})]);
                hasSplitEventsIdx = find(hasSplitEvents);
            else
                hasSplitEventsIdx = [];
            end


            if(~isempty(hasSplitEventsIdx))
                eventTracks(eventTrackIdx(1:end-1)) = [hasSplitEventsIdx(1);  diff(hasSplitEventsIdx)];
                eventTracks = cumsum(eventTracks);

                eventTimes = [eventTimes eventTimes+1];
                eventTracks = [eventTracks eventTracks];
                f = eventTimes >= dragtailWindows(1,eventTracks) & eventTimes <= dragtailWindows(2,eventTracks);
                eventTimes = eventTimes(f);
                eventTracks = eventTracks(f);
                idx = sub2ind(size(xSplitData),eventTimes - dragtailWindows(1,eventTracks) + 1,eventTracks);
                xSplitData(idx) = xData(idx);
                ySplitData(idx) = yData(idx);
            end
                        
            %% Initialize matrix for split events
            xMergeData = NaN(dLength, nTracks);
            yMergeData = NaN(dLength, nTracks);

            
            if(eventsExist)
                hasMergeEvents = ~cellfun('isempty',{tracks.mergeEvents})';
                eventTimes = [tracks(hasMergeEvents).mergeEvents];
                eventTracks = zeros(size(eventTimes));
                eventTrackIdx = cumsum([1 cellfun('length',{tracks(hasMergeEvents).mergeEvents})]);
                hasMergeEventsIdx = find(hasMergeEvents);
            else
                hasMergeEventsIdx = [];
            end

            if(~isempty(hasMergeEventsIdx))
                eventTracks(eventTrackIdx(1:end-1)) = [hasMergeEventsIdx(1);  diff(hasMergeEventsIdx)];
                eventTracks = cumsum(eventTracks);

                eventTimes = [eventTimes-1 eventTimes];
                eventTracks = [eventTracks eventTracks];
                f = eventTimes >= dragtailWindows(1,eventTracks) & eventTimes <= dragtailWindows(2,eventTracks);
                eventTimes = eventTimes(f);
                eventTracks = eventTracks(f);
                idx = sub2ind(size(xMergeData),eventTimes - dragtailWindows(1,eventTracks) + 1,eventTracks);
                xMergeData(idx) = xData(idx);
                yMergeData(idx) = yData(idx);
            end
                       
            %% Plot tracks
            hasLabels = isfield(tracks,'label');
            hlinesIn = findobj(hIn,'Type','Line');
            if hasLabels % If track is classified
                nColors = size(obj.Color,1);
                h = -ones(nColors+1,2);
                % Attempt to reuse line objects, delete the rest
                idx = 1:min(numel(hlinesIn),numel(h));
                h(idx) = double(hlinesIn(idx));
                % Delete excess lines we do not need
                delete(hlinesIn(numel(h)+1:end));
                for iColor = 1:nColors
                    iTracks = mod([tracks.label]-1, nColors) +1 == iColor;
                    h(iColor,1)=plotFast(h(iColor,1),xData(:,iTracks),yData(:,iTracks),'Linestyle',obj.Linestyle,...
                        'Linewidth', obj.Linewidth, 'Color',obj.Color(iColor,:),'Marker','none',varargin{:});
                    h(iColor,2)=plotFast(h(iColor,2),xGapData(:,iTracks),yGapData(:,iTracks),'Linestyle',obj.GapLinestyle,...
                        'Linewidth', obj.Linewidth, 'Color',obj.Color(iColor,:),'Marker','none',varargin{:});
                end
                
                % Plot merge / splits on top of classied tracks
                splitMarker = 'none';
                mergeMarker = 'none';
                if(obj.markMergeSplit)
                    splitMarker = obj.SplitMarker;
                    mergeMarker = obj.MergeMarker;
                end
                h(nColors+1,1) = plotFast(h(nColors+1),xSplitData, ySplitData, 'Linestyle', obj.Linestyle,...
                    'Linewidth', obj.Linewidth, 'Color', obj.SplitColor , 'Marker', splitMarker , varargin{:});
                h(nColors+1,2) = plotFast(h(nColors+2),xMergeData, yMergeData, 'Linestyle', obj.Linestyle,...
                    'Linewidth', obj.Linewidth, 'Color', obj.MergeColor, 'Marker', mergeMarker , varargin{:});
            else
                % Plot links and gaps
                h=-ones(4,1);
                % Attempt to reuse line objects, delete the rest
                idx = 1:min(numel(hlinesIn),numel(h));
                h(idx) = double(hlinesIn(idx));
                delete(hlinesIn(numel(h)+1:end));
                splitMarker = 'none';
                mergeMarker = 'none';
                if(obj.markMergeSplit)
                    splitMarker = obj.SplitMarker;
                    mergeMarker = obj.MergeMarker;
                end
                h(1) = plotFast(h(1),xData, yData, 'Linestyle', obj.Linestyle,...
                    'Linewidth', obj.Linewidth, 'Color',obj.Color, 'Marker', 'none', varargin{:});
                h(2) = plotFast(h(2),xGapData, yGapData, 'Linestyle', obj.GapLinestyle',...
                    'Linewidth', obj.Linewidth, 'Color',[1 1 1] - obj.Color, 'Marker', 'none', varargin{:});
                h(3) = plotFast(h(3),xSplitData, ySplitData, 'Linestyle', obj.Linestyle,...
                    'Linewidth', obj.Linewidth, 'Color', obj.SplitColor , 'Marker', splitMarker , varargin{:});
                h(4) = plotFast(h(4),xMergeData, yMergeData, 'Linestyle', obj.Linestyle,...
                    'Linewidth', obj.Linewidth, 'Color', obj.MergeColor, 'Marker', mergeMarker , varargin{:});
                uistack(h(2:4),'top');
            end
            
            % Display track numbers if option is selected
            hTextIn = findobj(hIn,'Type','Text');
            if obj.showLabel
                % Convert track numbers to text
                trackNr = sprintf('%.0f\n',tracks.number);
                % Produces an extra string, but we ignore it
                [trackNr,~] = regexp(trackNr,'\n','split','match');
                % Offset the text by two pixels in x and y
                xDataOffset = xData + 2;
                yDataOffset = yData + 2;
                % Find last non-NaN coordinate
                isDataNaN = isnan(xData);
                lastIdx = zeros(1,size(xData,2));
                [r,c] = find(~isDataNaN);
                lastIdx(c) = r;
                % Like sub2ind
                lastLinIdx = lastIdx + (0:size(xDataOffset,2)-1)*size(xDataOffset,1);

                % Filter by those tracks present
                f = lastLinIdx ~= 0;
                lastLinIdx = lastLinIdx(f);
                trackNr = trackNr(f);
                if(numel(trackNr) > numel(hTextIn))
                    hlabels = hTextIn;
                    idx = 1:numel(hTextIn);
                    % Reposition existing labels
                    set(hlabels(idx),{'Position','String'}, ...
                        [num2cell( ...
                            [xDataOffset(lastLinIdx(idx)) ; ...
                            yDataOffset(lastLinIdx(idx))]', ...
                            2) ...
                         trackNr(idx)'] ...
                        );
                    set(hlabels(idx),'Visible','on');
                    % Create new labels
                    idx = numel(hTextIn)+1:numel(trackNr);
                    hlabels(idx) = text( ...
                        xDataOffset(lastLinIdx(idx)), ...
                        yDataOffset(lastLinIdx(idx)), ...
                        trackNr(idx), ...
                        'Clipping','on', ...
                        'HitTest','off' ...
                    );
                else
                    hlabels = hTextIn;
                    % Reposition existing labels
                    idx = 1:numel(trackNr);
                    set(hlabels(idx),{'Position','String'}, ...
                        [num2cell( ...
                            [xDataOffset(lastLinIdx(idx)) ; ...
                            yDataOffset(lastLinIdx(idx))]', ...
                            2) ...
                         trackNr(idx)'] ...
                        );
                    set(hlabels(idx),'Visible','on');
                    % Make unsused labels invisible
                    set(hlabels(numel(trackNr)+1:end),'Visible','off');
                end

                if hasLabels
                    iColors = mod([tracks(f).label]-1, nColors) + 1;
                    uiColors = unique(iColors);
                    % Set once per color
                    for iColor = uiColors
                        s = iColors == iColor;
                        set(hlabels(s),'Color',obj.Color(iColor,:));
                    end
                else
                    set(hlabels,'Color',obj.Color);
                end
                h = [h(:) ; double(hlabels(:)) ];
            else
                delete(hTextIn);
            end
            
            % Set tag
            set(h(ishandle(h)), 'Tag', tag, 'ButtonDownFcn', obj.ButtonDownFcn);
        end
        
        function updateDraw(obj, h, data)
            tag=get(h(1),'Tag');
            obj.initDraw(data,tag,h);
            return;
            
        end
    end
    
    methods (Static)
        function params=getParamValidators()
            params(1).name='Color';
            params(1).validator=@(x)ischar(x) ||isvector(x);
            params(2).name='Linestyle';
            params(2).validator=@ischar;
            params(3).name='GapLinestyle';
            params(3).validator=@ischar;
            params(4).name='dragtailLength';
            params(4).validator=@isscalar;
            params(5).name='showLabel';
            params(5).validator=@isscalar;
            params(6).name='useDragtail';
            params(6).validator=@islogical;
            params(7).name='MergeColor';
            params(7).validator=@(x)ischar(x) ||isvector(x);
            params(8).name='SplitColor';
            params(8).validator=@(x)ischar(x) ||isvector(x);
            params(9).name='ButtonDownFcn';
            params(9).validator=@(x) isempty(x) || isa(x, 'function_handle');
            params(10).name='markMergeSplit';
            params(10).validator=@isscalar;
        end
        
        function f=getDataValidator()
            f=@isstruct;
        end
    end
end
function h = plotFast(h,xData,yData,varargin)
%plotFast Uses the low-level plotting function line to plot matrices by
%separate lines into a single line with gaps while also trying to avoid
%creating a new object.
%
% INPUT
% h - handle to update if valid
% xData - same as given to plot, can be a matrix where columns are separate
% lines
% yData - same as given to plot, can be a matrix where columns are separate
% lines
% varargin - passes through extra parameters to the line function
%
% OUTPUT
%
% h - handle to the line object created or used
%

% Mark Kittisopikul, 2014/11/24

    xData(end+1,:) = NaN;
    yData(end+1,:) = NaN;
    if(ishandle(h))
        set(h,'XData',xData(:),'YData',yData(:),varargin{:});
    else
        h = line(xData(:),yData(:),varargin{:});
        if(isempty(h))
            h = -1;
        end
    end
end
