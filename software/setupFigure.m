%[ha, hi, hf] = setupFigure(varargin) generates a multi-panel figure
%
% Optional inputs (first, second, third arguments):
%      nh: number of rows
%      nw: number of columns
%      na: total number of axes
%
% Options (specifier, value pairs):
%     SameAxes: true|{false} omits unnecessary tick labels if the data range
%               is the same in all panels
%    AxesWidth: width of each panel, in cm
%   AxesHeight: height of each panel, in cm
%
% Example:
% setupFigure(2,2,3) generates 3 panels in a 2x2 arrangement
%
% Francois Aguet, 2013
% Andrew R. Jamieson, 2016 - updated axes arrays to gobjects
%
% Copyright (C) 2019, Danuser Lab - UTSouthwestern 
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

function [ha, hi, hf] = setupFigure(varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addOptional('nh', 1, @isposint);
ip.addOptional('nw', 1, @isposint);
ip.addOptional('na', [], @isposint);
ip.addParameter('SameAxes', false, @islogical);
ip.addParameter('AspectRatio', []);
ip.addParameter('AxesWidth', []);
ip.addParameter('AxesHeight', []);
ip.addParameter('XSpace', []);
ip.addParameter('YSpace', []);
ip.addParameter('DisplayMode', 'print',  ... 
        @(x) any(strcmpi(x, {   'print', 'screen','multiple', ...
                                'u-track-3D-paper', 'JCB-Gerlich-paper','JCB-Gerlich-WIPS'})));
ip.addParameter('InsetPosition', []);
ip.addParameter('Name', '');
ip.addParameter('Box', 'off', @(x) any(strcmpi(x, {'on', 'off'})));
ip.parse(varargin{:});
nh = ip.Results.nh;
nw = ip.Results.nw;
na = ip.Results.na;
if isempty(na) || na>nh*nw
    na = nh*nw;
end

aw0 = ip.Results.AxesWidth;
ah0 = ip.Results.AxesHeight;
XSpace = ip.Results.XSpace;
YSpace = ip.Results.YSpace;

switch ip.Results.DisplayMode
    case 'print'
        if isempty(aw0)
            aw0 = 6;
        end
        if isempty(ah0)
            ah0 = 3.5;
        end
        if isempty(XSpace)
            XSpace = [1.5 0.75 0.5];
        end
        if isempty(YSpace)
            YSpace = [1.5 0.75 0.5];
        end
        axesFont = {'FontName', 'Helvetica', 'FontSize', 10};
    case 'screen'
        if isempty(aw0)
            aw0 = 12;
        end
        if isempty(ah0)
            ah0 = 7;
        end
        if isempty(XSpace)
            XSpace = [2 1.5 1];
        end
        if isempty(YSpace)
            YSpace = [2 1.5 1];
        end
        axesFont = {'FontName', 'Helvetica', 'FontSize', 12};
    case 'multiple'
        if isempty(aw0)
            aw0 = 12;
        end
        if isempty(ah0)
            ah0 = 7;
        end
        if isempty(XSpace)
            XSpace = [6 2.5 1];
        end
        if isempty(YSpace)
            YSpace = [6 2.5 1];
        end
        axesFont = {'FontName', 'Helvetica', 'FontSize', 8};
    case 'JCB-Gerlich-paper'
        if isempty(aw0)
            aw0 = 6;
        end
        if isempty(ah0)
            ah0 = 3.5;
        end
        if isempty(XSpace)
            XSpace = [2 2 1];
        end
        if isempty(YSpace)
            YSpace = [1.2 1 0.5]; % Bottom center Top
        end
        axesFont = {'FontName', 'Arial', 'FontSize', 8};
    case 'u-track-3D-paper'
        if isempty(aw0)
            aw0 = 6;
        end
        if isempty(ah0)
            ah0 = 3.5;
        end
        if isempty(XSpace)
            XSpace = [2 2 1];
        end
        if isempty(YSpace)
            YSpace = [1.2 1 0.5]; % Bottom center Top
        end
        axesFont = {'FontName', 'Helvetica', 'FontSize', 6};
   case 'JCB-Gerlich-WIPS'
        if isempty(aw0)
            aw0 = 6;
        end
        if isempty(ah0)
            ah0 = 3.5;
        end
        if isempty(XSpace)
            XSpace = [6 2 2];
        end
        if isempty(YSpace)
            YSpace = [4 2 2]; % Bottom center Top
        end
        axesFont = {'FontName', 'Arial', 'FontSize', 6};
end
tickLength = [0.015 0.025]*6/aw0;


w0 = aw0 + sum(XSpace);

h0 = ah0 + sum(YSpace);
if ~isempty(ip.Results.AspectRatio)
    ah0 = ip.Results.AspectRatio*aw0;
end

% default proportions: left/bottom: 1.5, width: 6, height: 3.5, top/right: 0.5
aw = aw0/w0;
xl = XSpace(1)/w0; % left spacing (relative to single axes)
if ~ip.Results.SameAxes && isempty(ip.Results.XSpace)
    xc = xl;
else
    xc = XSpace(2)/w0;  % spacing btw axes
end
xr = XSpace(3)/w0;

ah = ah0/h0;
yb = YSpace(1)/h0;
if ip.Results.SameAxes
    yc = YSpace(2)/h0;
else
    yc = yb;
end
yt = YSpace(3)/h0;

% width (relative to normalized single axes)
w = xl + nw*aw + (nw-1)*xc + xr;
% height
h = yb + nh*ah + (nh-1)*yc + yt;


% convert to normalized units
aw = aw/w;
xl = xl/w;
xc = xc/w;

ah = ah/h;
yb = yb/h;
yc = yc/h;

% resize figure window
fpos = get(0, 'DefaultFigurePosition')/get(0,'ScreenPixelsPerInch')*2.54;
fpos(3) = w*w0;
fpos(4) = h*h0;

fposPx = fpos/2.54*get(0,'ScreenPixelsPerInch');
fpos0 = get(0, 'ScreenSize');
c = min(0.8*fpos0(3:4)./fposPx(3:4));
if c<1
    fpos(3:4) = c*fpos(3:4);
end
fpos(2) = 5;
hf = figure('PaperPositionMode', 'auto', 'Color', 'w', 'InvertHardcopy', 'off',...
    'Units', 'centimeters', 'Position', fpos, 'Units', 'pixels', 'Name', ip.Results.Name);

ha = gobjects(na,1);
x0 = zeros(na,1);
y0 = zeros(na,1);
hi = gobjects(na,1);
ipos = ip.Results.InsetPosition;
if numel(ipos)==2
    %ipos = [ipos 0.95-ipos(1) 0.95-ipos(2)];
    ipos = [ipos 1-ipos(1) 1-ipos(2)];
end
for i = 1:na
    y0(i) = nh-ceil(i/nw);
    x0(i) = mod(i-1,nw);
    pos = [xl+x0(i)*(aw+xc) yb+y0(i)*(ah+yc) aw ah];
    ha(i) = axes('Position', pos);
    hold(ha(i), 'on');
    
    if ~isempty(ipos)
        hi(i) = axes('Position', [pos(1)+pos(3)*ipos(3) pos(2)+pos(4)*ipos(4)...
            ipos(1)*pos(3) ipos(2)*pos(4)]);
        hold(hi(i), 'on');
    end
end

if ip.Results.SameAxes
%     hayy = ha(y0>0);
%     haxx = ha(y0>0);
    set(ha(y0>0), 'XTickLabel', []);
    set(ha(x0>0), 'YTickLabel', []);
end

set(ha, 'TickDir', 'out', 'TickLength', tickLength,...
    'LineWidth', 1, 'Layer', 'top', axesFont{:});

if strcmpi(ip.Results.Box, 'on') % for a closed box w/o ticks, re-plot axes on top
    for i = 1:na
        hb = axes('Position', get(ha(i), 'Position'), 'Box', 'on', 'XTick', [],...
            'YTick', [], 'Color', 'none', 'LineWidth', 1);
        linkaxes([ha(i) hb]);
    end
end

if ~isempty(ipos)
    set(hi, 'TickDir', 'out', 'TickLength', tickLength*max(aw,ah)/max(ipos(1)*pos(3), ipos(2)*pos(4)),...
    'LineWidth', 1, 'Layer', 'top');
end
