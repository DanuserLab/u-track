function [varargout] = imseriesshow(im, varargin ) 
% [varargout] = imseriesshow(im, varargin)
%
% A simple function to show series of 2D images from a volume.
%
%
%  Required Input Arguments: 
% 
%                 im - input grayscale image
% 
%  Optional arguments: param-name/value pairs
% 
%       displayColor - allows you to specify the RGB color in which to
%       			   dislpay the image should be a vector of length 3.
%					   Ex: [1,0,0] for red colore. 
%
%				       The image will be displayed in shades of the 
%                      specified color. Default is shades of gray.
% 
%           spacing - pixel spacing of the grayscale image
% 
%      displayRange - intensity display range of the input image
% 
%                     Specify this if you want to the grayscale image
%                     to be displayed in a specific intensity range. 
%                     All intensity values outside this range will be     
%                     made equal to the nearest border in the specified
%                     intensity range.                           
%
% Example:
%    
%     % 2D case
%     im2d = zeros(500,500);
%     im2d(100:400,100:400) = rand(301,301);
% 
%     imseriesshow( im2d ); % 2D image display
%     set( gcf, 'Name', '2D Example: Image display' );
% 
%     % 3D case
%     im3d = zeros(500,500,10);
%     im3d(100:400,100:400,:) = rand(301,301,10);
% 
%     imseriesshow( im3d ); % 3D image overlay
%     set( gcf, 'Name', '3D Example: 3D image display' ); 
% 
%     imseriesshow( im3d, 'displayColor', [0,1,0] ); % 3D colored image display
%     set( gcf, 'Name', '3D Example: 3D image display in a user specified color' ); 
%
% Author: Deepak Roy Chittajallu 
%
%% Parameter Parsing
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

% get required parameters
p = inputParser;
p.addRequired( 'im', @(x)( isnumeric(x) && ~isscalar(x) && ismember( ndims(x), [2,3] ) ) );
p.parse(im);

im = double(p.Results.im);
volSize = size(im);
volSize(3) = size(im, 3);

% get optional parameters
defaultDisplayRange = ComputeImageDynamicRange( im, 99.0 );
p.addParamValue( 'spacing', [1,1,1], @(x) ( isnumeric(x) && ~isscalar(x) && numel(x) == ndims(im) ) );
p.addParamValue( 'displayRange', defaultDisplayRange, @(x) (isempty(x) || (isnumeric(x) && numel(x) == 2)) );
p.addParamValue( 'displayColor', [1, 1, 1], @(x) ( isnumeric(x) && ~isscalar(x) && numel(x) == 3 ) );
p.parse(im, varargin{:});

spacing = p.Results.spacing;
displayrange = p.Results.displayRange;
displaycolor = p.Results.displayColor;

if isempty(displayrange)
    displayrange = [min(im(:)), max(im(:))];
end

%% Data
data.im = im;
data.volSize = volSize;
data.displayrange = displayrange;
data.displaycolor = displaycolor;

% Slice and planes
data.sliceno = 1;
data.curPlane = 3;
data.plane(1).slice = 1;
data.plane(2).slice = 1;
data.plane(3).slice = 1;

% Planes volume sizes
data.plane(1).volSize = [volSize(3) volSize(1) volSize(2)];
data.plane(2).volSize = [volSize(3) volSize(2) volSize(1)];
data.plane(3).volSize = [volSize(1) volSize(2) volSize(3)];

% Default x-limits and y-limits for each plane
data.plane(1).xlim = [ 0.5 volSize(1)+0.5 ];
data.plane(1).ylim = [ 0.5 volSize(3)+0.5 ];
data.plane(2).xlim = [ 0.5 volSize(2)+0.5 ];
data.plane(2).ylim = [ 0.5 volSize(3)+0.5 ];
data.plane(3).xlim = [ 0.5 volSize(2)+0.5 ];
data.plane(3).ylim = [ 0.5 volSize(1)+0.5 ];

% Data spacing
data.spacingUse = ~ismember( 'spacing', p.UsingDefaults );
data.plane(1).spacing = [spacing(3) spacing(1) 1];
data.plane(2).spacing = [spacing(3) spacing(2) 1];
data.plane(3).spacing = [spacing(2) spacing(1) 1];

% Data log
data.logUse = 0;
data.imLog = ComputeImageLogTransformForDisplay( data.im );
data.imLogDisplayRange = [ min(data.imLog(:)), max(data.imLog(:)) ];

%% Create UI controls
hMainFigure = figure;
clf(hMainFigure);

% Display figure toolbar
set(hMainFigure, 'ToolBar', 'figure')

% Default Axis
data.ui.hAxes = axes('Position', [0.0 , 0.18 , 1.0 , 0.8], ...
                     'YDir','reverse',...
                     'TickDir', 'out', ...
                     'XGrid', 'off', ...
                     'YGrid', 'off', ...
                     'DataAspectRatio', [1 1 1], ...
                     'PlotBoxAspectRatio', [volSize(2) volSize(1) 1], ...
                     'Parent', hMainFigure, ...
                     'Visible', 'off');

% Default image
data.ui.hImage = image('BusyAction', 'cancel', ...
           'Parent', data.ui.hAxes, ...
           'Interruptible', 'off');

% Slice navigation controls
data.ui.pbh_dec = uicontrol(hMainFigure, 'Style', 'pushbutton', 'String', '<<', ...
    'Units', 'normalized', 'Position', [0.3 0.11 0.066 0.05], ...
    'Callback', {@pushFirstSlice_Callback});

data.ui.pbh_dec = uicontrol(hMainFigure, 'Style', 'pushbutton', 'String', '<', ...
    'Units', 'normalized', 'Position', [0.366 0.11 0.066 0.05], ...
    'Callback', {@pushdec_Callback});

data.ui.eth_sno = uicontrol(hMainFigure, 'Style', 'edit', ...
    'String', '0',...
    'Units', 'normalized', 'Position', [0.433 0.11 0.133 0.05]);

data.ui.pbh_inc = uicontrol(hMainFigure, 'Style', 'pushbutton', 'String', '>', ...
    'Units', 'normalized', 'Position', [0.566 0.11 0.066 0.05], ...
    'Callback', {@pushinc_Callback});

data.ui.pbh_inc = uicontrol(hMainFigure, 'Style', 'pushbutton', 'String', '>>', ...
    'Units', 'normalized', 'Position', [0.632 0.11 0.066 0.05], ...
    'Callback', {@pushLastSlice_Callback});

% Cursor point info controls
data.ui.eth_xloc = uicontrol(hMainFigure, 'Style', 'edit', ...
    'String', 'X: INV', ...
    'Units', 'normalized', 'Position', [0.3 0.06 0.133 0.05]);

data.ui.eth_yloc = uicontrol(hMainFigure, 'Style', 'edit', ...
    'String', 'Y: INV', ...
    'Units', 'normalized', 'Position', [0.433 0.06 0.133 0.05]);

data.ui.eth_Imval = uicontrol(hMainFigure, 'Style', 'edit', ...
    'String', 'I: INV', ...
    'Units', 'normalized', 'Position', [0.566 0.06 0.133 0.05]);

% Planes ui
data.ui.ch_planes = uibuttongroup('visible', 'on', 'Units', 'normalized', ...
    'Position', [0.3 0.01 0.4 0.05], ...
    'SelectionChangeFcn', {@planes_Callback});

if ndims(data.im) == 3
    otherPlaneEnableState = 'on';
else
    otherPlaneEnableState = 'off';
end

data.ui.planeSP = uicontrol('Style', 'Radio', 'String', 'YZ',...
    'Units', 'normalized', 'Position', [0.05 0.15 0.20 0.75], ...
    'parent', data.ui.ch_planes, 'HandleVisibility', 'off','Enable',otherPlaneEnableState);
data.ui.planeCP = uicontrol('Style', 'Radio', 'String', 'XZ',...
    'Units', 'normalized', 'Position', [0.25 0.15 0.20 0.75], ...
    'parent', data.ui.ch_planes, 'HandleVisibility', 'off','Enable',otherPlaneEnableState);
data.ui.planeTP = uicontrol('Style', 'Radio', 'String', 'XY',...
    'Units', 'normalized', 'Position', [0.45 0.15 0.20 0.75], ...
    'parent', data.ui.ch_planes, 'HandleVisibility', 'off', 'Value', 1);
data.ui.planeLog = uicontrol('Style', 'checkbox', 'String', 'log',...
    'Units', 'normalized', 'Position', [0.65 0.15 0.20 0.75], ...
    'parent', data.ui.ch_planes, 'HandleVisibility', 'off', 'Value', data.logUse, ...
    'Callback',{@checkLog_Callback});
data.ui.planeSpacing = uicontrol('Style', 'checkbox', 'String', 'Spacing',...
    'Units', 'normalized', 'Position', [0.85 0.15 0.10 0.75], ...
    'parent', data.ui.ch_planes, 'HandleVisibility', 'off', 'Value', data.spacingUse, ...
    'Callback',{@checkSpacing_Callback});

% Set callbacks
set(hMainFigure, 'WindowScrollWheelFcn', @FnSliceScroll_Callback)
set(hMainFigure, 'WindowButtonMotionFcn', @FnMainFig_MouseMotionFunc);

% Pan and Zoom callback
hZoom = zoom(hMainFigure);
set(hZoom, 'ActionPostCallback', @postPanZoom_Callback);
hPan = pan(hMainFigure);
set(hPan, 'ActionPostCallback', @postPanZoom_Callback);

% Define default output and return it if it is requested by users
mOutputArgs{1} = hMainFigure;
if nargout > 0
    [varargout{1:nargout}] = mOutputArgs{:};
end
guidata(hMainFigure,data);
imsliceshow(data);
%close(gcbf);

end

%% Show slice
function imsliceshow(data)

% Compute rgb image and display it
cdata = getRGBSlice(data);
set(data.ui.hImage, 'cdata', cdata);

% Slice number
volSize = data.plane(data.curPlane).volSize;
strtmp = sprintf('%d / %d', data.sliceno, volSize(3));
set(data.ui.eth_sno, 'String', strtmp);

end

%% Compute a rgb slice
function [ slice ] = getRGBSlice(data)

if data.logUse      
    im = data.imLog;
    displayrange = data.imLogDisplayRange;    
else
    im = data.im;
    displayrange = data.displayrange;
end

% Get the slice, depending on the plane showed
if data.curPlane == 1
    slice = squeeze(im(:,data.sliceno,:))';
elseif data.curPlane == 2
    slice = squeeze(im(data.sliceno,:,:))';
else
    slice = squeeze(im(:,:,data.sliceno));
end    

slice = mat2gray(slice, displayrange);    

% Transform to gray
slice = repmat(slice, [1 1 3]);
color_hsv = rgb2hsv(data.displaycolor);
slice(:,:,1) = color_hsv(1);
slice(:,:,2) = color_hsv(2);
slice = hsv2rgb( slice );

end

%% Get real coordinates, depending on the plane
function [ ind coord ] = getRealCoord(data, coord)

if data.curPlane == 1
    coord = [coord(2) coord(3) coord(1)];
elseif data.curPlane == 2
    coord = [coord(3) coord(2) coord(1)];
end
ind = sub2ind(data.volSize, coord(1), coord(2), coord(3));

end

%% Zoom post callback
function postPanZoom_Callback(hSrc, eventdata) %#ok<INUSD>

data = guidata(hSrc);

% Save x-limits and y-limits
data.plane(data.curPlane).xlim = get(data.ui.hAxes, 'XLim');
data.plane(data.curPlane).ylim = get(data.ui.hAxes, 'YLim');

% Check the plot ratio is good
volSize = data.plane(data.curPlane).volSize;
set(data.ui.hAxes, 'PlotBoxAspectRatio', [volSize(2) volSize(1) 1]);

% Spacing
if data.spacingUse
    set(data.ui.hAxes, 'DataAspectRatio', data.plane(data.curPlane).spacing);
else
    set(data.ui.hAxes, 'DataAspectRatio', [1 1 1]);
end

guidata(hSrc, data);

end

%% Check/Uncheck Spacing use
function checkSpacing_Callback(hSrc, eventdata) %#ok<INUSD>

data = guidata(hSrc);

% Modify spacing
data.spacingUse = get(data.ui.planeSpacing, 'Value');
if data.spacingUse
    set(data.ui.hAxes, 'DataAspectRatio', data.plane(data.curPlane).spacing);
else
    set(data.ui.hAxes, 'DataAspectRatio', [1 1 1]);
end

guidata(hSrc, data);
imsliceshow(data);

end

%% Check/Uncheck Log use
function checkLog_Callback(hSrc, eventdata) %#ok<INUSD>

    data = guidata(hSrc);

    % Modify spacing
    data.logUse = get(data.ui.planeLog, 'Value');

    guidata(hSrc, data);
    imsliceshow(data);

end

%% Planes
function planes_Callback(hSrc, eventdata) %#ok<INUSD>

data = guidata(hSrc);

% Save slice number, x-limits and y-limits
data.plane(data.curPlane).slice = data.sliceno;
data.plane(data.curPlane).xlim = get(data.ui.hAxes, 'XLim');
data.plane(data.curPlane).ylim = get(data.ui.hAxes, 'YLim');

% Change to new plane
if get(data.ui.planeSP, 'Value')        % Sagittal
    data.curPlane = 1;
elseif get(data.ui.planeCP, 'Value')    % Coronal
    data.curPlane = 2;
elseif get(data.ui.planeTP, 'Value')    % Tranverse
    data.curPlane = 3;
end
data.sliceno = data.plane(data.curPlane).slice;

% Plot ratio
volSize = data.plane(data.curPlane).volSize;
set(data.ui.hAxes, 'PlotBoxAspectRatio', [volSize(2) volSize(1) 1]);

% Set zoom
set(data.ui.hAxes, 'XLim', data.plane(data.curPlane).xlim);
set(data.ui.hAxes, 'YLim', data.plane(data.curPlane).ylim);

% Spacing
if data.spacingUse
    set(data.ui.hAxes, 'DataAspectRatio', data.plane(data.curPlane).spacing);
else
    set(data.ui.hAxes, 'DataAspectRatio', [1 1 1]);
end

% Save and show slice
guidata(hSrc, data);
imsliceshow(data);

end

%% First Slice
function pushFirstSlice_Callback(hSrc, eventdata) %#ok<INUSD>

data = guidata(hSrc);
data.sliceno = 1;

guidata(hSrc, data);
imsliceshow(data);

end

%% Last Slice
function pushLastSlice_Callback(hSrc, eventdata) %#ok<INUSD>

data = guidata(hSrc);
volSize = data.plane(data.curPlane).volSize;
data.sliceno = volSize(3);

guidata(hSrc, data);
imsliceshow(data);

end

%% Dec Slice
function pushdec_Callback(hSrc, eventdata) %#ok<INUSD>

data = guidata(hSrc);

if data.sliceno > 1
    data.sliceno = data.sliceno - 1;
end

guidata(hSrc, data);
imsliceshow(data);

end

%% Inc Slice
function pushinc_Callback(hSrc, eventdata) %#ok<INUSD>

data = guidata(hSrc);

volSize = data.plane(data.curPlane).volSize;
if data.sliceno < volSize(3)
    data.sliceno = data.sliceno+1;
end

guidata(hSrc, data);
imsliceshow(data);

end

%% Slice scrolling
function FnSliceScroll_Callback(hSrc, eventdata)

data = guidata(hSrc);
volSize = data.plane(data.curPlane).volSize;

if eventdata.VerticalScrollCount > 0
    if data.sliceno < volSize(3)
        data.sliceno = data.sliceno + 1;
    end
elseif eventdata.VerticalScrollCount < 0
    if data.sliceno > 1 
        data.sliceno = data.sliceno - 1;
    end
end

guidata(hSrc, data);
imsliceshow(data);
UpdateCursorPointInfo(data);

end

%% Update cursor point info -- xloc, yloc, int_val
function UpdateCursorPointInfo(data)

cp = get(gca, 'CurrentPoint');
cp = round(cp(1, 1:2));
volSize = data.plane(data.curPlane).volSize;

% Display pointer coordinates and value
if IsPointInsideImage( cp(1,1:2) , data )
    set(data.ui.eth_xloc, 'String', sprintf('X: %d / %d', cp(1,1), volSize(2)));
    set(data.ui.eth_yloc, 'String', sprintf('Y: %d / %d', cp(1,2), volSize(1)));
    
    ind = getRealCoord(data, [cp(2) cp(1) data.sliceno]);
    set(data.ui.eth_Imval,'String',sprintf('I: %.1f' , data.im(ind)));
else
    set(data.ui.eth_xloc,'String',sprintf('X: INV') );
    set(data.ui.eth_yloc,'String',sprintf('Y: INV') );
    set(data.ui.eth_Imval,'String',sprintf('I: INV') );
end

end

%% Figure pointer
function FnMainFig_MouseMotionFunc(hSrc, eventdata) %#ok<INUSD>

% No change if use zoom or pan
hZoom = zoom(hSrc);
hPan = pan(hSrc);
data = guidata(hSrc);

if strcmp(get(hZoom, 'Enable'), 'off') && strcmp(get(hPan, 'Enable'), 'off')
    cp = get(gca, 'CurrentPoint');
    if IsPointInsideImage(cp(1,1:2), data)
        set(hSrc, 'Pointer', 'crosshair');
    else
        set(hSrc, 'Pointer', 'arrow');
    end
end

% Cursor info
UpdateCursorPointInfo(data);

end

%%
function [ blnInside ] = IsPointInsideImage(cp, data)

% Point inside figure limits
volInfLim = ceil([ data.plane(data.curPlane).xlim(1) data.plane(data.curPlane).ylim(1) ]);
volSupLim = floor([ data.plane(data.curPlane).xlim(2) data.plane(data.curPlane).ylim(2) ]);
blnInside = all( cp <= volSupLim ) && all( cp >= volInfLim );

end

%%
function [ imCorrected ] = AdjustImageIntensityRange( im, ImageIntensityRange )
    
    imCorrected = mat2gray( im, ImageIntensityRange ) * range(ImageIntensityRange);
    
end

%%
function [ intensityRange ] = ComputeImageDynamicRange( im, cover_percent )

    [p,x] = hist( im, 255 );   
    p = p / sum(p);
    pcum = cumsum(p);
    
    min_xlow = [];
    min_xhigh = [];
    min_xdiff = [];
    
    for i = 1:numel(x)
        for j = i+1:numel(x)
    
            if (pcum(j) - pcum(i) + p(i)) < 0.01 * cover_percent
                continue;
            end
            
            if isempty(min_xdiff) || (x(j) - x(i)) < min_xdiff
                min_xlow = x(i);
                min_xhigh = x(j);
                min_xdiff = x(j) - x(i);              
            end
        end
    end
    
    w = 0.5 * (x(2) - x(1));
    intensityRange = [min_xlow-w, min_xhigh+w];
    
end

%%
function [ imLog ] = ComputeImageLogTransformForDisplay( im )

    imLog = im - min( im(:) );
    ImageIntensityRange = ComputeImageDynamicRange( imLog, 99.0 );
    log_bottom = ImageIntensityRange(1) + range(ImageIntensityRange)/256.0 + eps; % just to give log a bottom
    imLog = log_bottom + AdjustImageIntensityRange( imLog, ImageIntensityRange );
    imLog = log( imLog );
    
end