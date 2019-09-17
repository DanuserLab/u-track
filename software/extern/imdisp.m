%IMDISP  Display one or more images nicely
%
% Examples:
%   imdisp
%   imdisp(I)
%   imdisp(I, map)
%   imdisp(I, lims)
%   imdisp(I, map, lims)
%   imdisp(..., param1, value1, param2, value2, ...)
%   h = imdisp(...)
%
% This function displays one or more images nicely. Images can be defined
% by arrays, filenames or imstream objects (which can be videos or image
% sequences). Multiple images can be input in a cell array or stacked along
% the fourth dimension (as well as in an imstream object), and are
% displayed as a grid of subplots (an improvement over MONTAGE). The size
% of grid is calculated or user defined. The figure size is set so that
% images are magnified by an integer value.
%
% If the image grid size is user defined, images not fitting in the grid
% can be scrolled through using the following key presses:
%    Up - Back a row.
%    Down - Forward a row.
%    Left - Back a page (or column if there is only one row).
%    Right - Forward a page (or column if there is only one row).
%    Shift - 2 x speed.
%    Ctrl - 4 x speed.
%    Shift + Ctrl - 8 x speed.
% Other keypresses available are:
%    'f' - Print out the index of the current frame.
%    'g' - Go to frame input by user. Values <= 0 offset from the end.
%    'q' - End figure interactivity.
%    Esc - Close the figure.
%
% This allows fast scrolling through a movie or image stack, e.g. 
%    imdisp(imstream('xylophone.mp4'), 'Size', 1)
% The function can be used as a visual DIR, e.g. 
%    imdisp()
% to display all images in the current directory on a grid, or 
%    imdisp({}, 'Size', 1)
% to scroll through them one at a time.
%
% IN:
%   I - MxNxCxP array of images, 1xP cell array or imstream object. C is 1
%       for indexed images or 3 for RGB images. P is the number of images.
%       If I is a cell array then each cell must contain an image. Images
%       can equally be defined by filenames. If I is an empty cell array
%       then all the images in the current directory are used. Default: {}.
%   map - Kx3 colormap to be used with indexed images. Default: gray(256).
%   lims - [LOW HIGH] display range for indexed images. Default: [min(I(:))
%          max(I(:))].
%   Optional parameters - name, value parameter pairs for the following:
%      'Size' - [H W] size of grid to display image on. If only H is given
%               then W = H. If either H or W is NaN then the number of rows
%               or columns is chosen such that all images fit. If both H
%               and W are NaN or the array is empty then the size of grid
%               is chosen to fit all images in as large as possible.
%               Default: [].
%      'Indices' - 1xL list of indices of images to display. Default: 1:P.
%      'Border' - [TB LR] borders to give each image top and bottom (TB)
%                 and left and right (LR), to space out images. Borders are
%                 normalized to the subplot size, i.e. TB = 0.01 gives a
%                 border 1% of the height of each subplot. If only TB is
%                 given, LR = TB. Default: 0.01.
%      'DisplayRange' - Same as lims input.
%      'Map' - Kx3 colormap or (additionally from above) name of MATLAB
%              colormap, for use with indexed images. Default: gray(256).
%      'FigureSize' - [W H] size of the figure used to set the montage 
%                     layout.
%      'FigureHandle' - 1x1 figure handle. The dimensions of this figure
%                       will be used to set the montage layout.
%
% OUT:
%   h - HxW array of handles to images.
%
%   See also IMAGE, IMAGESC, IMSHOW, MONTAGE.
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

% Copyright: Oliver Woodford 2010-2014

function hIm = imdisp(I, varargin)
% Parse inputs
[map, layout, gap, indices, lims, figSize] = parse_inputs(varargin);

if nargin == 0 || (iscell(I) && isempty(I))
    % Read in all the images in the directory
    I = get_im_names();
    if isempty(I)
        % No images found
        if nargout > 0
            hIm = [];
        end
        return
    end
end

% Check if input is filenames
if ischar(I)
    [x, y, c] = size(I);
    if (x > 1 && y > 1) || c > 1 
        I = num2cell(I, 2);
    else
        I = {I(:)'};
    end
end

% Get limits, etc.
cell_or_stream = false;
if isnumeric(I) || islogical(I)
    [y, x, c, n] = size(I);
    if isempty(lims)
        lims = min_max(I);
    elseif isequal(0, lims)
        lims = default_limits(I);
    elseif c == 3
        % Rescale
        if ~isfloat(I)
            I = single(I);
        end
        I = min(max((I - lims(1)) ./ (lims(2) - lims(1)), 0), 1);
        lims = [0 1];
    end
    if isfloat(I) && c == 3 && n > 1
        I = uint8(I * 256 - 0.5);
        lims = round(lims * 256 - 0.5);
    end
elseif is_cell_or_stream(I)
    cell_or_stream = true;
    n = numel(I);
    A = I{1};
    if ischar(A)
        % Read in the image (or images for multi-frame files)
        if n == 1
            cell_or_stream = false;
            I = imread_rgb_multi(A);
            if iscell(I)
                n = numel(I);
                A = I{1};
                [y, x, c] = size(A);
            else
                [y, x, c, n] = size(I);
                A = I;
            end
        else
            A = imread_rgb(A);
            I{1} = A;
            [y, x, c] = size(A);
        end
    else
        [y, x, c] = size(A);
    end
    % Assume all images are the same size and type as the first
    if isempty(lims) || isequal(0, lims)
        lims = default_limits(A);
    end
else
    error('I not of recognized type.');
end

% Select indexed images
if ~isequal(indices, -1)
    if iscell(I)
        I = I(indices);
        n = numel(I);
    else
        I = I(:,:,:,indices);
        n = size(I, 4);
    end
end

% Get the current figure
hFig = get(0, 'CurrentFigure');
if isempty(hFig)
    % Create a new figure
    hFig = figure();
end

if n < 2
    % Set the colormap
    set(hFig, 'Colormap', map);
end

% Display the image(s)
if n == 0
    hIm = display_image([], gca, [0 1]);
    
    if nargout == 0
        clear hIm % Avoid printing this out
    end
    return
elseif n == 1
    % IMSHOW mode
    % Display the single image
    hAx = gca;
    if cell_or_stream
        I = I{1};
    end
    hIm = display_image(I, hAx, lims);
    
    if nargout == 0
        clear hIm % Avoid printing this out
    end
    
    % Only resize image if it is alone in the figure
    if numel(findobj(get(hFig, 'Children'), 'Type', 'axes')) > 1
        return
    end
    % Could still be the first subplot - do another check
    axesPos = get(hAx, 'Position');
    newAxesPos = [gap(1) gap(end) 1-2*gap(1) 1-2*gap(end)];
    if max(abs(axesPos - get(hFig, 'DefaultAxesPosition'))) < 1e-15
        % Default position => not a subplot
        % Fill the window
        set(hAx, 'Units', 'normalized', 'Position', newAxesPos);
        axesPos = newAxesPos;
    end
    if ~isequal(axesPos, newAxesPos)
        % Figure not alone, so don't resize.
        return
    end
    layout = [1 1];
else
    % MONTAGE mode
    % Compute a good layout
    layout = choose_layout(n, y, x, layout, figSize);

    % Create a data structure to store the data in
    num = prod(layout);
    state.num = num * ceil(n / num);
    hIm = zeros(layout);
    hAx = zeros(layout);
    
    % Clear the figure
    if n > num
        hFig = clf(hFig, 'reset');
    else
        hFig = clf(hFig);
    end
    
    % Set the colormap
    set(hFig, 'Colormap', map);

    % Set the first lot of images
    index = mod(0:num-1, state.num) + 1;
    hw = 1 ./ layout;
    gap = gap ./ layout;
    dims = hw - 2 * gap;
    dims = dims([2 1]);
    for a = 1:layout(1)
        for b = 1:layout(2)
            c = index(b + (layout(1) - a) * layout(2));
            if c > n
                A = [];
            elseif cell_or_stream
                A = I{c};
                if ischar(A)
                    A = imread_rgb(A);
                    I{c} = A;
                end
            else
                A = I(:,:,:,c);
            end
            hAx(a,b) = axes('Position', [(b-1)*hw(2)+gap(2) (a-1)*hw(1)+gap(1) dims], 'Units', 'normalized');
            hIm(a,b) = display_image(A, hAx(a,b), lims);
        end
    end
    
    % Check if we need to be able to scroll through images
    if n > num
        % Intialize rest of data structure
        state.hIm = hIm;
        state.hAx = hAx;
        state.index = 1;
        state.layout = layout;
        state.lims = lims;
        state.n = n;
        state.I = I;
        state.is_cell_or_stream = cell_or_stream;
        % Set the callback for image navigation, and save the image data in the figure
        set(hFig, 'KeyPressFcn', @keypress_callback, 'Interruptible', 'off', 'BusyAction', 'cancel', 'UserData', state, 'Name', 'Image index: 1.');
    end
    
    % Flip hIm so it matches the layout
    hIm = hIm(end:-1:1,:);
    
    if nargout == 0
        clear hIm % Avoid printing this out
    end
end

if strcmp(get(hFig, 'WindowStyle'), 'docked') || x == 0 || y == 0
    % Figure is docked or image is empty, so can't resize
    return
end

% Set the figure size well
% Compute the image size
ImSz = layout([2 1]) .* [x y] ./ (1 - 2 * gap([end 1]));
    
% Get the size of the monitor we're on
figPosCur = get(hFig, 'Position');
% Monitor sizes
MonSz = get(0, 'MonitorPositions');
if ~ishg2(hFig)
    % Correct the size
    MonSz(:,3:4) = MonSz(:,3:4) - MonSz(:,1:2) + 1;
end
MonOn = size(MonSz, 1);
if MonOn > 1
    % Make the origin the top left corner of the primary monitor
    correction = 0;
    if ispc
        for a = 1:MonOn
            if isequal(MonSz(a,1:2), [1 1])
                correction = MonSz(a,4);
                break
            end
        end
    end
    % Determine which monitor the centre of the image is on
    figCenter = figPosCur(1:2) + figPosCur(3:4) / 2;
    figCenter = MonSz(:,1:2) - repmat(figCenter, [MonOn 1]);
    figCenter = [figCenter figCenter+MonSz(:,3:4)];
    MonOn = all(sign(figCenter) == repmat([-1 -1 1 1], [MonOn 1]), 2);
    MonOn(1) = MonOn(1) | ~any(MonOn);
    MonSz = MonSz(MonOn,:);
    % Correct the origin
    if correction
        MonSz(2) = correction - MonSz(4) - MonSz(2) + 2;
    end
end

% Check if the window is maximized
% This is a hack which may only work on Windows! No matter, though.
if isequal(MonSz([1 3]), figPosCur([1 3]))
    % Leave maximized
    return
end

% Compute the size to set the window
MaxSz = MonSz(3:4) - [20 120];
RescaleFactor = min(MaxSz ./ ImSz);
if RescaleFactor > 1
    % Integer scale for enlarging, but don't make too big
    MaxSz = min(MaxSz, [1200 800]);
    RescaleFactor = max(floor(min(MaxSz ./ ImSz)), 1);
end
figPosNew = ceil(ImSz * RescaleFactor);

% Don't move the figure if the size isn't changing
if isequal(figPosCur(3:4), figPosNew)
    return
end

% Keep the centre of the figure stationary
figPosNew = [floor(figPosCur(1:2)+(figPosCur(3:4)-figPosNew)/2) figPosNew];

% Ensure the figure is in bounds
figPosNew(1:2) = min(max(figPosNew(1:2), MonSz(1:2)+6), MonSz(1:2)+MonSz(3:4)-[6 101]-figPosNew(3:4));

% Set the figure size and position
set(hFig, 'Position', figPosNew);
end

%% Keypress callback
% The function which does all the display stuff
function keypress_callback(fig, event_data)
% Check what key was pressed and update the image index as necessary
switch event_data.Character
    case 28 % Left
        up = -1; % Back a page
    case 29 % Right
        up = 1; % Forward a page
    case 30 % Up
        up = -0.1; % Back a row
    case 31 % Down
        up = 0.1; % Forward a row
    case 'f'
        state = get(fig, 'UserData');
        fprintf('Image index: %d\n', state.index); % Print out the index
        return;
    case 'g'
        % Get the user to input a frame number
        up = input('Enter frame number to go to: ');
        state = get(fig, 'UserData');
        if up <= 0
            up = state.num + up;
        end
        up = up - state.index;
    case 'q'
        set(fig, 'KeyPressFcn', []); % Quit the widget
        return;
    case 27 % Escape
        close(fig); % Close the figure
        return;
    otherwise
        % Another key was pressed - ignore it
        return
end
if event_data.Character <= 31
    % Get the state data
    state = get(fig, 'UserData');
    % Use control and shift for faster scrolling
    if ~isempty(event_data.Modifier)
        up = up * (2 ^ (strcmpi(event_data.Modifier, {'shift', 'control'}) * [1; 2]));
    end
end
% Get the current index
index = state.index;
% Get number of images
n = prod(state.layout);
% Generate valid indices
if abs(up) < 1
    % Increment by row, or by 10 if only one image
    index = index + state.layout(2) * (up * (10 ^ (2 - (prod(state.layout) > 1)))) - 1;
else
    if state.layout(1) == 1
        % Increment by column
        index = index + up - 1;
    else
        % Increment by page
        index = index + n * up - 1;
    end
end
index = mod(index:index+n, state.num) + 1;
% Plot the images
figure(fig);
for a = 1:state.layout(1)
    for b = 1:state.layout(2)
        % Get the image
        c = index(b + (state.layout(1) - a) * state.layout(2));
        if c > state.n
            % Set the image data
            set(state.hIm(a,b), 'CData', []);
        elseif state.is_cell_or_stream
            A = state.I{c};
            if ischar(A)
                % Filename - read the image from disk
                A = imread_rgb(A);
                state.I{c} = A;
            end
            % Set the image data
            set(state.hIm(a,b), 'CData', rescale_rgb(A, state.lims));
            % Reset the axes limits
            if ~isempty(A)
                set(state.hAx(a,b), 'XLim', [0.5 size(A, 2)+0.5], 'YLim', [0.5 size(A, 1)+0.5]);
            end
        else
            % Set the image data
            set(state.hIm(a,b), 'CData', state.I(:,:,:,c));
        end
    end
end
drawnow;
% Save the current index
state.index = index(1);
set(fig, 'UserData', state, 'Name', sprintf('Image index: %d.', state.index));
end

%% Display the image
function hIm = display_image(A, hAx, lims)
if isempty(A)
    hIm = image(zeros(1, 1, 3));
    set(hIm, 'CData', []);
else
    hIm = image(rescale_rgb(A, lims));
end
set(hAx, 'Visible', 'off', 'DataAspectRatio', [1 1 1], 'CLim', lims);
try
    set(hAx, 'SortMethod', 'childorder');
catch
    % Support older versions
    set(hAx, 'DrawMode', 'fast');
end
set(get(hAx, 'XLabel'), 'Visible', 'on');
set(get(hAx, 'YLabel'), 'Visible', 'on');
set(get(hAx, 'Title'), 'Visible', 'on');
set(hIm, 'CDataMapping', 'scaled');
end

%% Choose a good layout for the images
function layout = choose_layout(n, y, x, layout, sz)
v = numel(layout);
N = isnan(layout);
if v == 0 || all(N)
    % Compute approximate layout
    sz = sz ./ [x y];
    layout = ceil(sz([2 1]) ./ sqrt(prod(sz) / n));
    % Remove superfluous rows or columns
    while 1
        switch ([prod(layout - [1 0]) prod(layout - [0 1])] >= n) * [2; 1]
            case 0
                break;
            case 1
                layout = layout - [0 1];
            case 2
                layout = layout - [1 0];
            case 3
                if min(sz .* (layout - [0 1])) > min(sz .* (layout - [1 0]))
                    layout = layout - [0 1];
                else
                    layout = layout - [1 0];
                end
        end
    end
elseif v == 1
    layout = layout([1 1]);
elseif any(N)
    layout(N) = ceil(n / layout(~N));
end
layout = reshape(layout, 1, 2);
end

%% Read image to uint8 rgb array
function A = imread_rgb(name)
try
    [A, map, alpha] = imread(name);
catch
    % Format not recognized by imread, so create a red cross (along diagonals)
    A = eye(101) | diag(ones(100, 1), 1) | diag(ones(100, 1), -1);
    A = (uint8(1) - uint8(A | flipud(A))) * uint8(255);
    A = cat(3, zeros(size(A), 'uint8')+uint8(255), A, A);
    return
end
A = A(:,:,:,1); % Keep only first frame of multi-frame files
if ~isempty(map)
    map = uint8(map * 256 - 0.5); % Convert to uint8 for storage
    A = reshape(map(uint32(A)+1,:), [size(A) size(map, 2)]); % Assume indexed from 0
elseif size(A, 3) == 4
    if lower(name(end)) == 'f'
        % TIFF in CMYK colourspace - convert to RGB
        if isfloat(A)
            A = A * 255;
        else
            A = single(A);
        end
        A = 255 - A;
        A(:,:,4) = A(:,:,4) / 255;
        A = uint8(A(:,:,1:3) .* A(:,:,[4 4 4]));
    else
        % Assume 4th channel is an alpha matte
        alpha = A(:,:,4);
        A = A(:,:,1:3);
    end
end
if ~isempty(alpha)
    % Apply transprency over a grey checkerboard pattern
    if isa(alpha, 'uint8')
        alpha = double(alpha) / 255;
    end
    A = double(A) .* alpha(:,:,ones(1, size(A, 3)));
    sqSz = max(size(alpha));
    sqSz = floor(max(log(sqSz / 100), 0) * 10 + 1 + min(sqSz, 100) / 20);
    grid = repmat(85, ceil(size(alpha) / sqSz));
    grid(2:2:end,1:2:end) = 171;
    grid(1:2:end,2:2:end) = 171;
    grid = kron(grid, ones(sqSz));
    alpha = grid(1:size(A, 1),1:size(A, 2)) .* (1 - alpha);
    A = uint8(A + alpha(:,:,ones(1, size(A, 3))));
end
end

%% Read (potentially) multi-frame image to uint8 rgb array
function A = imread_rgb_multi(name)
try
    % Get file info
    info = imfinfo(name);
catch
    % Revert to standard case
    A = imread_rgb(name);
    return
end
if numel(info) < 2
    % Single image
    A = imread_rgb(name);
else
    % Multi-frame image
    switch lower(info(1).Format)
        case 'gif'
            [A, map] = imread(name, 'frames', 'all');
            if ~isempty(map)
                map = uint8(map * 256 - 0.5); % Convert to uint8 for storage
                A = reshape(map(uint32(A)+1,:), [size(A) size(map, 2)]); % Assume indexed from 0
                A = permute(A, [1 2 5 4 3]);
            end
        case {'tif', 'tiff'}
            A = cell(numel(info), 1);
            for a = 1:numel(A)
                [A{a}, map] = imread(name, 'Index', a, 'Info', info);
                if ~isempty(map)
                    map = uint8(map * 256 - 0.5); % Convert to uint8 for storage
                    A{a} = reshape(map(uint32(A{a})+1,:), [size(A) size(map, 2)]); % Assume indexed from 0
                end
                if size(A{a}, 3) == 4
                    % TIFF in CMYK colourspace - convert to RGB
                    if isfloat(A{a})
                        A{a} = A{a} * 255;
                    else
                        A{a} = single(A{a});
                    end
                    A{a} = 255 - A{a};
                    A{a}(:,:,4) = A{a}(:,:,4) / 255;
                    A{a} = uint8(A(:,:,1:3) .* A{a}(:,:,[4 4 4]));
                end
            end
        otherwise
            % Multi-frame not supported for this format
            A = imread_rgb(name);
    end
end
end

%% Get the names of all images in a directory
function L = get_im_names()
D = dir;
n = 0;
L = cell(size(D));
% Go through the directory list
for a = 1:numel(D)
    % Check if file is a supported image type
    if numel(D(a).name) > 4 && ~D(a).isdir && (any(strcmpi(D(a).name(end-3:end), {'.png', '.tif', '.jpg', '.bmp', '.ppm', '.pgm', '.pbm', '.gif', '.ras'})) || any(strcmpi(D(a).name(end-4:end), {'.tiff', '.jpeg'})))
        n = n + 1;
        L{n} = D(a).name;
    end
end
L = L(1:n);
end

%% Parse inputs
function [map, layout, gap, indices, lims, figSize] = parse_inputs(inputs)

% Set defaults
map = [];
layout = [];
gap = 0;
indices = -1;
lims = 0;
figSize = get(0, 'ScreenSize');
figSize = figSize(3:4);

% Check for map and display range
for b = 1:numel(inputs)
    if ~isnumeric(inputs{b})
        b = b - 1;
        break;
    end
    if size(inputs{b}, 2) == 3
        map = inputs{b};
    elseif numel(inputs{b}) < 3
        lims = inputs{b};
    end
end

% Go through option pairs
for a = b+1:2:numel(inputs)
    switch lower(inputs{a})
        case 'map'
            map = inputs{a+1};
            if ischar(map)
                map = feval(map, 256);
            end
        case {'size', 'grid'}
            layout = inputs{a+1};
        case {'gap', 'border'}
            gap = inputs{a+1};
        case 'indices'
            indices = inputs{a+1};
        case {'lims', 'displayrange'}
            lims = inputs{a+1};
        case 'figuresize'
            figSize = inputs{a+1};
        case 'figurehandle'
            figSize = get(inputs{a+1}, 'Position');
            figSize = figSize(3:4);            
        otherwise
            error('Input option %s not recognized', inputs{a});
    end
end

if isempty(map)
   map = gray(256);
end
end

%% Rescale RGB images to the correct limits
function A = rescale_rgb(A, lims)
if size(A, 3) == 3 && ~isequal(lims, default_limits(A))
    A = rescale(A, lims);
end
end

%% Return default limits for the image type
function lims = default_limits(A)
if size(A, 3) == 1
    lims = min_max(A);
else
    lims = [0 1];
    if ~isfloat(A)
        lims = lims * double(intmax(class(A)));
    end
end
end

%% Return minimum and maximum values
function lims = min_max(A)
M = isfinite(A);
if all(M(:))
    lims = double([min(A(:)) max(A(:))]);
else
    lims = double([min(A(M)) max(A(M))]);
end
if isempty(lims)
    lims = [0 1];
elseif lims(1) == lims(2)
    lims(2) = lims(1) + 1;
end
end

%% Determine if images are in a cell array or stream object (e.g. imstream)
function tf = is_cell_or_stream(I)
tf = true;
if iscell(I)
    return;
end
try
    % Stream objects must overload nueml() and {} indexing
    assert(any(strcmp('numel', methods(I))));
    assert(isnumeric(I{1}));
    return;
catch
end
tf = false;
end

%% Determine if using hg2
function tf = ishg2(fig)
tf = isa(fig, 'matlab.ui.Figure');
end
