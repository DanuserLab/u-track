%SC  Display/output truecolor images with a range of colormaps
%
% Examples:
%   sc(image)
%   sc(image, limits)
%   sc(image, map)
%   sc(image, limits, map)
%   sc(image, map, limits)
%   sc(..., col1, mask1, col2, mask2,...)
%   out = sc(...)
%   [out, clim map] = sc(...)
%   sc
%
% Generates a truecolor RGB image based on the input values in 'image' and
% any maximum and minimum limits specified, using the colormap specified.
% The image is displayed on screen if there is no output argument.
% 
% SC has these advantages over MATLAB image rendering functions:
%   - images can be displayed or output; makes combining/overlaying images
%     simple.
%   - images are rendered/output in truecolor (RGB [0,1]); no nasty
%     discretization of the input data.
%   - many special, built-in colormaps for viewing various types of data.
%   - linearly interpolates user defined linear and non-linear colormaps.
%   - no border and automatic, integer magnification (unless figure is
%     docked or maximized) for better display.
%   - multiple images can be generated simultaneously.
%
% For a demonstration, simply call SC without any input arguments.
%
% If you don't like the display behaviour of SC, e.g. the fixed aspect
% ratio, no border and axes, then you can replicate the behaviour of
% IMAGESC using the wrapper function IMSC.
%
% IN:
%   image - MxNxCxP or 3xMxNxP image array. MxN are the dimensions of the
%           image(s), C is the number of channels, and P the number of
%           images. If P > 1, images can only be exported, not displayed.
%   limits - [min max] sets the colorbar limits; values in image are
%            clamped to this range.
%   map - Kx3 or Kx4 user defined colormap matrix, where the optional 4th
%         column is the relative distance between colours along the scale,
%         or a string containing the name of the colormap to use to create
%         the output image. Default: 'none', which is RGB for 3-channel
%         images, grayscale otherwise.  Conversion of multi-channel images
%         to intensity for intensity-based colormaps is done using the L2
%         norm. Most MATLAB colormaps are supported. All named colormaps
%         can be reversed by prefixing '-' to the string. This maintains
%         integrity of the colorbar. Special, non-MATLAB colormaps are:
%      'hicontrast' - a high contrast colormap for intensity images that
%                     maintains intensity scale when converted to
%                     grayscale, for example when printing in black &
%                     white.
%      'prob' - first channel is plotted as hue, and the other channels
%               modulate intensity. Useful for laying probabilites over
%               images.
%      'prob_jet' - first channel is plotted as jet colormap, and the other
%                   channels modulate intensity.
%      'gray_jet' - first channel is plotted as jet colormap, second
%                   channel as gray, and third channel as alpha between the
%                   two.
%      'diff' - intensity values are marked blue for > 0 and red for < 0.
%               Darker colour means larger absolute value. For multi-
%               channel images, the L2 norm of the other channels sets
%               green level. 3 channel images are converted to YUV and
%               images with more that 3 channels are projected onto the
%               principle components first.
%      'compress' - compress many channels to RGB while maximizing
%                   variance.
%      'contrast' - apply the CONTRAST function to create an approximately
%                   equal intensity distribution.
%      'complex' - display a complex matrix as hue for argument and
%                  intensity for magnitude (darker = larger).
%      'flow' - display two channels representing a 2d Cartesian vector as
%               hue for angle and intensity for magnitude (darker colour
%               indicates a larger magnitude).
%      'phase' - first channel is intensity, second channel is phase in
%                radians. Darker colour means greater intensity, hue
%                represents phase from 0 to 2 pi.
%      'stereo' - pair of concatenated images used to generate a red/cyan
%                 anaglyph.
%      'stereo_col' - pair of concatenated RGB images used to generate a
%                     colour anaglyph.
%      'segment' - given concatenated image and segmentation (index image),
%                  colour each segment with the mean colour of the image in
%                  that segment.
%      'rand' - gives an index image a random colormap. Useful for viewing
%               segmentations.
%      'rgb2gray' - converts an RGB image to grayscale in the same fashion
%                   as MATLAB's rgb2gray (in the image processing toolbox).
%   col/mask pairs - Pairs of parameters for coloring specific parts of the
%                    image differently. The first (col) parameter can be
%                    a MATLAB color specifier, e.g. 'b' or [0.5 0 1], or
%                    one of the colormaps named above, or an MxNx3 RGB
%                    image. The second (mask) paramater should be an MxN
%                    logical array indicating those pixels (true) whose
%                    color should come from the specified color parameter.
%                    If there is only one col parameter, without a mask
%                    pair, then mask = any(isnan(I, 3)), i.e. the mask is
%                    assumed to indicate the location of NaNs. Note that
%                    col/mask pairs are applied in order, painting over
%                    previous pixel values.
%
% OUT:
%   out - MxNx3xP truecolour (double) RGB image array in range [0, 1].
%   clim - 1x2 vector of colorbar limits used. [] if this mapping is not
%          available, e.g. if output image is multi-spectral.
%   map - 256x3 colormap of the output. [] if clim is []. 
%
% See also IMAGE, IMAGESC, IMSC, COLORMAP, COLORBAR.
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

% Copyright: Oliver Woodford, 2007-2013

function [I, limits, map] = sc(I, varargin)
%% Check for arguments
if nargin == 0
    % If there are no input arguments then run the demo
    if nargout > 0
        error('Output expected from no inputs!');
    end
    demo; % Run the demo
    return
end

%% Size our image(s)
[y, x, c, n] = size(I);
I = reshape(I, y, x, c, n);

%% Check if image is given with RGB colour along the first dimension
if y == 3 && c > 3
    % Flip colour to 3rd dimension
    I = permute(I, [2 3 1 4]);
    [y, x, c, n] = size(I);
end

%% Don't do much if I is empty
if isempty(I)
    if nargout == 0
        % Clear the current axes if we were supposed to display the image
        cla; axis off;
    else
        % Create an empty array with the correct dimensions
        I = zeros(y, x, (c~=0)*3, n);
    end
    return
end

%% Check for multiple images
% If we have a non-singleton 4th dimension we want to display the images in
% a 3x4 grid and use buttons to cycle through them
if n > 1
    % Return transformed images in an YxXx3xN array
    A = zeros(y, x, 3, n);
    for a = 1:n
        A(:,:,:,a) = sc(I(:,:,:,a), varargin{:});
    end
    I = A;
    if nargout == 0
        clear I
        if n > 12
            sz = [3 4];
        else
           sz = []; 
        end
        imdisp(A, 'Size', sz, 'Border', 0.01); 
    end
    return
end

%% Parse the input arguments coming after I (1st input)
[I, map, limits, mask] = parse_inputs(I, varargin, y, x);

%% Call the rendering function
if isnumeric(map) && (size(map, 2) == 3 || size(map, 2) == 4)
    % Table-based colormap
    [I, limits, map] = table_wrapper(I, map, limits);
elseif ischar(map)    
    % Predefined colormap
    [I, limits, map] = colormap_switch(I, map, limits);
else
    error('map input invalid.');
end

%% Update any masked pixels
I = reshape(I, y*x, 3);
for a = 1:size(mask, 2)
    I(mask{2,a},1) = mask{1,a}(:,1);
    I(mask{2,a},2) = mask{1,a}(:,2);
    I(mask{2,a},3) = mask{1,a}(:,3);
end
I = reshape(I, [y x 3]); % Reshape to correct size


%% Only display if the output isn't used
if nargout == 0
    % Display the image
    hIm = imdisp(I, 'Border', 0);
    % Set the colormap (if valid)
    if ~isempty(limits)
        try
            colormap(gca(), map);
        catch
            set(gcf(), 'Colormap', map);
        end
        if limits(1) == limits(2)
            limits(2) = limits(2) + 1;
        end
        set(gca(), 'CLim', limits);
        set(hIm, 'CDataMapping', 'scaled');
    end
    % Don't print out the matrix if we've forgotten the ";"
    clear I
end
return

%% Colormap switch
function [I, limits, map] = colormap_switch(A, cmap, limits)
% Convert complex planes to new colour planes
[y, x, c] = size(A);
A = full(double(A));
if ~isreal(A)
    A = reshape(A, y, x, 1, c);
    A = cat(3, real(A), imag(A));
    c = c * 2;
end
% Reshape
I = reshape(A, y*x, 1, c);
% If map starts with a '-' sign, invert the colormap
reverseMap = cmap(1) == '-';
% If the map ends with a '*', attempt to make map convert linearly to
% grayscale
grayMap = cmap(end) == '*';
% Large switch statement for all the colormaps
switch lower(cmap(1+reverseMap:end-grayMap))
%% Prism
    case 'prism'
        % Similar to the MATLAB internal prism colormap, but only works on
        % index images, assigning each index (or rounded float) to a
        % different colour
        [I, limits] = index_im(I);
        % Generate prism colormap
        map = prism(6);
        if reverseMap
            map = map(end:-1:1,:); % Reverse the map
        end
        % Lookup the colours
        I = mod(I, 6) + 1;
        I = map(I,:);
%% Rand
    case 'rand'
        % Assigns a random colour to each index
        [I, limits, num_vals] = index_im(I);
        % Generate random colormap
        map = rand(num_vals, 3);
        % Lookup the colours
        I = map(I,:);
%% Diff
    case 'diff'
        % Show positive as blue and negative as red, white is 0
        switch c
            case 1
                I(:,2:3) = 0;
            case 2
                % Second channel can only have absolute value
                I(:,3) = abs(I(:,2));
            case 3
                % Diff of RGB images - convert to YUV first
                I = rgb2yuv(I);
                I(:,3) = sqrt(sum(I(:,2:end) .^ 2, 2)) ./ sqrt(2);
            otherwise
                % Use difference along principle component, and other
                % channels to modulate second channel
                I = calc_prin_comps(I);
                I(:,3) = sqrt(sum(I(:,2:end) .^ 2, 2)) ./ sqrt(c - 1);
                I(:,4:end) = [];
        end
        % Generate limits
        if isempty(limits)
            limits = [min(I(:,1)) max(I(:,1))];
        end
        limits = max(abs(limits));
        if limits > 0
            % Scale
            if c > 1
                I(:,[1 3]) = I(:,[1 3]) / limits;
            else
                I = I / (limits * 0.5);
            end
        else
            limits = 1;
        end
        % Colour
        M = I(:,1) > 0;
        I(:,2) = -I(:,1) .* ~M;
        I(:,1) = I(:,1) .* M;
        if reverseMap
            % Swap first two channels
            I = I(:,[2 1 3]);
        end
        %I = 1 - I * [1 0.4 1; 0.4 1 1; 1 1 0.4]; % (Green/Red)
        I = 1 - I * [1 1 0.4; 0.4 1 1; 1 0.4 1]; % (Blue/Red)
        I = min(max(I, 0), 1);
        limits = [-limits limits]; % For colorbar
%% Flow and Complex
    case {'flow', 'complex'}
        % Calculate amplitude and phase, and use 'phase'
        if c ~= 2
            error('''%s'' requires two channels', cmap);
        end
        % Compute magnitude
        A = sqrt(sum(I .^ 2, 3));
        % Colormap dependent stuff
        if lower(cmap(1+reverseMap)) == 'c'
            % Switch axes for 'complex'
            I = cat(3, I(:,:,2), -I(:,:,1));
            % Limits for 'complex'
            if isempty(limits)
                limits = [min(A) max(A)];
            end
        else
            % Limits for 'flow'
            if isempty(limits)
                limits = [min(A) max(A)*2];
            else
                limits = [0 max(abs(limits)*sqrt(2))*2];
            end
        end
        I(:,1) = atan2(I(:,2), I(:,1));
        I(:,2) = A;
        if reverseMap
            % Invert the amplitude
            I(:,2) = -I(:,2);
            limits = -limits([2 1]);
        end
        I = phase_helper(I, limits, 2); % Last parameter tunes how saturated colors can get
        % Set NaNs (unknown flow) to 0
        I(isnan(I)) = reverseMap;
        limits = []; % This colormap doesn't have a valid colorbar
%% Phase
    case 'phase'
        % Plot amplitude as intensity and angle as hue
        if c < 2
            error('''phase'' requires two channels');
        end
        if isempty(limits)
            limits = [min(I(:,1)) max(I(:,1))];
        end
        if reverseMap
            % Invert the phase
            I(:,2) = -I(:,2);
        end
        I = I(:,[2 1]);
        if diff(limits)
            I = phase_helper(I, limits, 1.3); % Last parameter tunes how saturated colors can get
        else
            % No intensity - just cycle hsv
            I = real2rgb(mod(I(:,1) / (2 * pi), 1), 'hsv');
        end
        limits = []; % This colormap doesn't have a valid colorbar
%% RGB2Gray
    case {'rgb2gray', 'rgb2grey'}
        % Compress RGB to greyscale
        [I, limits] = rgb2gray(I, limits, reverseMap);
%% RGB2YUV
    case 'rgb2yuv'
        % Convert RGB to YUV - not for displaying or saving to disk!
        [I, limits] = rgb2yuv(I);
%% YUV2RGB
    case 'yuv2rgb'
        % Convert YUV to RGB - undo conversion of rgb2yuv
        if c ~= 3
            error('''yuv2rgb'' requires a 3 channel image');
        end
        I = I(:,:) * [1 1 1; 0, -0.39465, 2.03211; 1.13983, -0.58060  0];
        I = rescale(I, limits);
        limits = []; % This colormap doesn't have a valid colorbar
%% Prob
    case 'prob'
        % Plot first channel as grey variation of 'bled' and modulate
        % according to other channels
        if c > 1
            A = rgb2gray(I(:,:,2:end), [], false);
            I = I(:,:,1);
        else
            A = 0.5;
        end
        cmap2 = cmap;
        cmap2((1:4)+reverseMap) = 'bled';
        [I, limits] = real2rgb(I, cmap2, limits);
        I = rescale(A + I, [-0.1 1.3]);
%% Prob_jet
    case 'prob_jet'
        % Plot first channel as 'jet' and modulate according to other
        % channels
        if c > 1
            A = rgb2gray(I(:,:,2:end), [], false);
            I = I(:,:,1);
        else
            A = 0.5;
        end
        cmap2 = cmap;
        cmap2(reverseMap+(1:5)) = [];
        M = isnan(I);
        [I, limits] = real2rgb(I, cmap2, limits);
        I(repmat(M, [1 1 3])) = 0;
        I = rescale(A + I, [0.2 1.8]);
%% Gray_jet
    case 'gray_jet'
        % Plot first channel as 'jet', second channel as 'gray', and third
        % channel as alpha between the two.
        if c ~= 3
            error('gray_jet requires a 3 channel image');
        end
        J = real2rgb(I(:,:,1), 'jet', limits);
        G = real2rgb(I(:,:,2), 'gray', limits);
        A = real2rgb(I(:,:,3), 'gray', limits);
        I = J .* A + G .* (1 - A);
        limits = [];
%% Compress
    case 'compress'
        % Compress to RGB, maximizing variance
        % Determine and scale to limits
        I = rescale(I, limits);
        if reverseMap
            % Invert after everything
            I = 1 - I;
        end
        % Zero mean
        meanCol = mean(I, 1);
        isBsx = exist('bsxfun', 'builtin');
        if isBsx
            I = bsxfun(@minus, I, meanCol);
        else
            I = I - meanCol(ones(x*y, 1, 'uint8'),:);
        end
        % Calculate top 3 principle components
        I = calc_prin_comps(I, 3);
        % Normalize each channel independently
        if isBsx
            I = bsxfun(@minus, I, min(I, [], 1));
            I = bsxfun(@times, I, 1./max(I, [], 1));
        else
            for a = 1:3
                I(:,a) = I(:,a) - min(I(:,a));
                I(:,a) = I(:,a) / max(I(:,a));
            end
        end
        % Put components in order of human eyes' response to channels
        I = I(:,[2 1 3]);
        limits = []; % This colormap doesn't have a valid colorbar
%% Contrast
    case 'contrast'
        % Use MATLAB's CONTRAST function to spread the gray levels evenly
        [I, limits] = intensity(I, limits, 0);
        map = contrast(I, 256);
        if reverseMap
            map = map(end:-1:1,:);
        end
        I = real2rgb(I, map, [0 1]);
%% Stereo (anaglyph)
    case 'stereo'
        % Convert 2 colour images to intensity images
        % Show first channel as red and second channel as cyan
        A = rgb2gray(I(:,1:floor(end/2)), limits, false);
        I = rgb2gray(I(:,floor(end/2)+1:end), limits, false);
        if reverseMap
            I(:,2:3) = A(:,1:2); % Make first image cyan
        else
            I(:,1) = A(:,1); % Make first image red
        end
        limits = []; % This colormap doesn't have a valid colorbar
%% Coloured anaglyph
    case 'stereo_col'
        if c ~= 6
            error('''stereo_col'' requires a 6 channel image');
        end
        I = rescale(I, limits);
        % Red channel from one image, green and blue from the other
        if reverseMap
            I(:,1) = I(:,4); % Make second image red
        else
            I(:,2:3) = I(:,5:6); % Make first image red
        end
        I = I(:,1:3);
        limits = []; % This colormap doesn't have a valid colorbar
%% Coloured segments
    case 'segment'
        % Last channel assumed to be an index image. Output each segment as
        % a mean of the other channels in the segment.
        if c < 2
            error('''segment'' requires at least 2 channels');
        end
        [S, num_vals, num_vals] = index_im(I(:,end));
        % Generate the means
        R = 1 ./ (accumarray(S, 1, [num_vals 1]) + 1e-300);
        for a = c-1:-1:1
            J(:,1,a) = accumarray(S, I(:,a), [num_vals 1]) .* R;
        end
        if c == 4
            J = rescale(J, limits);
            limits = [];
        else
            [J, limits, map] = table_wrapper(J, 'gray', limits);
        end
        % Generate the image
        I = J(S,:);
%% None
    case 'none'
        % No colormap - just output the image
        if c == 3
            I = rescale(I, limits);
            limits = [];
            if reverseMap
                I = 1 - I;
            end
        else
            cmap((1:4)+reverseMap) = 'gray';
            [I, limits, map] = table_wrapper(I, cmap, limits);
        end
%% Grey
    case 'grey'
        cmap(3+reverseMap) = 'a';
        [I, limits, map] = table_wrapper(I, cmap, limits);
%% Colourcube
    case {'colorcube2', 'colourcube2'}
        % Psychedelic colormap inspired by MATLAB's version
        [I, limits] = intensity(I, limits, reverseMap); % Intensity map
        step = 4;
        I = I * (step * (1 - eps));
        J = I * step;
        K = floor(J);
        I = cat(3, mod(K, step)/(step-1), J - floor(K), mod(floor(I), step)/(step-1));
%% Functional colormap
    otherwise
        [I, limits, map] = table_wrapper(I, cmap, limits);
end
I = reshape(I, y, x, 3);
% Fix up the colormap
if nargout > 2
    if isempty(limits)
        map = [];
    elseif ~exist('map', 'var')
        map = (0:255) * ((limits(2) - limits(1)) / 255) + limits(1);
        map = squeeze(colormap_switch(map, cmap, limits));
    end
end
return

%% Parse input variables
function [I, map, limits, mask] = parse_inputs(I, inputs, y, x)

% Check the first two arguments for the colormap and limits
ninputs = numel(inputs);
map = 'none';
limits = [];
mask = 1;
for a = 1:min(2, ninputs)
    if ischar(inputs{a}) && numel(inputs{a}) > 1
        % Name of colormap
        map = inputs{a};
    elseif isnumeric(inputs{a})
        [p, q, r] = size(inputs{a});
        if (p * q * r) == 2
            % Limits
            limits = double(inputs{a});
        elseif p > 1 && (q == 3 || q == 4) && r == 1
            % Table-based colormap
            map = inputs{a};
        else
            break;
        end
    else
        break;
    end
    mask = mask + 1;
end
% Check for following inputs
if mask > ninputs
    mask = cell(2, 0);
    return
end
% Following inputs must either be colour/mask pairs, or a colour for NaNs
if ninputs - mask == 0
    mask = cell(2, 1);
    mask{1} = inputs{end};
    mask{2} = ~all(isfinite(I), 3);
elseif mod(ninputs-mask, 2) == 1
    mask = reshape(inputs(mask:end), 2, []);
else
    error('Error parsing inputs');
end
% Go through pairs and generate
M = false(y, x);
for a = 1:size(mask, 2)
    % Generate any masks from functions
    if isa(mask{2,a}, 'function_handle')
        mask{2,a} = mask{2,a}(I);
    end
    if ~islogical(mask{2,a})
        error('Mask is not a logical array');
    end
    if ~isequal(size(mask{2,a}), [y x])
        error('Mask does not match image size');
    end
    if ischar(mask{1,a})
        if numel(mask{1,a}) == 1
            % Generate colours from MATLAB colour strings
            mask{1,a} = str2color(mask{1,a});
        else
            % Assume it's a colormap name
            mask{1,a} = sc(I, mask{1,a});
        end
    end
    mask{1,a} = reshape(mask{1,a}, [], 3);
    if size(mask{1,a}, 1) ~= y*x && size(mask{1,a}, 1) ~= 1
        error('Replacement color/image of unexpected dimensions');
    end
    if size(mask{1,a}, 1) ~= 1
        mask{1,a} = mask{1,a}(mask{2,a},:);
    end
    M = M | mask{2,a};
end
% Set all masked pixels to NaN, so they're ignored later, e.g. during limit
% computation
I(M(:,:,ones(1, size(I, 3)))) = NaN;
return

%% RGB2gray
function [B, limits] = rgb2gray(A, limits, reverseMap)
% Compress RGB to grayscale
c = size(A, 3);
B = A;
if c == 3
    B = reshape(reshape(B, [], 3) * [0.299; 0.587; 0.114], size(B, 1), size(B, 2));
end
[B, limits] = intensity(B, limits, reverseMap);
% Replicate channels
B = B(:,:,[1 1 1]);
return
        
%% RGB2YUV
function [B, limits] = rgb2yuv(A)
% Convert RGB to YUV - not for displaying or saving to disk!
c = size(A, 3);
if c ~= 3
    error('rgb2yuv requires a 3 channel image');
end
B = reshape(A, [], c) * [0.299, -0.14713, 0.615; 0.587, -0.28886, -0.51498; 0.114, 0.436, -0.10001];
limits = []; % This colormap doesn't have a valid colorbar
return

%% Phase helper
function I = phase_helper(I, limits, n)
I = squeeze(I);
I(:,1) = mod(I(:,1)/(2*pi), 1);
I(:,2) = I(:,2) - limits(1);
I(:,2) = I(:,2) * (n / (limits(2) - limits(1)));
I(:,3) = n - I(:,2);
I(:,[2 3]) = min(max(I(:,[2 3]), 0), 1);
I = hsv2rgb(reshape(I, [], 1, 3));
return

%% Real2rgb wrapper
function [B, limits, map] = table_wrapper(A, cmap, limits)
% Convert to intensity
[B, limits] = intensity(A, limits, 0);
% Convert to colour using the specified colormap
[B, map, map] = real2rgb(B, cmap, [0 1]);
return

%% Intensity maps
function [I, limits] = intensity(I, limits, reverseMap)
% Squash to 1d using L2 norm
if size(I, 3) > 1
    I = sqrt(sum(I .^ 2, 3));
end
% Determine and scale to limits
[I, limits] = rescale(I, limits);
if limits(1) == limits(2)
    limits(2) = limits(1) + 1;
end
if reverseMap
    % Invert after everything
    I = 1 - I;
end
return

%% Index images
function [J, limits, num_vals] = index_im(I)
% Returns an index image
if size(I, 3) ~= 1
    error('Index maps only work on single channel images');
end
J = round(I);
rescaled = any(abs(I - J) > 0.01);
if rescaled
    % Appears not to be an index image. Rescale over 256 indices
    m = min(I);
    m = m * (1 - sign(m) * eps);
    I = I - m;
    I = I * (256 / max(I(:)));
    J = ceil(I);
    num_vals = 256;
elseif nargout > 2
    % Output the number of values
    J = J - (min(J) - 1);
    num_vals = max(J);
end
% These colormaps don't have valid colorbars
limits = [];    
return

%% Calculate principle components
function I = calc_prin_comps(I, numComps)
c = size(I, 3);
I = reshape(I, [], c);
if nargin < 2
    numComps = c;
end
% Do SVD
[I, S] = svd(I, 0);
% Calculate projection of data onto components
S = diag(S);
S = S(1:numComps)';
if exist('bsxfun', 'builtin')
    I = bsxfun(@times, I(:,1:numComps), S);
else
    I = I(:,1:numComps) .* S(ones(size(I, 1), 1, 'uint8'),:);
end
return

%% Demo function to show capabilities of sc
function demo
%% Demo gray & lack of border
figure; fig = gcf; Z = peaks(256); sc(Z);
display_text([...
' Lets take a standard, MATLAB, real-valued function:\n\n    peaks(256)\n\n'...
' Calling:\n\n    figure\n    Z = peaks(256);\n    sc(Z)\n\n'...
' gives (see figure). SC automatically scales intensity to fill the\n'...
' truecolor range of [0 1].\n\n'...
' If your figure isn''t docked, then the image will have no border, and\n'...
' will be magnified by an integer factor (in this case, 2) so that the\n'...
' image is a reasonable size.']);

%% Demo colour image display 
figure(fig); clf;
load mandrill; mandrill = ind2rgb(X, map); sc(mandrill);
display_text([...
' That wasn''t so interesting. The default colormap is ''none'', which\n'...
' produces RGB images given a 3-channel input image, otherwise it produces\n'...
' a grayscale image. So calling:\n\n    load mandrill\n'...
'    mandrill = ind2rgb(X, map);\n    sc(mandrill)\n\n gives (see figure).']);

%% Demo discretization
figure(fig); clf;
subplot(121); sc(Z, 'jet'); label(Z, 'sc(Z, ''jet'')');
subplot(122); imagesc(Z); axis image off; colormap(jet(64)); % Fix the fact we change the default depth
label(Z, 'imagesc(Z); axis image off; colormap(''jet'');');
display_text([...
' However, if we want to display intensity images in color we can use any\n'...
' of the MATLAB colormaps implemented (most of them) to give truecolor\n'...
' images. For example, to use ''jet'' simply call:\n\n'...
'    sc(Z, ''jet'')\n\n'...
' The MATLAB alternative, shown on the right, is:\n\n'...
'    imagesc(Z)\n    axis equal off\n    colormap(jet)\n\n'...
' which generates noticeable discretization artifacts.']);

%% Demo intensity colormaps
figure(fig); clf;
subplot(221); sc(Z, 'hsv'); label(Z, 'sc(Z, ''hsv'')');
subplot(222); sc(Z, 'colorcube'); label(Z, 'sc(Z, ''colorcube'')');
subplot(223); sc(Z, 'hicontrast'); label(Z, 'sc(Z, ''hicontrast'')');
subplot(224); sc(Z-round(Z), 'diff'); label(Z, 'sc(Z-round(Z), ''diff'')');
display_text([...
' There are several other intensity colormaps to choose from. Calling:\n\n'...
'    help sc\n\n'...
' will give you a list of them. Here are several others demonstrated.']);

%% Demo saturation limits & colormap reversal
figure(fig); clf;
subplot(121); sc(Z, [0 max(Z(:))], '-hot'); label(Z, 'sc(Z, [0 max(Z(:))], ''-hot'')');
subplot(122); sc(mandrill, [-0.5 0.5]); label(mandrill, 'sc(mandrill, [-0.5 0.5])');
display_text([...
' SC can also rescale intensity, given an upper and lower bound provided\n'...
' by the user, and invert most colormaps simply by prefixing a ''-'' to the\n'...
' colormap name. For example:\n\n'...
'    sc(Z, [0 max(Z(:))], ''-hot'');\n'...
'    sc(mandrill, [-0.5 0.5]);\n\n'...
' Note that the order of the colormap and limit arguments are\n'...
' interchangable.']);

%% Demo prob
load gatlin;
gatlin = X;
figure(fig); clf; im = cat(3, abs(Z)', gatlin(1:256,end-255:end)); sc(im, 'prob');
label(im, 'sc(cat(3, prob, gatlin), ''prob'')');
display_text([...
' SC outputs the recolored data as a truecolor RGB image. This makes it\n'...
' easy to combine colormaps, either arithmetically, or by masking regions.\n'...
' For example, we could combine an image and a probability map\n'...
' arithmetically as follows:\n\n'...
'    load gatlin\n'...
'    gatlin = X(1:256,end-255:end);\n'...
'    prob = abs(Z)'';\n'...
'    im = sc(prob, ''hsv'') .* sc(prob, ''gray'') + sc(gatlin, ''rgb2gray'');\n'...
'    sc(im, [-0.1 1.3]);\n\n'...
' In fact, that particular colormap has already been implemented in SC.\n'...
' Simply call:\n\n'...
'    sc(cat(3, prob, gatlin), ''prob'');']);

%% Demo colorbar
colorbar;
display_text([...
' SC also makes possible the generation of a colorbar in the normal way, \n'...
' with all the colours and data values correct. Simply call:\n\n'...
'    colorbar\n\n'...
' The colorbar doesn''t work with all colormaps, but when it does,\n'...
' inverting the colormap (using ''-map'') maintains the integrity of the\n'...
' colorbar (i.e. it works correctly) - unlike if you invert the input data.']);

%% Demo combine by masking
figure(fig); clf;
sc(Z, [0 max(Z(:))], '-hot', sc(Z-round(Z), 'diff'), Z < 0);
display_text([...
' It''s just as easy to combine generated images by masking too. Here''s an\n'...
' example:\n\n'...
'    im = cat(4, sc(Z, [0 max(Z(:))], ''-hot''), sc(Z-round(Z), ''diff''));\n'...
'    mask = repmat(Z < 0, [1 1 3]);\n'...
'    mask = cat(4, mask, ~mask);\n'...
'    im = sum(im .* mask, 4);\n'...
'    sc(im)\n\n'...
' In fact, SC can also do this for you, by adding image/colormap and mask\n'...
' pairs to the end of the argument list, as follows:\n\n'...
'    sc(Z, [0 max(Z(:))], ''-hot'', sc(Z-round(Z), ''diff''), Z < 0);\n\n'...
' A benefit of the latter approach is that you can still display a\n'...
' colorbar for the first colormap.']);

%% Demo user defined colormap
figure(fig); clf; sc(abs(Z), rand(10, 3)); colorbar;
display_text([...
' SC can also use user defined colormaps to display indexed images.\n'...
' These can be defined as a linear colormap. For example:\n\n'...
'    sc(abs(Z), rand(10, 3))\n    colorbar;\n\n'...
' Note that the colormap is automatically linearly interpolated.']);

%% Demo non-linear user defined colormap
figure(fig); clf; sc(abs(Z), [rand(10, 3) exp((1:10)/2)']); colorbar;
display_text([...
' Non-linear colormaps can also be defined by the user, by including the\n'...
' relative distance between the given colormap points on the colormap\n'...
' scale in the fourth column of the colormap matrix. For example:\n\n'...
'    sc(abs(Z), [rand(10, 3) exp((1:10)/2)''])\n    colorbar;\n\n'...
' Note that the colormap is still linearly interpolated between points.']);

%% Demo compress
load mri;
mri = D;
figure(fig); clf;
sc(squeeze(mri(:,:,:,1:6)), 'compress');
display_text([...
' SC has many colormaps for displaying multi-dimensional data, including\n'...
' ''flow'' for optic flow fields, ''phase'' for edge maps with phase, and\n'...
' ''stereo'' for generating anaglyphs.\n\n'...
' For images with more than 3 channels, SC can compress these images to RGB\n'...
' while maintaining the maximum amount of variance in the data. For\n'...
' example, this 6 channel image:\n\n'...
'    load mri\n    mri = D;\n    sc(squeeze(mri(:,:,:,1:6), ''compress'')']);

%% Demo texture map
figure(fig); clf;
surf(Z, sc(Z, 'hicontrast'), 'edgecolor', 'none');
display_text([...
' Further benefits of SC outputting the image as an array (on top of being\n'...
' able to combine images) are that the image can be saved straight to disk\n'...
' using imwrite(), or can be used to texture map a surface, thus:\n\n'...
'    tex = sc(Z, ''hicontrast'');\n'...
'    surf(Z, tex, ''edgecolor'', ''none'');']);

%% Demo multiple images
close(fig); % Only way to get round loss of focus (bug?)
fig = figure(); clf; sc(mri, 'bone');
display_text([...
' Finally, SC can process multiple images when passed in as a 4d array,\n'...
' either for export or for display. In the latter case, images are\n'...
' displayed as a montage of up to 12 images. For example:\n\n'...
'    sc(mri, ''bone'')']);

display_text([...
' Any additional images can then be viewed by scrolling through them using\n'...
' using the scroll (arrow) keys. Try it and see.']);

clc; fprintf('End of demo.\n');
return

%% Some helper functions for the demo
function display_text(str)
clc;
fprintf([str '\n\n']);
fprintf('Press a key to go on.\n');
figure(gcf);
waitforbuttonpress;
return

function label(im, str)
text(size(im, 2)/2, size(im, 1)+12, str,...
    'Interpreter', 'none', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
return
