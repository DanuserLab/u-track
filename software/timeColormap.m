function cmap = timeColormap(nTimepoints,colorseq)
%TIMECOLORMAP is Khuloud's time color map
%
% SYNOPSIS: cmap = timeColormap(nTimepoints)
%
% INPUT nTimepoints: (opt) length of colormap. Default: 64
%       colorseq   : (opt) string indicating sequence of colors to use. 
%                    Any permutation of red,green and blue, e.g. rgb, gbr, etc.
%                    Default: gbr.
%
% OUTPUT cmap : nTimepoints-by-3 colormap that changes from green to blue to red
%
% REMARKS This code has been lifted out of Khuloud's plotTracks2D
%         In plotTracks2D, the colormap returns one fewer entry than
%         timepoints, because tracks are drawn between timepoints. This has
%         been changed.
%
% created with MATLAB ver.: 7.6.0.324 (R2008a) on Windows_NT
%
% created by: Jonas Dorn
% DATE: 16-May-2008
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

%% input
if nargin < 1 || isempty(nTimepoints)
    nTimepoints = 64;
end

%allow variation in sequence of colors - KJ
if nargin < 2 || isempty(colorseq)
    colorseq = 'gbr';
end

%determine color sequence based on input
colorCode = zeros(nTimepoints,3,3);
for i=1:3
    colorTmp = colorseq(i);
    switch colorTmp
        case 'r'
            colorCode(:,:,i) = repmat([1 0 0],nTimepoints,1);
        case 'g'
            colorCode(:,:,i) = repmat([0 1 0],nTimepoints,1);
        case 'b'
            colorCode(:,:,i) = repmat([0 0 1],nTimepoints,1);
    end
end

%% make colormap

%old code replaced with new code to allow flexibility in color sequence - KJ
%
% %get the fraction of each color in each time interval to be plotted
% numTimePlotOver2 = ceil((nTimepoints-1)/2); %needed to change blue color over time
% redVariation = (0:nTimepoints-1)'/(nTimepoints-1);
% greenVariation = (nTimepoints-1:-1:0)'/(nTimepoints-1);
% blueVariation = [(0:numTimePlotOver2-1)'/(numTimePlotOver2);...
%     (nTimepoints-numTimePlotOver2-1:-1:0)'/(nTimepoints-numTimePlotOver2)];
% 
% %get the overall color per time interval
% cmap = [redVariation greenVariation blueVariation];

%get the fraction of each color in each time interval to be plotted
colorVariation = zeros(nTimepoints,3,3);
numTimePlotOver2 = ceil((nTimepoints-1)/2); %needed to change the middle color over time
colorVariation(:,:,1) = repmat((nTimepoints-1:-1:0)'/(nTimepoints-1),1,3);
colorVariation(:,:,2) = repmat([(0:numTimePlotOver2-1)'/(numTimePlotOver2);...
    (nTimepoints-numTimePlotOver2-1:-1:0)'/(nTimepoints-numTimePlotOver2)],1,3);
colorVariation(:,:,3) = repmat((0:nTimepoints-1)'/(nTimepoints-1),1,3);

%get the overall color per time interval
cmap = sum(colorCode .* colorVariation, 3);


%% ~~~ the end ~~~