function [posAlongDir,deltaPosAlongDir,prefDir] = projectCoordOntoDir(...
    trackCoord,deltaCoord,centerCoord,prefDir)
%PROJECTCOORDALONGDIR reduces coordinates into their 1D projection along a line
%
%SYNPOSIS [posAlongDir,deltaPosAlongDir,prefDir] = projectCoordOntoDir(...
%    trackCoord,deltaCoord,centerCoord,prefDir)
%
%INPUT  trackCoord      : # time points x dimensionality array of
%                         coordinates.
%       deltaCoord      : # time points x dimensionality array of
%                         coordinate standard deviations.
%                         Optional. Default: zeros.
%       centerCoord     : Coordinates of point used as direction
%                         reference. Optional. Default: []. Ignored if
%                         prefDir is input.
%       prefDir         : Direction to project coordinates onto.
%                         Optional. Default: preferred direction of motion
%                         as calculated from scatter of positions.
%
%OUTPUT posAlongDir     : 1D projected coordinates along specified
%                         direction.
%       deltaPosAlongDir: Standard deviation of projected coordinates.
%       prefDir         : Preferred direction onto which coordinates are
%                         projected, = input prefDir if input.
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

%% input

if nargin < 1 || isempty(trackCoord)
    disp('--projectCoordOntoDir: Please input coordinates to be projected');
    return
end

if nargin < 2 || isempty(deltaCoord)
    deltaCoord = zeros(size(trackCoord));
end

if nargin < 3 || isempty(centerCoord)
    centerCoord = [];
end

if nargin < 4 || isempty(prefDir)
    calcDir = 1;
else
    calcDir = 0;
    prefDir = prefDir / norm(prefDir);
end

%% preferred direction

if calcDir

    %perform eigenvalue decomposition on scatter of positions
    [eigenVec,eigenVal] = eig(nancov(trackCoord));
    eigenVal = diag(eigenVal);

    %find track's direction of motion as the eigenvector with the largest
    %eigenvalue
    prefDir = eigenVec(:,eigenVal==max(eigenVal));
%     prefDir = eigenVec(:,eigenVal==min(eigenVal));

    if ~isempty(centerCoord)

        %calculate center of mass of track
        trackCenterCoord = nanmean(trackCoord);

        %calculate vector from nucleus center to track center
        vecNucCenter2Track = trackCenterCoord - centerCoord;

        %modify direction such that +ve movement points away from the center
        prefDir = sign(vecNucCenter2Track * prefDir) * prefDir;

    end

    %normalize the direction vector
    prefDir = prefDir / norm(prefDir);

end

%% projection

%calculate the dot product of positions with the direction vector
%calculate the projections' standard deviations
posAlongDir = trackCoord * prefDir;
deltaPosAlongDir = sqrt( deltaCoord.^2 * prefDir.^2 );

%% ~~~ the end ~~~

