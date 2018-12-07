function gaussList = GaussListND(coordList,sigma,center,intNorm,rotation)
%GAUSSLISTND calculates the value of a N-D Gaussian at specific pixel/voxel coordinates
%
% SYNOPSIS gaussList = GaussListND(coordList,sigma,center,intNorm,rotation)
%
% INPUT    coordList : m-by-n list of coordinates, where m is the number of
%                      coordinates and n the number of dimensions
%          sigma     : 1-by-n (or scalar): sigma of Gaussian
%          center    : (opt) 1-by-n vector of center of Gaussian.
%                      Default: zeros(1,n)
%          intNorm   : (opt) switch for how the Gaussian should be normed
%                      Default: 0
%                      0 - no norming. Max of Gaussian == 1
%                      1 - normed so that integral of infinite Gaussian = 1
%          rotation  : (opt) Equal to the number of degree you want the
%                            coordinate to be rotate for. If rotation is
%                            equal to 1, rotation will be random.
%                            Default: 0; 
%                            Rotation is only supported for 2D and 3D case
%
% OUTPUT   gaussList : m-by-1 list of intensities. Intensity is the
%                      integral of the Gaussian over the pixel/voxel
%
% REMARKS  The code assumes that a pixel has the edge length 1!
%
% c: 2/05 jonas
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

%======================
% TEST INPUT
%======================

% check number of input arguments
% nIn = nargin;
% the following doesn't work with Matlab 6.5.0
% error(nargchk(2,4,nIn,'struct'));
% if nargin < 2 || nargin > 5
if nargin < 3
    error('wrong number of input arguments!')
end

% check dimensionality of coordList.
% if isempty(coordList)
%     error('you have to supply a list of coordinates for GaussList23D')
% else
%     [nCoords,nDims] = size(coordList);
% end
[nCoords,nDims] = size(coordList);
if nCoords == 0 || nDims == 0
    error('you have to supply a list of coordinates for GaussList23D')
end

% sigma
% ls = length(sigma);
% switch ls
%     case nDims
%         % make as long as coords
% %         sigma = repmat(sigma,[nCoords,1]);
%     case 1
% %         sigma = repmat(sigma,[nCoords,nDims]);
%     otherwise
%         error('sigma has to be a scalar or a 1-by-n vector!')
% end

% center
if nargin < 3 || isempty(center)
%     center = zeros(nCoords,nDims);
    center = 0;
% else
%     lc = length(center);
%     switch lc
%         case nDims
% %             center = repmat(center, [nCoords,1]);
%         case 1
% %             center = repmat(center, [nCoords,3]);
%         otherwise
%             error('center has to be a scalar or a 1-by-n vector!')
%     end
end

% intNorm
if nargin < 4 || isempty(intNorm)
    intNorm = 0;
end

%rotation
% coordDim = size(coordList,2);
doRotation = true;
if nargin < 5 || isempty(rotation) || rotation == 0
    rotation = 0;
    doRotation = false;
%     alp = 0;
%     bet = 0;
%     delt = 0;
% elseif rotation == 1 && nDims <= 2
elseif nDims == 2
    rotation = floor(rand(1) * 360);
% elseif rotation == 1 && nDims == 3
elseif nDims == 3
    alp = floor(rand(1) * 180);
    bet = floor(rand(1) * 180);
    delt = floor(rand(1) * 360);
% end
% if rotation && (nDims < 2 || nDims > 3)
else
    error('rotation is only supported for 2-3 dimensions')
end
%======================

%======================
% CALC GAUSSLIST
%======================


% instead of calculating Gauss-values for very complicated geometries, we
% make a coordinate transformation so that we can use sigma=1 in all
% dimensions

if doRotation
    
    %Translate center to origin.
%     coordList = coordList - center;
    coordList = bsxfun(@minus,coordList,center);
    
    if nDims == 2
        
        % 2 Dimension rotation.
        %Rotation.
        %Rotation of the coordinate. x' = xcos@ - ysin@. y' = xsin@ + ycos@.
%         tmpX = coordList(:,1) .* cosd(rotation) - coordList(:,2) .* sind(rotation);
%         tmpY = coordList(:,1) .* sind(rotation) + coordList(:,2) .* cosd(rotation);
        c1 = cosd(rotation);
        s1 = sind(rotation);
        coordList = coordList * [c1 s1; -s1 c1];
        
        %Translation back to original center.
%         coordList(:,1) = tmpX(:,1) + center(:,1);
%         coordList(:,2) = tmpY(:,1) + center(:,2);
%         coordList = bsxfun(@plus,[tmpX tmpY],center);
        coordList = bsxfun(@plus,coordList,center);
        
    elseif nDims == 3
        
        % 3 Dimension rotation.
        %Rotation of the coordinate.
        c1 = cos(alp); c2 = cos(bet); c3 = cos(delt);
        s1 = sin(alp); s2 = sin(bet); s3 = sin(delt);
        
        l1 = (c2 * c3) - (c1*s2*s3); l2 = -(c2 * s3) - (c1 * s2 * c3);
        l3 = s1*s2;
        m1 = (s2*s3 + c1*c2*s3); m2 = -(s2*s3) + (c1*c2*c3);
        m3 = -(s1*c2);
        n1 = s1*s3; n2 = s1*c3; n3 = c1;
        %Calculation of my new coordinate in function of the rotation.
        tmpX = coordList(:,1) .* l1 + coordList(:,2) .* l2 + coordList(:,3) .* l3;
        tmpY = coordList(:,1) .* m1 + coordList(:,2) .* m2 + coordList(:,3) .* m3;
        tmpZ = coordList(:,1) .* n1 + coordList(:,3) .* n2 + coordList(:,3) .* n3;
        
        %Translation back to original center - KJ addition to make
        %consistent with 2D case
        %otherwise the code returns nonsense
%         coordList(:,1) = tmpX(:,1) + center(:,1);
%         coordList(:,2) = tmpY(:,1) + center(:,2);
%         coordList(:,3) = tmpZ(:,1) + center(:,3);
        coordList = bsxfun(@plus,[tmpX tmpY tmpZ],center);
        
    end
 
    
end

% 0.5*erfc(-(x+0.5)/sqrt(2))-0.5*erfc(-(x-0.5)/sqrt(2)) gives the integral on the
% pixel at 1 of a Gaussian with mean 0 and sigma 1

% convert coordList to 0/1
% coordList = (coordList - center)./sigma;
coordList = bsxfun(@minus,coordList,center);
coordList = bsxfun(@rdivide,coordList,sigma);

% double coordList as preparation for erfc
%fixed bug: must divide the 0.5 by sigma - KJ
% coordList = cat(3,coordList-0.5./sigma, coordList+0.5./sigma);
halfInvSigma = 0.5./sigma;
halfInvSigma = cat(3,-halfInvSigma,halfInvSigma);
coordList = bsxfun(@plus,coordList,halfInvSigma);

% calculate gaussList
%Jonas was missing the minus sign in erfc. I corrected that - KJ
gaussList = diff(0.5 * erfc(-coordList/sqrt(2)),1,3);
gaussList = prod(gaussList,2);

% norm gaussList
switch intNorm
    case 0
        % "un-norm" Gaussian
        gaussList = gaussList*((2*pi)^(0.5*nDims)*prod(sigma(1,:)));
%     case 1
        % gaussList is already normed
end
