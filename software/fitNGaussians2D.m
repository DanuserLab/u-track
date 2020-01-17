function [F,J] = fitNGaussians2D(x0,image,index,psfSigma)
%FITNGAUSSIANS2D yields F, the difference between an image and a theoretical image produced by N Gaussians, and J, the Jacobian of F.
%
%SYNOPSIS [F,J] = fitNGaussians2D(x0,image,index,psfSigma)
%
%INPUT  x0      : initial guess of PSF positions and amplitudes and
%                 background noise.
%       image   : Image part being analyzed.
%       index   : x,y-indices of pixels considered.
%       psfSigma: Standard deviation of point spread function (in pixels).
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

%OUTPUT F       : Residuals from fitting an image with supplied
%                 Gaussians.
%       J       : The Jacobian matrix of F.
%       errFlag : 0 if function executes normally, 1 otherwise.
%
%REMARKS F = model image - real image, important to know if the sign of the
%residuals matters.
%
%Khuloud Jaqaman, August 2005

%% Output

F = [];
J = [];

%% Input

%check whether correct number of input arguments was used
if nargin ~= 4
    disp('--fitNGaussians2D: Incorrect number of input arguments!');
    return
end

%% Calculating F & J

%extract background intensity from x0 and remove from vector
bgAmp = x0(end);
x0 = x0(1:end-1);

%get number of PSFs considered
numPSF = length(x0)/3;

%reshape 3nx1 vector x0 into nx3 matrix
x0 = reshape(x0,3,numPSF);
x0 = x0';

%extract PSF center positions and amplitudes
psfPos = x0(:,1:2);
psfAmp = x0(:,3);

%find minimum and maximum pixel indices
minIndxX = min(index(:,1));
maxIndxX = max(index(:,1));
minIndxY = min(index(:,2));
maxIndxY = max(index(:,2));

%determine the contribution of each PSF (assuming amplitude 1) to a 
%pixel based on its x-coordinate (needed to calculate F & J)
psfIntegX = zeros(maxIndxX-minIndxX+1,numPSF);
for i=1:numPSF
    psfIntegX(:,i) = GaussListND((minIndxX:maxIndxX)',...
        psfSigma,psfPos(i,1));
end

%determine the contribution of each PSF (assuming amplitude 1) to a 
%pixel based on its y-coordinate (needed to calculate F & J)
psfIntegY = zeros(maxIndxY-minIndxY+1,numPSF);
for i=1:numPSF
    psfIntegY(:,i) = GaussListND((minIndxY:maxIndxY)',...
        psfSigma,psfPos(i,2));
end

%calculate the value of each PSF (assuming amplitude 1) at the 
%x-coordinates of the corners of all pixels (needed to calculate J)
psfValueX = zeros(maxIndxX-minIndxX+2,numPSF);
for i=1:numPSF
    psfValueX(:,i) = exp(-((minIndxX-0.5:maxIndxX+0.5)'...
        -psfPos(i,1)).^2/2/psfSigma^2);
end

%calculate the value of each PSF (assuming amplitude 1) at the
%y-coordinates of the corners of all pixels (needed to calculate J)
psfValueY = zeros(maxIndxY-minIndxY+2,numPSF);
for i=1:numPSF
    psfValueY(:,i) = exp(-((minIndxY-0.5:maxIndxY+0.5)'...
        -psfPos(i,2)).^2/2/psfSigma^2);
end

%get number of pixels in image
numPixel = length(image);

%get xy-indices relative to minimum
relIndxX = index(:,1) - minIndxX + 1;
relIndxY = index(:,2) - minIndxY + 1;

%calculate the value of F at all pixels
F = (sum(repmat(psfAmp,1,numPixel).*psfIntegX(relIndxX,:)'.*psfIntegY(relIndxY,:)',1))' ...
    + repmat(bgAmp,numPixel,1) - image;

%remove pixels with NaN (which means they are out of the cropped image
%area)
indxPixel = find(~isnan(image));

F = F(indxPixel);

if nargout > 1
    %calculate the derivative at all pixels
    J = ones(numPixel,3*numPSF+1); %(last column for background amplitude)
    J(:,1:3:3*numPSF) = repmat(psfAmp',numPixel,1).*(psfValueX(relIndxX,:)-...
        psfValueX(relIndxX+1,:)).*psfIntegY(relIndxY,:); %w.r.t. x
    J(:,2:3:3*numPSF) = repmat(psfAmp',numPixel,1).*(psfValueY(relIndxY,:)-...
        psfValueY(relIndxY+1,:)).*psfIntegX(relIndxX,:); %w.r.t. y
    J(:,3:3:3*numPSF) = psfIntegX(relIndxX,:).*psfIntegY(relIndxY,:); %w.r.t. amp

    %remove pixels with NaN (which means they are out of the cropped image
    %area)
    J = J(indxPixel,:);
end

%% ~~ the end ~~ 

%% OLD CODE

% % J = ones(numPixel,3*numPSF+1);
% % F = ones(numPixel,1);
% % 
% % for i=1:numPixel %for each pixel
% % 
% %     %get xy-indices relative to minimum
% %     relIndxX = index(i,1) - minIndxX + 1;
% %     relIndxY = index(i,2) - minIndxY + 1;
% %     
% %     %calculate the value of F
% %     F(i) = sum(psfAmp.*psfIntegX(relIndxX,:)'.*psfIntegY(relIndxY,:)') ...
% %         + bgAmp - image(i);
% % 
% %     %calculate the derivative wrt x-coordinate
% %     J(i,1:3:3*numPSF) = psfAmp'.*(psfValueX(relIndxX,:)-...
% %         psfValueX(relIndxX+1,:)).*psfIntegY(relIndxY,:)/psfSigma^2;
% % 
% %     %calculate the derivative wrt y-coordinate
% %     J(i,2:3:3*numPSF) = psfAmp'.*(psfValueY(relIndxY,:)-...
% %         psfValueY(relIndxY+1,:)).*psfIntegX(relIndxX,:)/psfSigma^2;
% % 
% %     %calculate the derivative wrt amplitude
% %     J(i,3:3:3*numPSF) = psfIntegX(relIndxX,:).*psfIntegY(relIndxY,:);
% % 
% %     %since the derivative wrt background intensity = 1, this is already
% %     %accounted for in the initial assignment of J.
% %     
% % end



