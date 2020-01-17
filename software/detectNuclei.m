function [movieInfo, dPix]=detectNuclei(I, radius, varargin)
% detectNuclei detect nuclei from a fluorescent image
%
% SYNOPSIS detectNuclei(I, radius)
%
% INPUT
%   I - input image
%
%   radius - radius of the nuclei (in pixels)
%
%   confluent - optional. Boolean specifying wether the nuclei to be
%   detected are organized in a confluent epithelial sheet. Default: false.
%
%   Possible Parameter Structure Field Names:
%       ('FieldName' -> possible values)
%
%       edgeFilter - a string giving the type of edge filter to be applied
%       to the image. Can be chosen between sobel', 'canny', 'prewitt'.
%       Default: sobel.
%
%       sigma - for the canny detection, a scalar giving the standard
%       deviation of the Gaussian kernel on which the filters are based
%
%       p - a scalar
%
%
%       useDblLog - optional.
%
%       doPlot -  a boolean value to display intermediate graphs. Default:
%       false.
%
%
% For more information, see:
% M. Rosa Ng, A. Besser, G. Danuser, and J. S. Brugge, J Cell Biol 2012 199:545-563
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

% Achim Besser, Dec 2010 (last modified Aug 2011)
% Adapted by Sebastien Besson, Nov 2012

% Input check
ip =inputParser;
ip.addRequired('I', @isnumeric)
ip.addRequired('radius', @isposint)
ip.addOptional('confluent', false, @isscalar)
edgeFilters = {'sobel', 'canny', 'prewitt','none'};
ip.addParamValue('edgeFilter', 'sobel', @(x) ismember(x, edgeFilters));
ip.addParamValue('sigma', 2, @isscalar);
ip.addParamValue('p', .01, @isscalar);
ip.addParamValue('useDblLog', true, @isscalar);
ip.addParamValue('doPlot', false, @isscalar);
ip.parse(I, radius, varargin{:});

% Round radius value to next odd value
radius = 2 * floor(ip.Results.radius/2) +1;
doPlot = ip.Results.doPlot;

[rowsOrg, colsOrg]=size(I);

%**********************************************************************
% 1.Step: crop non-zero part of the image                             *
%**********************************************************************

% Crop an inner part of the image, the boundary might contain zero
% pixels. The following lines find the largest rectangle with non-zero
% values that fits into the image (with equal spacing to the bundary).
realIm=(I~=0);
bwBD=bwboundaries(realIm,8,'noholes');
% one should first finde the maximum, anyways:
bwBD=bwBD{1};

if doPlot==1
    figure, imshow(I,[]), title('Gradient magnitude (gradmag)')
    colormap gray;
    hold on
    plot(bwBD(:,2),bwBD(:,1),'*b')
    hold off;
end

distVals=[bwBD(:,1),bwBD(:,2),rowsOrg-bwBD(:,1),colsOrg-bwBD(:,2)];
dPix=max(min(distVals,[],2));

Icrop =I((1+dPix):(rowsOrg-dPix),(1+dPix):(colsOrg-dPix));
I=Icrop;

%**********************************************************************
% 2. Step: Preprocessing: gradient subtraction, doublelog, G-filter   *
%**********************************************************************

% calculate the gradient information in the image:
if strcmp(ip.Results.edgeFilter,'sobel') || strcmp(ip.Results.edgeFilter,'prewitt')
    hy = fspecial(ip.Results.edgeFilter); %'prewitt' % perform very much the same!
    hx = hy';
    Iy = imfilter(double(I), hy, 'replicate');
    Ix = imfilter(double(I), hx, 'replicate');
    gradmag = sqrt(Ix.^2 + Iy.^2);
elseif strcmp(ip.Results.edgeFilter,'canny')
    [gradmag] = steerableDetector(I,1,ip.Results.sigma);
elseif strcmp(edgeFilter,'none')
    gradmag=zeros(size(I));
end

% Normalize the gradient according to the magnitude of the image intensity:
gradmag = gradmag*max(I(:))/max(vertcat(eps,gradmag(:)));
if doPlot==1
    figure, imshow(gradmag,[]), title('Gradient magnitude (gradmag)')
    colormap gray;
end

minI=min(I(:));
maxI=max(I(:));

numPix=numel(I);
binEdges=linspace(floor(minI),ceil(maxI),1000);
n = histc(I(:),binEdges);

ncum=cumsum(n);
binID=find(ncum>ip.Results.p*numPix,1,'first');
bgVal=binEdges(binID);

% substract the gradient from the original image to enhance the edges:
ImGrad=I-gradmag;
ImGrad(ImGrad<bgVal)=bgVal;  % bring all neg. or very small values back to bg-level!

if doPlot==1
    figure, imagesc(ImGrad), title('Gradient substracted image')
end

% Equalize neighboring maxima by applying a double-log:
bgValScaled=log(log(bgVal));
if ip.Results.useDblLog && isreal(bgValScaled) && bgValScaled>-Inf
    ImGrad=log(log(ImGrad));
else
    display('Couldnt use log(log(.)) to suppress large next to small maxima!');
end

% filter with a gaussian:
Iflt = filterGauss2D(ImGrad, radius);
if doPlot==1
    figure, imagesc(Iflt), title('Gradient substracted and Gauss-filtered image')
    colormap gray;
end

%**********************************************************************
% 4. Step: Find the local maxima in the preprocessed image            *
%**********************************************************************

% Find the local maxima:
se   = strel('disk', radius);
Imax=locmax2d(Iflt,getnhood(se),1);
if doPlot==1
    figure, imagesc(Imax), title('Maxima in the image')
    colormap gray;
    figure;
    hist(Imax(Imax(:)>0),1000)
end

% cut off the maxima in the noise:
try
    if ip.Results.confluent
        % This only works for dense sheets!
        level1 = thresholdFluorescenceImageFewBg(Imax(Imax(:)>0), doPlot);
    else
        % This only works for significant amount of Bg, but is then more reliable than the above method!
        level1 = thresholdFluorescenceImage(Imax(Imax(:)>0), doPlot, 1);
    end
catch
    % This shouldn't be used anymore. Instead, one should use the
    % algorithm above!
    display('!!!switched to cutFirstHistMode!!!')
    [~, level1]=cutFirstHistMode(Imax(Imax(:)>0),0);
end
Imax(Imax(:)<level1)=0;
Imax(Imax(:)>0)=1;

% Show the first set of maxima:
se   = strel('disk', round(radius/4));
ImaxDil = imdilate(Imax,se);
if doPlot==1
    Idspl=Icrop;
    Idspl(ImaxDil == 1) = 0;
    figure, imagesc(Idspl), title('Extended detected nuclei')
end

%**********************************************************************
% 4. Step: Find a second set of max using a simple Gauss-filter       *
%**********************************************************************

% Simply filter the original image with a Gaussian, some of these might
% have been cut-off if they are small in size and have a huge gradient!
IfltSimple = filterGauss2D(Icrop, radius);
ImaxSimple=locmax2d(IfltSimple,getnhood(se),1);
[~, level2]=cutFirstHistMode(ImaxSimple(ImaxSimple(:)>0),0);
ImaxSimple(ImaxSimple(:)<level2)=0;
ImaxSimple(ImaxSimple(:)>0)=1;

se   = strel('disk', radius);
ImaxComb = imdilate(Imax+ImaxSimple,se);
labelsMaxROI = bwlabel(ImaxComb>0, 4);
% label each Imax. Find doubled values. These maxima have been falsely
% merged by regiongrowing:
maxID=find(Imax(:)==1);
maxLabels=labelsMaxROI(maxID);

n = histc(maxLabels,1:max(maxLabels));
dblLables=find(n>1);
falseROI=ismember(labelsMaxROI,dblLables);
if doPlot==1
    falseROI(ImaxDil==1)=0;
    figure, imagesc(falseROI), title('Falsely fused regions!')
end
% Clean Labels:
labelsMaxROI(falseROI)=0;

% Set the max that we want to keep into the Labels matrix:
dblLabelInd=ismember(maxLabels,dblLables);
maxIdKeep=maxID(dblLabelInd);

maxL=max(labelsMaxROI(:));
% first fill up the Labels and then append at the end:
newLabelList=[dblLables',(maxL+1):(maxL+length(maxIdKeep)-length(dblLables))];
for k=1:length(maxIdKeep)
    labelsMaxROI(maxIdKeep(k))=newLabelList(k);
end

%figure, imagesc(labelsMaxROI), title('Labels')

regions = regionprops(labelsMaxROI, 'centroid');

nucPos=round(vertcat((regions(:).Centroid)));
ImaxCombClean=zeros(size(Icrop));
ImaxCombClean(sub2ind(size(Icrop),nucPos(:,2),nucPos(:,1)))=1;

if doPlot==1
    Idspl=Icrop;
    se   = strel('disk', round(radius/4));
    ImaxCombCleanDil = imdilate(ImaxCombClean,se);
    Idspl(ImaxCombCleanDil > 0) = 0;
    figure, imagesc(Idspl), title('Extended nuclei')
end


%**********************************************************************
% 6. Step: Prepare the output                                         *
%**********************************************************************


linId=sub2ind(size(I),nucPos(:,2),nucPos(:,1));
pos=nucPos+dPix;
movieInfo.xCoord(:,1)=pos(:,1);
movieInfo.xCoord(:,2)=radius/2;
movieInfo.yCoord(:,1)=pos(:,2);
movieInfo.yCoord(:,2)=radius/2;
movieInfo.amp(:,1)=I(linId);
movieInfo.amp(:,2)=0;