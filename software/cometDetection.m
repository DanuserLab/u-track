% featuresInfo = cometDetection(img, mask, psfSigma, mode)
%
% Inputs :      img : input image
%              mask : cell mask
%             psfSigma : standard deviation of the Gaussian PSF
%            {mode} : parameters to estimate, default 'xyArtc'
%           {alpha} : alpha-values for statistical test
%           {kSigma} : alpha-values for statistical test
%           {minDist} : minimum distance betwen detected features
%           {filterSigma} : sigma for the steerable filter
%
% Outputs:  featuresInfo : output structure with anisotropic Gaussian
%                          parameters, standard deviations (compatible with
%                          Khuloud's tracker.
%
% Sylvain Berlemont, April 2011
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

function featuresInfo = cometDetection(img, mask, psfSigma, varargin)

% Parse inputs
ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('img', @isnumeric);
ip.addRequired('mask', @islogical);
ip.addRequired('psfSigma', @isscalar);
ip.addParamValue('mode', 'xyArtc', @ischar);
ip.addParamValue('alpha', 0.05, @isscalar);
ip.addParamValue('kSigma', 4, @isscalar);
ip.addParamValue('minDist', .25, @isscalar);
ip.addParamValue('filterSigma',psfSigma*sqrt(2), @isscalar);

ip.parse(img, mask, psfSigma, varargin{:});
mode = ip.Results.mode;
alpha = ip.Results.alpha;
kSigma = ip.Results.kSigma;
minDist = ip.Results.minDist;
filterSigma =ip.Results.filterSigma;

img = double(img);
[nrows ncols] = size(img);

% Filter image with laplacian
bandPassIso = filterLoG(img,psfSigma);
bandPassIso(bandPassIso < 0) = 0;
bandPassIso(~mask) = 0;

% Filter image with steerable filter
[R,T] = steerableDetector(img,2,filterSigma);
 
% Compute the local maxima of the bandpass filtered images
locMaxIso = locmax2d(R, [5 5]);
 
bw = blobSegmentThreshold(bandPassIso,0,0,mask);
labels=bwlabel(bw);

locMaxIso(~bw) = 0;

indMax = find(locMaxIso);
[y x] = ind2sub(size(img), indMax);

P = zeros(size(y, 1), 7);
P(:,1) = x;
P(:,2) = y;
P(:,3) = img(indMax);
P(:,4) = 2*psfSigma;       % sigmaX
P(:,5) = psfSigma;         % sigmaY
P(:,6) = T(indMax)+pi/2;

% % Subresolution detection
% hside = ceil(kSigma * psfSigma);
% npx = (2 * hside + 1)^2;
% xmin = x - hside;
% xmax = x + hside;
% ymin = y - hside;
% ymax = y + hside;

PP =num2cell(P,1);

[xRange,yRange,nzIdx] = arrayfun(@(x0,y0,sigmaX,sigmaY,theta)...
    anisoGaussian2DSupport(x0,y0,sigmaX,sigmaY,theta,kSigma,[ncols nrows]),...
    PP{[1 2 4 5 6]},'UniformOutput',false);
% hside = ceil(kSigma * psfSigma);
% npx = (2 * hside + 1)^2;
xmin = cellfun(@min,xRange);
xmax = cellfun(@max,xRange);
ymin = cellfun(@min,yRange);
ymax = cellfun(@max,yRange);
npx = cellfun(@numel,nzIdx);

isValid = find(xmin >= 1 & xmax <= ncols & ymin >= 1 & ymax <= nrows);

xmin = xmin(isValid);
xmax = xmax(isValid);
ymin = ymin(isValid);
ymax = ymax(isValid);
P = P(isValid,:);

stdP = zeros(size(P));
stdR = zeros(size(P,1),1);

kLevel = norminv(1 - alpha / 2.0, 0, 1); % ~2 std above background

success = false(numel(xmin),1);

for iFeature = 1:numel(xmin)
    mask =labels(ymin(iFeature):ymax(iFeature), xmin(iFeature):xmax(iFeature));
    mask(mask==labels(yRange{iFeature}(fix(end/2)),xRange{iFeature}(fix(end/2)))) = 0;
    
    anisoMask = false(length(yRange{iFeature}),length(xRange{iFeature}));
    anisoMask(nzIdx{iFeature})=true;
    anisoMask(mask~=0)=false;
    npx(iFeature)=nnz(anisoMask);
    crop = img(ymin(iFeature):ymax(iFeature), xmin(iFeature):xmax(iFeature));
    crop(~anisoMask)=NaN;

    
    P(iFeature,7) = min(crop(:)); % background
    P(iFeature,3) = P(iFeature,3) - P(iFeature,7); % amplitude above background
        
    [params, stdParams, ~, res] = fitAnisoGaussian2D(crop, ...
        [0, 0, P(iFeature,3), 3 * P(iFeature,4), P(iFeature,5), ...
        P(iFeature,6), P(iFeature,7)], mode);
        
    % TEST: position must remain in a confined area
    px = floor(floor(size(crop)/2)+1+params(1:2));
    isValid = all(px>=1) & all(px<=size(crop));
    if isValid
        isValid = anisoMask(px(1),px(2));
    end
    
    % TEST: sigmaX > 1
    isValid = isValid & params(4) > 1;
    
%     TEST: goodness-of-fit
%     stdRes = std(res.data(~isnan(res.data)));
%     [~, pval] = kstest(res.data(~isnan(res.data)) ./ stdRes, [], alpha);
    isValid = isValid & res.pval > alpha;

    % TEST: amplitude
    SE_psfSigma_r = (res.std / sqrt(2*(npx(iFeature)-1))) * kLevel;
    psfSigma_A = stdParams(3);
    A_est = params(3);
    df2 = (npx(iFeature) - 1) * (psfSigma_A.^2 + SE_psfSigma_r.^2).^2 ./ (psfSigma_A.^4 + SE_psfSigma_r.^4);
    scomb = sqrt((psfSigma_A.^2 + SE_psfSigma_r.^2) / npx(iFeature));
    T = (A_est - res.std * kLevel) ./ scomb;    
    isValid = isValid & (1 - tcdf(T, df2)) < alpha;
  
    % TEST: extreme value of 
    isValid = isValid & params(4) < 10 * psfSigma;
    
    success(iFeature) = isValid;
    
    P(iFeature,1) = P(iFeature,1) + params(1);
    P(iFeature,2) = P(iFeature,2) + params(2);
    P(iFeature,3) = params(3);
    P(iFeature,4) = params(4);
    P(iFeature,5) = params(5);
    P(iFeature,6) = params(6);
    P(iFeature,7) = params(7);
    
    stdP(iFeature,1) = stdParams(1);
    stdP(iFeature,2) = stdParams(2);
    stdP(iFeature,3) = stdParams(3);
    stdP(iFeature,4) = stdParams(4);
    stdP(iFeature,6) = stdParams(5);
    stdP(iFeature,7) = stdParams(6);
    
    stdR(iFeature) = res.std;
end

P = P(success,:);
stdP = stdP(success,:);

% Remove any detection which has been localised at the same position
isValid = true(size(P,1),1);
idxKD = KDTreeBallQuery(P(:,1:2), P(:,1:2), repmat(minDist, size(P,1), 1));
idxKD = idxKD(cellfun(@(x) length(x)>1, idxKD));
    
for k = 1:length(idxKD);
    stdRes = stdR(idxKD{k});
    isValid(idxKD{k}(stdRes ~= min(stdRes))) = false;
end

P = P(isValid,:);
stdP = stdP(isValid,:);

featuresInfo.xCoord = [P(:,1), stdP(:,1)];
featuresInfo.yCoord = [P(:,2), stdP(:,2)];
featuresInfo.amp = [P(:,3), stdP(:,3)];
featuresInfo.sigmaX = [P(:,4), stdP(:,4)];
featuresInfo.sigmaY = [P(:,5), stdP(:,5)];
featuresInfo.theta = [P(:,6), stdP(:,6)];
featuresInfo.bkg = [P(:,7), stdP(:,7)];
