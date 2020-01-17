%[pstruct, mask, imgLM, imgLoG] = pointSourceDetection3D(vol, sigma, varargin)
%
% Inputs :   
%                 vol : input volume
%               sigma : standard deviation of the Gaussian PSF
%                       If the PSF is anisotropic, 'sigma' should be a two-element vector:
%                       [sigma_xy sigma_z]
%
% Options (as 'specifier'-value pairs): 
%
%              'mode' : parameters to estimate. Default: 'xyzAc'.
%             'alpha' : alpha value used in the statistical tests. Default: 0.05.
%  'alphaLocalMaxima' : alpha value used in selecting candidate local maxima for fitting. Default: 0.05.
%              'mask' : mask of pixels (i.e., cell mask) to include in the detection. Default: all.
%       'FitMixtures' : true|{false}. Toggles mixture-model fitting.
%       'MaxMixtures' : maximum number of mixtures to fit. Default: 5.
%   'RemoveRedundant' : {true}|false. Discard localizations that coincide within 'RedundancyRadius'.
%  'RedundancyRadius' : Radius for filtering out redundant localizatios. Default: 0.25
%     'RefineMaskLoG' : true|{false}. Apply threshold to LoG-filtered img to refine mask of significant pixels.
%   'RefineMaskValid' : true|{false}. Return only mask regions where a significant signal was localized.
%        'ConfRadius' : Confidence radius for positions, beyond which the fit is rejected. Default: 2*sigma
%        'WindowSize' : Window size for the fit. Default: 2*sigma, i.e., [-2*sigma ... 2*sigma]^2
%'LocalMaxWindowSize' : Window size for locmax3d. Default: max(3,roundOddOrEven(ceil(2*sigma([1 1 2])),'odd'))
%
% Outputs:  
%             pstruct : output structure with Gaussian parameters, standard deviations, p-values
%                mask : mask of significant (in amplitude) pixels
%               imgLM : image of local maxima
%              imgLoG : Laplacian of Gaussian-filtered image
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

% Francois Aguet, August 2013 (last modified: 11/05/2013)

function [pstruct, mask, imgLM, imgLoG] = pointSourceDetection3D(vol, sigma, varargin)

% Parse inputs
ip = inputParser;
ip.CaseSensitive = false;
ip.KeepUnmatched = true;
ip.addRequired('vol', @isnumeric);
ip.addRequired('sigma', @isnumeric);
ip.addParamValue('Mode', 'xyzAc', @ischar);
ip.addParamValue('AlphaLocalMaxima', [], @isscalar);%Alpha value used in selection of candidate local maxima
ip.addParamValue('Alpha', 0.05, @isscalar);
ip.addParamValue('Mask', [], @(x) isnumeric(x) || islogical(x));
ip.addParamValue('FitMixtures', false, @islogical);
ip.addParamValue('MaxMixtures', 5, @(x) numel(x)==1 && x>0 && round(x)==x);
ip.addParamValue('RemoveRedundant', true, @islogical);
ip.addParamValue('RedundancyRadius', 0.25, @isscalar);
ip.addParamValue('RefineMaskLoG', false, @islogical);
ip.addParamValue('RefineMaskValid', false, @islogical);
ip.addParamValue('ConfRadius', []); % Default: 2*sigma, see fitGaussians3D.
ip.addParamValue('WindowSize', []); % Default: 2*sigma, see fitGaussians3D.
ip.addParamValue('LocalMaxWindowSize',[]); % Default: max(3,roundOddOrEven(ceil(2*sigma([1 1 2])),'odd'))
ip.addParamValue('ClearMaskBorder',true,@islogical); % 

ip.parse(vol, sigma, varargin{:});

if isempty(ip.Results.AlphaLocalMaxima)
    %Default is to use same as in fit
    alphaLocalMaxima = ip.Results.Alpha;    
else
    alphaLocalMaxima = ip.Results.AlphaLocalMaxima;
end

if ~isa(vol, 'double')
    vol = double(vol);
end

if numel(sigma)==1
    sigma = [sigma sigma];
end

ws = ip.Results.WindowSize;
if isempty(ws)
    ws = ceil(2*sigma);
elseif numel(ws)==1
    ws = [ws ws];
end

localMaxWindowSize = ip.Results.LocalMaxWindowSize;
if(isempty(localMaxWindowSize))
    localMaxWindowSize=max(3,roundOddOrEven(ceil(2*sigma([1 1 2])),'odd'));
end


%-------------------------------------------------------------------------------------------
% Convolutions
%-------------------------------------------------------------------------------------------
% right-hand side of symmetric kernels
gx = exp(-(0:ws(1)).^2/(2*sigma(1)^2));
gz = exp(-(0:ws(2)).^2/(2*sigma(2)^2));
fg = conv3fast(vol, gx, gx, gz);
fu =  conv3fast(vol,    ones(1,ws(1)+1), ones(1,ws(1)+1), ones(1,ws(2)+1));
fu2 = conv3fast(vol.^2, ones(1,ws(1)+1), ones(1,ws(1)+1), ones(1,ws(2)+1));

% Laplacian of Gaussian-filtered input
gx2 = (0:ws(1)).^2 .*gx;
gz2 = (0:ws(2)).^2 .*gz;
fgx2 = conv3fast(vol, gx2, gx, gz);
fgy2 = conv3fast(vol, gx, gx2, gz);
fgz2 = conv3fast(vol, gx, gx, gz2);
imgLoG = (2/sigma(1)^2+1/sigma(2)^2)*fg - ((fgx2+fgy2)/sigma(1)^4 + fgz2/sigma(2)^4);
clear fgx2 fgy2 fgz2;

% Gaussian kernel (spatial)
[x,y,z] = meshgrid(-ws(1):ws(1),-ws(1):ws(1),-ws(2):ws(2));
g = exp(-(x.^2+y.^2)/(2*sigma(1)^2)) .* exp(-z.^2/(2*sigma(2)^2));
n = numel(g);
gsum = sum(g(:));
g2sum = sum(g(:).^2);

% solution to linear system
A_est = (fg - gsum*fu/n) / (g2sum - gsum^2/n);
c_est = (fu - A_est*gsum)/n;

J = [g(:) ones(n,1)]; % g_dA g_dc
C = inv(J'*J);

f_c = fu2 - 2*c_est.*fu + n*c_est.^2; % f-c
RSS = A_est.^2*g2sum - 2*A_est.*(fg - c_est*gsum) + f_c;
clear fg fu2;
RSS(RSS<0) = 0; % negative numbers may result from machine epsilon/roundoff precision
sigma_e2 = RSS/(n-3);

sigma_A = sqrt(sigma_e2*C(1,1));

% standard deviation of residuals
sigma_res = sqrt(RSS/(n-1));
clear fu;

kLevel = norminv(1-alphaLocalMaxima/2.0, 0, 1);

SE_sigma_c = sigma_res/sqrt(2*(n-1)) * kLevel;
df2 = (n-1) * (sigma_A.^2 + SE_sigma_c.^2).^2 ./ (sigma_A.^4 + SE_sigma_c.^4);
scomb = sqrt((sigma_A.^2 + SE_sigma_c.^2)/n);
T = (A_est - sigma_res*kLevel) ./ scomb;

% mask of admissible positions for local maxima
mask = tcdf(-T, df2) < 0.05;

% clear mask borders (change border conditions for conv3fast to 'zero')
if(ip.Results.ClearMaskBorder)
    mask([1 2 end-1 end],:,:) = 0;
    mask(:,[1 2 end-1 end],:) = 0;
    mask(:,:,[1 2 end-1 end]) = 0;
end

% all local max
allMax = locmax3d(imgLoG, localMaxWindowSize, 'ClearBorder', false);

% local maxima above threshold in image domain
imgLM = allMax .* mask;

pstruct = [];
if sum(imgLM(:))~=0 % no local maxima found, likely a background image
    
    if ip.Results.RefineMaskLoG
        % -> set threshold in LoG domain
        logThreshold = min(imgLoG(imgLM~=0));
        logMask = imgLoG >= logThreshold;
        
        % combine masks
        mask = mask | logMask;
    end
    
    % re-select local maxima
    imgLM = allMax .* mask;
    clear allMax;
    
    % apply exclusion mask
    if ~isempty(ip.Results.Mask)
        imgLM(ip.Results.Mask==0) = 0;
    end
    
    lmIdx = find(imgLM~=0);
    [lmy,lmx,lmz] = ind2sub(size(vol), lmIdx);
    
    if ~isempty(lmIdx)
        % run localization on local maxima
        if ~ip.Results.FitMixtures
            pstruct = fitGaussians3D(vol, [lmx lmy lmz], A_est(lmIdx), sigma,...
                c_est(lmIdx), ip.Results.Mode, 'mask', mask, 'alpha', ip.Results.Alpha,...
                'ConfRadius', ip.Results.ConfRadius, 'WindowSize', ws);
        else
            pstruct = fitGaussianMixtures3D(vol, [lmx lmy lmz], A_est(lmIdx), sigma,...
               c_est(lmIdx), 'mask', mask, 'alpha', ip.Results.Alpha,...
               'ConfRadius', ip.Results.ConfRadius, 'WindowSize', ws,...
               'maxM', ip.Results.MaxMixtures);
        end
    
        % remove NaN values
        idx = ~isnan([pstruct.x]);
        if sum(idx)~=0
            fnames = fieldnames(pstruct);           
            for k = 1:length(fnames)
                pstruct.(fnames{k}) = pstruct.(fnames{k})(:,idx);
            end
            
            % significant amplitudes
            idx = [pstruct.hval_Ar] == 1;
            
            % eliminate duplicate positions (resulting from localization)
            if ip.Results.RemoveRedundant
                pM = [pstruct.x' pstruct.y' pstruct.z'];
                idxKD = KDTreeBallQuery(pM, pM, ip.Results.RedundancyRadius*ones(numel(pstruct.x),1));
                idxKD = idxKD(cellfun(@numel, idxKD)>1);
                for k = 1:length(idxKD);
                    RSS = pstruct.RSS(idxKD{k});
                    idx(idxKD{k}(RSS ~= min(RSS))) = 0;
                end
            end
            
            if sum(idx)>0
                fnames = fieldnames(pstruct);
                for k = 1:length(fnames)
                    pstruct.(fnames{k}) = pstruct.(fnames{k})(:,idx);
                end
                pstruct.hval_Ar = logical(pstruct.hval_Ar);
                pstruct.hval_AD = logical(pstruct.hval_AD);
                pstruct.isPSF = ~pstruct.hval_AD;
                
                % adjust mixture index if only a single component remains
                if ip.Results.FitMixtures
                    mv = 0:max(pstruct.mixtureIndex);
                    multiplicity = getMultiplicity(pstruct.mixtureIndex);
                    pstruct.mixtureIndex(ismember(pstruct.mixtureIndex, mv(multiplicity==1))) = 0;
                end
            else
                pstruct = [];
            end
            else
            pstruct = [];
        end
    end
end

if ~isempty(pstruct) && ip.Results.RefineMaskValid
    CC = bwconncomp(mask);
    labels = labelmatrix(CC);
    loclabels = labels(sub2ind(size(vol), pstruct.y_init, pstruct.x_init));
    idx = setdiff(1:CC.NumObjects, loclabels);
    CC.PixelIdxList(idx) = [];
    CC.NumObjects = length(CC.PixelIdxList);
    mask = labelmatrix(CC)~=0;
end
