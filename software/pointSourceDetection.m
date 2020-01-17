%[pstruct, mask, imgLM, imgLoG] = pointSourceDetection(img, sigma, mode)
%
% Inputs :   
%                 img : input image
%               sigma : standard deviation of the Gaussian PSF
%
% Options (as 'specifier'-value pairs): 
%
%              'mode' : parameters to estimate. Default: 'xyAc'.
%             'alpha' : alpha value used in the statistical tests. Default: 0.05.
%              'mask' : mask of pixels (i.e., cell mask) to include in the detection. Default: all.
%       'FitMixtures' : true|{false}. Toggles mixture-model fitting.
%       'MaxMixtures' : maximum number of mixtures to fit. Default: 5.
%   'RemoveRedundant' : {true}|false. Discard localizations that coincide within 'RedundancyRadius'.
%  'RedundancyRadius' : Radius for filtering out redundant localizatios. Default: 0.25
%         'Prefilter' : {true}|false. Prefilter to calculate mask of significant pixels.
%     'RefineMaskLoG' : {true}|false. Apply threshold to LoG-filtered img to refine mask of significant pixels.
%   'RefineMaskValid' : {true}|false. Return only mask regions where a significant signal was localized.
%        'ConfRadius' : Confidence radius for positions, beyond which the fit is rejected. Default: 2*sigma
%        'WindowSize' : Window size for the fit. Default: 4*sigma, i.e., [-4*sigma ... 4*sigma]^2
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

% Francois Aguet, April 2011 (last modified: 09/30/2013)

function [pstruct, mask, imgLM, imgLoG] = pointSourceDetection(img, sigma, varargin)

% Parse inputs
ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('img', @isnumeric);
ip.addRequired('sigma', @isscalar);
ip.addParamValue('Mode', 'xyAc', @ischar);
ip.addParamValue('Alpha', 0.05, @isscalar);
ip.addParamValue('Mask', [], @(x) isnumeric(x) || islogical(x));
ip.addParamValue('FitMixtures', false, @islogical);
ip.addParamValue('MaxMixtures', 5, @isposint);
ip.addParamValue('RemoveRedundant', true, @islogical);
ip.addParamValue('RedundancyRadius', 0.25, @isscalar);
ip.addParamValue('Prefilter', true, @islogical);
ip.addParamValue('RefineMaskLoG', true, @islogical);
ip.addParamValue('RefineMaskValid', true, @islogical);
ip.addParamValue('ConfRadius', []); % Default: 2*sigma, see fitGaussians2D.
ip.addParamValue('WindowSize', []); % Default: 4*sigma, see fitGaussians2D.
ip.KeepUnmatched = true;
ip.parse(img, sigma, varargin{:});
mode = ip.Results.Mode;
alpha = ip.Results.Alpha;

if ~isa(img, 'double')
    img = double(img);
end

% Gaussian kernel
w = ceil(4*sigma);
x = -w:w;
g = exp(-x.^2/(2*sigma^2));
u = ones(1,length(x));

% convolutions
imgXT = padarrayXT(img, [w w], 'symmetric');
fg = conv2(g', g, imgXT, 'valid');
fu = conv2(u', u, imgXT, 'valid');
fu2 = conv2(u', u, imgXT.^2, 'valid');

% Laplacian of Gaussian
gx2 = g.*x.^2;
imgLoG = 2*fg/sigma^2 - (conv2(g, gx2, imgXT, 'valid')+conv2(gx2, g, imgXT, 'valid'))/sigma^4;
imgLoG = imgLoG / (2*pi*sigma^2);

% 2-D kernel
g = g'*g;
n = numel(g);
gsum = sum(g(:));
g2sum = sum(g(:).^2);

% solution to linear system
A_est = (fg - gsum*fu/n) / (g2sum - gsum^2/n);
c_est = (fu - A_est*gsum)/n;

if ip.Results.Prefilter
    J = [g(:) ones(n,1)]; % g_dA g_dc
    C = inv(J'*J);
    
    f_c = fu2 - 2*c_est.*fu + n*c_est.^2; % f-c
    RSS = A_est.^2*g2sum - 2*A_est.*(fg - c_est*gsum) + f_c;
    RSS(RSS<0) = 0; % negative numbers may result from machine epsilon/roundoff precision
    sigma_e2 = RSS/(n-3);
    
    sigma_A = sqrt(sigma_e2*C(1,1));
    
    % standard deviation of residuals
    sigma_res = sqrt(RSS/(n-1));
    
    kLevel = norminv(1-alpha/2.0, 0, 1);
    
    SE_sigma_c = sigma_res/sqrt(2*(n-1)) * kLevel;
    df2 = (n-1) * (sigma_A.^2 + SE_sigma_c.^2).^2 ./ (sigma_A.^4 + SE_sigma_c.^4);
    scomb = sqrt((sigma_A.^2 + SE_sigma_c.^2)/n);
    T = (A_est - sigma_res*kLevel) ./ scomb;
    pval = tcdf(-T, df2);
    
    % mask of admissible positions for local maxima
    mask = pval < 0.05;
else
    mask = true(size(img));
end

% all local max
allMax = locmax2d(imgLoG, 2*ceil(sigma)+1);

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
    
    % apply exclusion mask
    if ~isempty(ip.Results.Mask)
        imgLM(ip.Results.Mask==0) = 0;
    end
    
    
    [lmy, lmx] = find(imgLM~=0);
    lmIdx = sub2ind(size(img), lmy, lmx);    
    
    if ~isempty(lmIdx)
        % run localization on local maxima
        if ~ip.Results.FitMixtures            
            pstruct = fitGaussians2D(img, lmx, lmy, A_est(lmIdx), sigma*ones(1,length(lmIdx)),...
                c_est(lmIdx), mode, 'mask', mask, 'alpha', alpha,...
                'ConfRadius', ip.Results.ConfRadius, 'WindowSize', ip.Results.WindowSize);
        else
            pstruct = fitGaussianMixtures2D(img, lmx, lmy, A_est(lmIdx), sigma*ones(1,length(lmIdx)),...
                c_est(lmIdx), 'mask', mask, 'alpha', alpha, 'maxM', ip.Results.MaxMixtures);
        end
        
        % remove NaN values
        idx = ~isnan([pstruct.x]);
        if sum(idx)~=0
            fnames = fieldnames(pstruct);
            for k = 1:length(fnames)
                pstruct.(fnames{k}) = pstruct.(fnames{k})(idx);
            end
            
            % significant amplitudes
            idx = [pstruct.hval_Ar] == 1;
            
            % eliminate duplicate positions (resulting from localization)
            if ip.Results.RemoveRedundant
                pM = [pstruct.x' pstruct.y'];
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
                    pstruct.(fnames{k}) = pstruct.(fnames{k})(idx);
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
    loclabels = labels(sub2ind(size(img), pstruct.y_init, pstruct.x_init));
    idx = setdiff(1:CC.NumObjects, loclabels);
    CC.PixelIdxList(idx) = [];
    CC.NumObjects = length(CC.PixelIdxList);
    mask = labelmatrix(CC)~=0;
end
