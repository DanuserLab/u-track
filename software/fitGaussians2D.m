% pStruct = fitGaussians2D(img, x, y, A, sigma, c, mode, varargin)
%
% Inputs:   img : input image
%             x : initial (or fixed) x-positions
%             y : initial (or fixed) y-positions
%             A : initial (or fixed) amplitudes
%             s : initial (or fixed) Gaussian PSF standard deviations
%             c : initial (or fixed) background intensities
%
% Options:
%          mode : selector for optimization parameters, any of 'xyAsc'. Default: 'xyAc'
%
% Options ('specifier', value):
%          'Mask' : mask of spot locations
%    'ConfRadius' : Confinement radius for valid fits. Default: ceil(2*sigma)
%    'WindowSize' : Size of the support used for the fit, specified as half-width;
%                   i.e., for a window of 15x15, enter 7. Default: ceil(4*sigma)
%
% Output: pStruct: structure with fields:
%                  x : estimated x-positions
%                  y : estimated y-positions
%                  A : estimated amplitudes
%                  s : estimated standard deviations of the PSF
%                  c : estimated background intensities
%
%             x_pstd : standard deviations, estimated by error propagation
%             y_pstd : "
%             A_pstd : "
%             s_pstd : "
%             c_pstd : "
%            sigma_r : standard deviation of the background (residual)
%         SE_sigma_r : standard error of sigma_r
%            pval_Ar : p-value of an amplitude vs. background noise test (p < 0.05 -> significant amplitude)
%
%
% Usage for a single-channel image with mask and fixed sigma:
% fitGaussians2D(img, x, y, A, sigma, c, 'xyAc', 'mask', mask);
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

% Francois Aguet, March 28 2011 (last updated: Sep 30 2013)

function pStruct = fitGaussians2D(img, x, y, A, sigma, c, varargin)

% Parse inputs
ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('img', @isnumeric);
ip.addRequired('x');
ip.addRequired('y');
ip.addRequired('A');
ip.addRequired('sigma');
ip.addRequired('c');
ip.addOptional('mode', 'xyAc', @ischar);
ip.addParamValue('Alpha', 0.05, @isscalar);
ip.addParamValue('AlphaT', 0.05, @isscalar);
ip.addParamValue('Mask', [], @islogical);
ip.addParamValue('ConfRadius', []);
ip.addParamValue('WindowSize', []);
ip.parse(img, x, y, A, sigma, c, varargin{:});

np = length(x);
sigma = ip.Results.sigma;
if numel(sigma)==1
    sigma = sigma*ones(1,np);
end
mode = ip.Results.mode;
if ~isempty(ip.Results.Mask)
    labels = bwlabel(ip.Results.Mask);
else
    labels = zeros(size(img));
end

pStruct = struct('x', [], 'y', [], 'A', [], 's', [], 'c', [],...
    'x_pstd', [], 'y_pstd', [], 'A_pstd', [], 's_pstd', [], 'c_pstd', [],...
    'x_init', [], 'y_init', [],...
    'sigma_r', [], 'SE_sigma_r', [], 'RSS', [], 'pval_Ar', [], 'mask_Ar', [], 'hval_Ar', [], 'hval_AD', []);

xi = round(x);
yi = round(y);
[ny,nx] = size(img);

kLevel = norminv(1-ip.Results.Alpha/2.0, 0, 1); % ~2 std above background

iRange = [min(img(:)) max(img(:))];

estIdx = regexpi('xyAsc', ['[' mode ']']);


% initialize pStruct arrays
pStruct.x = NaN(1,np);
pStruct.y = NaN(1,np);
pStruct.A = NaN(1,np);
pStruct.s = NaN(1,np);
pStruct.c = NaN(1,np);

pStruct.x_pstd = NaN(1,np);
pStruct.y_pstd = NaN(1,np);
pStruct.A_pstd = NaN(1,np);
pStruct.s_pstd = NaN(1,np);
pStruct.c_pstd = NaN(1,np);

pStruct.x_init = reshape(xi, [1 np]);
pStruct.y_init = reshape(yi, [1 np]);

pStruct.sigma_r = NaN(1,np);
pStruct.SE_sigma_r = NaN(1,np);
pStruct.RSS = NaN(1,np);

pStruct.pval_Ar = NaN(1,np);
pStruct.mask_Ar = zeros(1,np);

pStruct.hval_AD = false(1,np);
pStruct.hval_Ar = false(1,np);


sigma_max = max(sigma);
w2 = ip.Results.ConfRadius;
if isempty(w2)
    w2 = ceil(2*sigma_max);
end
w4 = ip.Results.WindowSize;
if isempty(w4)
    w4 = ceil(4*sigma_max);
end

% for background estimation
if isempty(c)
    % mask template: ring with inner radius w3, outer radius w4
    [xm,ym] = meshgrid(-w4:w4);
    r = sqrt(xm.^2+ym.^2);
    annularMask = zeros(size(r));
    annularMask(r<=ceil(4*sigma_max) & r>=ceil(3*sigma_max)) = 1;
end

g = exp(-(-w4:w4).^2/(2*sigma_max^2));
g = g'*g;
g = g(:);

T = zeros(1,np);
df2 = zeros(1,np);
for p = 1:np
    
    % ignore points in border
    if (xi(p)>w4 && xi(p)<=nx-w4 && yi(p)>w4 && yi(p)<=ny-w4)
        
        % label mask
        maskWindow = labels(yi(p)-w4:yi(p)+w4, xi(p)-w4:xi(p)+w4);
        maskWindow(maskWindow==maskWindow(w4+1,w4+1)) = 0;
        
        window = img(yi(p)-w4:yi(p)+w4, xi(p)-w4:xi(p)+w4);

        % estimate background        
        if isempty(c)
            cmask = annularMask;
            cmask(maskWindow~=0) = 0;
            c_init = mean(window(cmask==1));
        else
            c_init = c(p);
        end
        
        % set any other components to NaN
        window(maskWindow~=0) = NaN;
        npx = sum(isfinite(window(:)));
        
        if npx >= 10 % only perform fit if window contains sufficient data points
        
            % fit
            if isempty(A)
                A_init = max(window(:))-c_init;
            else
                A_init = A(p);
            end
            
            [prm, prmStd, ~, res] = fitGaussian2D(window, [x(p)-xi(p) y(p)-yi(p) A_init sigma(p) c_init], mode);
            
            dx = prm(1);
            dy = prm(2);
            
            % exclude points where localization failed
            if (dx > -w2 && dx < w2 && dy > -w2 && dy < w2 && prm(3)<2*diff(iRange))
                
                pStruct.x(p) = xi(p) + dx;
                pStruct.y(p) = yi(p) + dy;
                pStruct.A(p) = prm(3);
                pStruct.s(p) = prm(4);
                pStruct.c(p) = prm(5);
                
                stdVect = zeros(1,5);
                stdVect(estIdx) = prmStd;
                
                pStruct.x_pstd(p) = stdVect(1);
                pStruct.y_pstd(p) = stdVect(2);
                pStruct.A_pstd(p) = stdVect(3);
                pStruct.s_pstd(p) = stdVect(4);
                pStruct.c_pstd(p) = stdVect(5);
                
                pStruct.sigma_r(p) = res.std;
                pStruct.RSS(p) = res.RSS;
                
                pStruct.SE_sigma_r(p) = res.std/sqrt(2*(npx-1));
                SE_sigma_r = pStruct.SE_sigma_r(p) * kLevel;
                
                pStruct.hval_AD(p) = res.hAD;
                
                % H0: A <= k*sigma_r
                % H1: A > k*sigma_r
                sigma_A = stdVect(3);
                A_est = prm(3);
                df2(p) = (npx-1) * (sigma_A.^2 + SE_sigma_r.^2).^2 ./ (sigma_A.^4 + SE_sigma_r.^4);
                scomb = sqrt((sigma_A.^2 + SE_sigma_r.^2)/npx);
                T(p) = (A_est - res.std*kLevel) ./ scomb;
                pStruct.mask_Ar(p) = sum(A_est*g>res.std*kLevel);
            end
        end
    end
end
% 1-sided t-test: A_est must be greater than k*sigma_r
pStruct.pval_Ar = tcdf(-T, df2);
pStruct.hval_Ar = pStruct.pval_Ar < ip.Results.AlphaT;
