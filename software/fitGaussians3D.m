% pStruct = fitGaussians3D(vol, X, A, sigma, c, mode, varargin) fits
% 3D Gaussians to the input volume at the specified locations.
%
% Inputs:   vol : input volume
%             X : Nx3 matrix of initial (or fixed) (x,y,z)-positions
%             A : initial (or fixed) amplitudes
%         sigma : initial (or fixed) Gaussian PSF standard deviation. If z-sampling is different
%                 from x,y sampling, this should be a 2-element vector
%             c : initial (or fixed) local background intensities
%
% Options:
%          mode : selector for optimization parameters, any of 'xyzAsrc'. Default: 'xyzAc'
%                 's' selects the (x,y) s.d., and 'r' the (z) s.d. 
%
% Options ('specifier', value):
%        'Mask' : mask of spot locations
%
% Output: pStruct: structure with fields:
%                  x : estimated x-positions
%                  y : estimated y-positions
%                  z : estimated z-positions
%                  A : estimated amplitudes
%                  s : estimated standard deviations of the PSF
%                      2xN matrix, with the 1st row containing (x,y) s.d. and 
%                      the 2nd row containing the (z) s.d.
%                  c : estimated background intensities
%
%             x_pstd : standard deviations, estimated by error propagation
%             y_pstd : "
%             z_pstd : "
%             A_pstd : "
%             s_pstd : "
%             c_pstd : "
%            sigma_r : standard deviation of the background (residual)
%         SE_sigma_r : standard error of sigma_r
%            pval_Ar : p-value of an amplitude vs. background noise test (p < 0.05 -> significant amplitude)
%
% Usage for a volume with know spot locations (mask) and fixed sigma:
% fitGaussians3D(vol, X, A, sigma, c, 'xyzAc', 'mask', mask);
%
% See also fitGaussian3D
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

% Francois Aguet, August 2013 (last updated: 10/27/2013)

function pStruct = fitGaussians3D(vol, X, A, sigma, c, varargin)

% Parse inputs
ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('vol', @isnumeric);
ip.addRequired('X');
ip.addRequired('A');
ip.addRequired('sigma');
ip.addRequired('c');
ip.addOptional('mode', 'xyzAc', @ischar);
ip.addParamValue('Alpha', 0.05, @isscalar);
ip.addParamValue('AlphaT', 0.05, @isscalar);
ip.addParamValue('Mask', [], @islogical);
ip.addParamValue('ConfRadius', []);
ip.addParamValue('WindowSize', []);
ip.parse(vol, X, A, sigma, c, varargin{:});

np = size(X,1);
sigma = ip.Results.sigma;
if numel(sigma)==1
    sigma = sigma*ones(np,2);
elseif numel(sigma)==2
    sigma = repmat(sigma(:)', [np 1]);
end
mode = ip.Results.mode;
if ~isempty(ip.Results.Mask)
    labels = double(labelmatrix(bwconncomp(ip.Results.Mask)));
else
    labels = zeros(size(vol));
end

pStruct = struct('x', [], 'y', [], 'z', [], 'A', [], 's', [], 'c', [],...
    'x_pstd', [], 'y_pstd', [], 'z_pstd', [], 'A_pstd', [], 's_pstd', [], 'c_pstd', [],...
    'x_init', [], 'y_init', [], 'z_init', [],...
    'sigma_r', [], 'SE_sigma_r', [], 'RSS', [], 'pval_Ar', [], 'hval_Ar', [], 'hval_AD', []);

[ny,nx,nz] = size(vol);
roundConstr = @(x,N) max(min(round(x),N),1);
xi = roundConstr(X(:,1), nx);
yi = roundConstr(X(:,2), ny);
zi = roundConstr(X(:,3), nz);

kLevel = norminv(1-ip.Results.Alpha/2.0, 0, 1); % ~2 std above background

iRange = [min(vol(:)) max(vol(:))];

estIdx = regexpi('xyzAsrc', ['[' mode ']']);


% initialize pStruct arrays
pStruct.x = NaN(1,np);
pStruct.y = NaN(1,np);
pStruct.z = NaN(1,np);
pStruct.A = NaN(1,np);
pStruct.s = NaN(2,np);
pStruct.c = NaN(1,np);

pStruct.x_pstd = NaN(1,np);
pStruct.y_pstd = NaN(1,np);
pStruct.z_pstd = NaN(1,np);
pStruct.A_pstd = NaN(1,np);
pStruct.s_pstd = NaN(2,np);
pStruct.c_pstd = NaN(1,np);

pStruct.x_init = reshape(xi, [1 np]);
pStruct.y_init = reshape(yi, [1 np]);
pStruct.z_init = reshape(zi, [1 np]);

pStruct.sigma_r = NaN(1,np);
pStruct.SE_sigma_r = NaN(1,np);
pStruct.RSS = NaN(1,np);

pStruct.pval_Ar = NaN(1,np);

pStruct.hval_AD = false(1,np);
pStruct.hval_Ar = false(1,np);

% if different sigma values are passed, use largest for filters
sigma_max = max(sigma,[],1);
w2 = ip.Results.ConfRadius;
if isempty(w2)
    w2 = ceil(2*sigma_max);
elseif numel(w2)==1
    w2 = [w2 w2];
end
ws = ip.Results.WindowSize;
if isempty(ws)
    ws = ceil(2*sigma_max);
elseif numel(ws)==1
    ws = [ws ws];
end

T = zeros(1,np);
df2 = zeros(1,np);
for p = 1:np
    
    % window boundaries
    xa = max(1,xi(p)-ws(1)):min(nx,xi(p)+ws(1));
    ya = max(1,yi(p)-ws(1)):min(ny,yi(p)+ws(1));
    za = max(1,zi(p)-ws(2)):min(nz,zi(p)+ws(2));
    
    % relative coordinates of (xi, yi, zi) in window. Origin at (0,0,0)
    ox = xi(p)-xa(1);
    oy = yi(p)-ya(1);
    oz = zi(p)-za(1);
    
    % label mask
    maskWindow = labels(ya, xa, za);
    maskWindow(maskWindow==maskWindow(oy+1,ox+1,oz+1)) = 0;
    
    window = vol(ya, xa, za);
    % set any other components to NaN
    window(maskWindow~=0) = NaN;
    npx = sum(isfinite(window(:)));
    
    if npx >= 20 % only perform fit if window contains sufficient data points
        
        % fit
        [prm, prmStd, ~, res] = fitGaussian3D(window, [X(p,1)-xi(p)+ox X(p,2)-yi(p)+oy X(p,3)-zi(p)+oz A(p) sigma(p,:) c(p)], mode);
        
        dx = prm(1)-ox;
        dy = prm(2)-oy;
        dz = prm(3)-oz;
        
        % exclude points where localization failed
        if (dx > -w2(1) && dx < w2(1) && dy > -w2(1) && dy < w2(1) && dz > -w2(2) && dz < w2(2) && prm(4)<2*diff(iRange))
            
            pStruct.x(p) = xi(p) + dx;
            pStruct.y(p) = yi(p) + dy;
            pStruct.z(p) = zi(p) + dz;
            pStruct.s(:,p) = prm(5:6);
            pStruct.A(p) = prm(4);
            pStruct.c(p) = prm(7);
            
            stdVect = zeros(1,7);
            stdVect(estIdx) = prmStd;
            
            pStruct.x_pstd(p) = stdVect(1);
            pStruct.y_pstd(p) = stdVect(2);
            pStruct.z_pstd(p) = stdVect(3);
            pStruct.A_pstd(p) = stdVect(4);
            pStruct.s_pstd(:,p) = stdVect(5:6);
            pStruct.c_pstd(p) = stdVect(7);
            
            pStruct.sigma_r(p) = res.std;
            pStruct.RSS(p) = res.RSS;
            
            pStruct.SE_sigma_r(p) = res.std/sqrt(2*(npx-1));
            SE_sigma_r = pStruct.SE_sigma_r(p) * kLevel;
            
            pStruct.hval_AD(p) = res.hAD;
            
            % H0: A <= k*sigma_r
            % H1: A > k*sigma_r
            sigma_A = stdVect(4);
            A_est = prm(4);
            df2(p) = (npx-1) * (sigma_A.^2 + SE_sigma_r.^2).^2 ./ (sigma_A.^4 + SE_sigma_r.^4);
            scomb = sqrt((sigma_A.^2 + SE_sigma_r.^2)/npx);
            T(p) = (A_est - res.std*kLevel) ./ scomb;
        end
    end
end
% 1-sided t-test: A_est must be greater than k*sigma_r
pStruct.pval_Ar = tcdf(-T, df2);
pStruct.hval_Ar = pStruct.pval_Ar < ip.Results.AlphaT;
