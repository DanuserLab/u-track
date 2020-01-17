% pStruct = fitGaussianMixtures2D(vol, X, A, sigma, c, mode, varargin) fits 
% 3D Gaussian mixtures to the input volume at the specified locations.
%
% Inputs:   
%           vol : input volume
%             X : Nx3 matrix of initial (or fixed) (x,y,z)-positions
%             A : initial (or fixed) amplitudes
%         sigma : initial (or fixed) Gaussian PSF standard deviations
%             c : initial (or fixed) background intensities
%
% Optional inputs : ('Mask', mask) pair with a mask of spot locations
%
% Output: 
%         pStruct: structure with fields:
%                  x : estimated x-positions
%                  y : estimated y-positions
%                  z : estimated z-positions
%                  A : estimated amplitudes
%                  s : estimated standard deviations of the PSF
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
%            pval_Ar : p-value of an amplitude vs. background noise test (p > 0.05 -> significant amplitude)
%
%
% Usage for a volume with know spot locations (mask) and fixed sigma:
% fitGaussianMixtures3D(vol, X, A, sigma, 'mask', mask);
%
% See also fitGaussian3D, fitGaussians3D
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

% Francois Aguet, March 28 2011 (last modified: Feb 5 2013)

function pStruct = fitGaussianMixtures3D(vol, X, A, sigma, c, varargin)

% Parse inputs
ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('vol', @isnumeric);
ip.addRequired('X');
ip.addRequired('A');
ip.addRequired('sigma');
ip.addRequired('c');
ip.addParamValue('Alpha', 0.05, @isscalar);
ip.addParamValue('AlphaT', 0.05, @isscalar);
ip.addParamValue('Mask', [], @islogical);
ip.addParamValue('ConfRadius', []);
ip.addParamValue('WindowSize', []);
ip.addParamValue('maxM', 5, @isscalar);
ip.parse(vol, X, A, sigma, c, varargin{:});

np = size(X,1);
sigma = ip.Results.sigma;
if numel(sigma)==1
    sigma = sigma*ones(np,2);
elseif numel(sigma)==2
    sigma = repmat(sigma(:)', [np 1]);
end
mode = 'xyzAc'; % sigma values are fixed to improve stability of the fit
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

% initialize pStruct arrays
pStruct.x = cell(1,np);
pStruct.y = cell(1,np);
pStruct.z = cell(1,np);
pStruct.A = cell(1,np);
pStruct.s = cell(2,np);
pStruct.c = cell(1,np);

pStruct.x_pstd = cell(1,np);
pStruct.y_pstd = cell(1,np);
pStruct.z_pstd = cell(1,np);
pStruct.A_pstd = cell(1,np);
pStruct.s_pstd = cell(2,np);
pStruct.c_pstd = cell(1,np);

pStruct.x_init = cell(1,np);
pStruct.y_init = cell(1,np);
pStruct.z_init = cell(1,np);

pStruct.sigma_r = cell(1,np);
pStruct.SE_sigma_r = cell(1,np);
pStruct.RSS = cell(1,np);

pStruct.pval_Ar = cell(1,np);

pStruct.hval_AD = cell(1,np);
pStruct.hval_Ar = cell(1,np);

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

mixtureIndex = 1;
% loop through initial points
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
    
    window = vol(ya,xa,za);
    % set any other mask components to NaN
    window(maskWindow~=0) = NaN;
    npx = sum(isfinite(window(:)));
    
    if npx >= 20 % only perform fit if window contains sufficient data points

        % initial fit with a single Gaussian 
        % Notation: reduced model: '_r', full model: '_f'
        [prm_f, prmStd_f, ~, res_f] = fitGaussian3D(window, [X(p,1)-xi(p)+ox X(p,2)-yi(p)+oy X(p,3)-zi(p)+oz A(p) sigma(p,:) c(p)], mode);

        % update standard deviations for xyzAc
        prmStd_f = [prmStd_f(1:4) 0 0 prmStd_f(5)];
        RSS_r = res_f.RSS;
        
        p_r = 5; % #parameters in the model (x,y,z,A,c)     
        i = 1; % iteration
        
        pval = 0;
        validBounds = true;
        while i<ip.Results.maxM && pval<0.05 && validBounds % models are significantly different
                
            i = i+1;
            prm_r = prm_f;
            prmStd_r = prmStd_f;
            res_r = res_f;
            
            % expanded model
            % new component: initial values given by max. residual point
            [maxRes, idx] = max(res_r.data(:));
            [y0, x0, z0] = ind2sub(size(window), idx);
            
            initV = [x0-1 y0-1 z0-1 maxRes prm_r];
            %initV(4:4:end-3) = sum(prm_r(4:4:end-3))/i; % may work better in some cases
            [prm_f, prmStd_f, ~, res_f] = fitGaussianMixture3D(window, initV, mode);
            
            RSS_f = res_f.RSS;
            p_f = p_r + 4; % 4 parameters (x,y,z,A) added at every iteration
            
            % test statistic (F-test)
            T = (RSS_r-RSS_f)/RSS_f * (npx-p_f-1)/(p_f-p_r);
            pval = 1-fcdf(T,p_f-p_r,npx-p_f-1);
            
            % update reduced model
            p_r = p_f;
            RSS_r = RSS_f;

            % restrict radius; otherwise neighboring signals are considered part of mixture
            dx = prm_f(1:4:end-3)-ox;
            dy = prm_f(2:4:end-3)-oy;
            dz = prm_f(3:4:end-3)-oz;
            if min(dx)<-w2(1) || max(dx)>w2(1) || min(dy)<-w2(1) || max(dy)>w2(1) || min(dz)<-w2(2) || max(dz)>w2(2)
                validBounds = false;
            end
        end
        ng = i-1; % # gaussians in final model
        
        % sigma, c are the same for each mixture
        x_est = prm_r(1:4:end-3)-ox;
        y_est = prm_r(2:4:end-3)-oy;
        z_est = prm_r(3:4:end-3)-oz;
        A_est = prm_r(4:4:end-3);
        
        % exclude points where localization failed
        if ng>1 || (x_est > -w2(1) && x_est < w2(1) && y_est > -w2(1) && y_est < w2(1) &&...
                z_est > -w2(2) && z_est < w2(2) && A_est<2*diff(iRange))
            pStruct.x{p} = xi(p) + x_est;
            pStruct.y{p} = yi(p) + y_est;
            pStruct.z{p} = zi(p) + z_est;
            pStruct.A{p} = A_est;
            % sigma and background offset are identical for all mixture components
            pStruct.s{p} = repmat(prm_r([end-2 end-2])', [1 ng]);
            pStruct.c{p} = repmat(prm_r(end), [1 ng]);
            
            pStruct.x_pstd{p} = prmStd_r(1:4:end-3);
            pStruct.y_pstd{p} = prmStd_r(2:4:end-3);
            pStruct.z_pstd{p} = prmStd_r(3:4:end-3);
            pStruct.A_pstd{p} = prmStd_r(4:4:end-3);
            pStruct.s_pstd{p} = repmat(prmStd_r([end-2 end-1])', [1 ng]);
            pStruct.c_pstd{p} = repmat(prmStd_r(end), [1 ng]);
            
            pStruct.x_init{p} = repmat(xi(p), [1 ng]);
            pStruct.y_init{p} = repmat(yi(p), [1 ng]);
            pStruct.z_init{p} = repmat(zi(p), [1 ng]);
            
            pStruct.sigma_r{p} = repmat(res_r.std, [1 ng]);
            pStruct.RSS{p} = repmat(res_r.RSS, [1 ng]);
            
            SE_sigma_r = res_r.std/sqrt(2*(npx-1));
            pStruct.SE_sigma_r{p} = repmat(SE_sigma_r, [1 ng]);
            SE_sigma_r = SE_sigma_r * kLevel;
            
            pStruct.hval_AD{p} = repmat(res_r.hAD, [1 ng]);
            
            % H0: A <= k*sigma_r
            % H1: A > k*sigma_r
            for i = 1:ng
                sigma_A = pStruct.A_pstd{p}(i);
                A_est = pStruct.A{p}(i);
                df2 = (npx-1) * (sigma_A.^2 + SE_sigma_r.^2).^2 ./ (sigma_A.^4 + SE_sigma_r.^4);
                scomb = sqrt((sigma_A.^2 + SE_sigma_r.^2)/npx);
                T = (A_est - res_r.std*kLevel) ./ scomb;
                % 1-sided t-test: A_est must be greater than k*sigma_r
                pStruct.pval_Ar{p}(i) = tcdf(-T, df2);
                pStruct.hval_Ar{p}(i) = pStruct.pval_Ar{p}(i) < ip.Results.AlphaT;
                %pStruct.mask_Ar{p}(i) = sum(A_est*g>res_r.std*kLevel); % # significant pixels
            end
            if ng>1
                pStruct.mixtureIndex{p} = mixtureIndex*ones(1,ng);
                mixtureIndex = mixtureIndex+1;
            else
                pStruct.mixtureIndex{p} = 0;
            end
        end
    end
end

% concatenate cell arrays
fnames = fieldnames(pStruct);
for f = 1:numel(fnames)
    pStruct.(fnames{f}) = [pStruct.(fnames{f}){:}];
end
