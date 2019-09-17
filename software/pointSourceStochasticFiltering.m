function [mask, dummy, imgLoG,imgLoGScaled,testMap] = pointSourceStochasticFiltering(vol, sigma, varargin)
% P. Roudot 2016. Credit to F. Aguet 2013
%
% Copyright (C) 2019, Danuser Lab - UTSouthwestern 
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

% Parse inputs
ip = inputParser;
ip.CaseSensitive = false;
ip.KeepUnmatched = true;
ip.addRequired('vol', @isnumeric);
ip.addRequired('sigma', @isnumeric);
% ip.addParamValue('Mode', 'xyzAc', @ischar);
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

tic;
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

gkernel=gx*gx'*reshape(gz,[1 1 numel(gz)]);
%imgLoGScaled=imgLoG./nansum(abs(gkernel(:)));
imgLoGScaled=imgLoG;%./nansum((gkernel(:).^2).^0.5;

clear fgx2 fgy2 fgz2;
disp('Convolutions');toc;

tic;
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
disp('Pixelic fit');toc;

tic;
f_c = fu2 - 2*c_est.*fu + n*c_est.^2; % f-c
RSS = A_est.^2*g2sum - 2*A_est.*(fg - c_est*gsum) + f_c;
RSS=RSS-(A_est*gsum+n*c_est-fu)/n;
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

testMap=A_est+c_est;

% clear mask borders (change border conditions for conv3fast to 'zero')
mask([1 2 end-1 end],:,:) = 0;
mask(:,[1 2 end-1 end],:) = 0;
mask(:,:,[1 2 end-1 end]) = 0;
disp('Stochastic thresholding');toc;

% all local max
% allMax = locmax3d(imgLoG, localMaxWindowSize, 'ClearBorder', false);
% local maxima above threshold in image domain

dummy = mask;