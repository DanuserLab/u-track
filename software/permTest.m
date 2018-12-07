function [H, pValue] = permTest(s1, s2, varargin)
% PERMTEST performs the two-sample permutation test for means (default) or medians
%
% Inputs:
%         s1, s2 : sample vectors
%
% Options: 
%          alpha : alpha value, default: 0.05.
%           tail : specifies the type of null hypothesis
%                  'both': mean(s1) != mean(s2)
%                  'right': mean(s1) > mean(s2)
%                   'left': mean(s1) < mean(s2)
%           nrep : # of permutations, default: 1900. This gives a coefficient of
%                  variation <=0.10 for alpha = 0.05. Calculated as
%                  nrep = (1-alpha)/cv^2/alpha. See [1] for details.
%
% Reference:
% [1] Efron, B. and Tibshirani, R., "An introduction to the Bootstrap," Ch. 15, 1993.
%
%  Example data from [1]
%   y = [10 27 31 40 46 50 52 104 146];
%   z = [16 23 38 94 99 141 197];
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

% Francois Aguet, 11/05/2012

ip = inputParser;
ip.addRequired('s1', @isnumeric);
ip.addRequired('s2', @isnumeric);
ip.addOptional('alpha', 0.05, @isscalar);
ip.addOptional('tail', 'both', @(x) any(strcmpi(x, {'both', 'left', 'right'})));
ip.addOptional('nrep', 1900, @isscalar);
ip.addParamValue('CmpFunction', @mean, @(x) isa(x, 'function_handle'));
ip.parse(s1, s2, varargin{:})
nrep = ip.Results.nrep;
fct = ip.Results.CmpFunction;

s1 = s1(:);
s2 = s2(:);
sAll = [s1; s2];

n1 = numel(s1);
n2 = numel(s2);
N = n1+n2;

% Calculate the number of permutations. If small, run exact test
w = warning('off', 'MATLAB:nchoosek:LargeCoefficient');
nperms = nchoosek(n1+n2, n1);
warning(w);
if nperms<=nrep % calculate all permutations
    P = false(N, nperms);
    pidx = nchoosek(1:N, n1); % returns row index of class 'sample 1'
    % convert to linear index
    pidx = pidx + repmat(N*(0:nperms-1)', [1 n1]);
    % category (1->sample1, 0->sample2) matrix for all permutations
    P(pidx) = true;
    delta = zeros(nperms,1);
    for i = 1:nperms
        delta(i) = fct(sAll(P(:,i))) - fct(sAll(~P(:,i)));
    end
    ns = nperms;
else % compute 'nrep' random permutations
    delta = zeros(nrep,1);
    for i = 1:nrep
        idx = randperm(N); % calculate random permutation of the samples
        delta(i) = fct(sAll(idx(1:n1))) - fct(sAll(idx(n1+1:end)));
    end
    ns = nrep;
end

deltaRef = fct(s1)-fct(s2);

switch ip.Results.tail
    case 'both'
        pValue = sum(abs(delta)>=abs(deltaRef))/ns;
    case 'right'
        pValue = sum(delta>=deltaRef)/ns;
    case 'left'
        pValue = sum(delta<=deltaRef)/ns;
end

H = pValue <= ip.Results.alpha;
