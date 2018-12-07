function [GaussGrad,Gauss] = GaussFitND_gradient(fitParameters,parameters,coordList,isNormed)
%GAUSSFITND_GRADIENT returns the gradient or [residuals,Jacobian] for a N-dimensional Gaussian
%
% SYNOPSIS: [GaussGrad,Gauss] = GaussFitND_gradient(fitParameters,parameters,coordList,isNormed)
%
% INPUT fitParameters: Cell array of strings of parameters with respect to
%        which the derivative should be taken. The derivatives will be
%        taken in the order in which they appear in the list. If you set
%        fitParameters = {'X2','A','S1'}, the first column of the gradient
%        will be the derivative with respect to X2, the second colum the
%        derivative with respect to the amplitude, and the third derivative
%        with respect to sigma1.
%           -'X#' position in dimension #
%           -'A' amplitude
%           -'S#' sigma in dimension # (Sxy if S1 and S2 are identical)
%           -'B' background
%		parameters X1...n,A,S1...n,B. You have to specify the parameters
%		 whether they will be fitted or not! If S1==S2, you specify SXY
%		 instead of S1 and S2 (will result in one fewer parameters). If you
%		 want to fit multiple Gaussians, the parameters should be given in
%		 multiple rows, with the same background everywhere.
%
%       coordList: list of coordinates where the gradient should be
%        evaluated
%
%		isNormed: (opt) 1 if Gaussian(s) should each be normed to integral
%		 1, {0} otherwise (max of Gaussian will be 1).
%
% OUTPUT GaussGrad: Gradient of the gaussian with respect to the parameters
%         in 'fitParameters'
%        Gauss: Gaussian of which the gradient has been taken.
%
% REMARKS If multiple Gaussians should be calculated, the derivatives will
%          be ordered [row1][row2][...][B], i.e. with the same background
%          for all the fits.
%
% created with MATLAB ver.: 7.1.0.246 (R14) Service Pack 3 on Windows_NT
%
% created by: jdorn
% DATE: 03-Mar-2006
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

%==========================
% TEST INPUT
%==========================

% defaults
def_isNormed = 0;


% nargin
if nargin < 3 || isempty(fitParameters) || isempty(parameters) || isempty(coordList)
    error('not enough input arguments for GaussFitND_gradient!')
end

% get dimensionality
nDims = size(coordList,2);
nCoords = size(coordList,1);

% check number of parameters. It can either be d+1+d+1 or d+1+(d-1)+1,
% depending on whether we have sxy or not
nParameters = size(parameters,2);
if nParameters == 2*nDims+2
    isSxy = 0;
elseif nParameters == 2*nDims+1
    isSxy = 1;
else
    error('You have to specify all the parameters for the %i-d Gaussian!',nDims)
end
% get the number of kernels
nKernels = size(parameters,1);


% check the number of fitParameters
nFitParameters = length(fitParameters);
if ~iscell(fitParameters) || nFitParameters > nParameters
    error(['FitParameters has to be a cell array of strings',...
        'no longer than the number of parameters of the Gaussian!'])
end


% check optional input
if nargin < 4 || isempty(isNormed);
    isNormed = def_isNormed;
end

anyBackground = any(strncmpi('b',fitParameters,1));

%==========================



%==========================
% TAKE THE DERIVATIVE
%==========================

% We loop through the kernels. For every one of them, we first
% calculate the Gaussian, then take the necessary derivatives and finally,
% add them to the output.

% preassign output
GaussGrad = zeros(nCoords,(nFitParameters-anyBackground)*nKernels+anyBackground);
Gauss = zeros(nCoords,1);

for iKernel = 1:nKernels

    %---- Read data

    % read center coords
    center = parameters(iKernel,1:nDims);

    % read sigmas
    sigma = parameters(iKernel,nDims+2:2*nDims+1-isSxy);
    if isSxy % repeat first sigma twice
        sigma = sigma([1,1:end]);
    end

    % read amplitude and background
    ampAndBg = parameters(iKernel,[nDims+1,end]);



    %---- Calculate Gaussian

    % calculate Gaussian. Don't ask for normed, b/c we'll be multiplying
    % with "normed" amplitude later
    tmpGauss = GaussListND(coordList,sigma,center,0);

    % set amplitude
    tmpGauss = tmpGauss * ampAndBg(1);

    %---- Calculate Gradients

    for iParm = 1:nFitParameters

        % count variables with respect to which the gradient is taken.
        % Background will always be the last one, so we can count to
        % nFitParameters-1 before we go to the next kernel
        varCt = (iKernel-1)*(nFitParameters-1) + iParm;

        % switch according to fitParameter
        currentParm = fitParameters{iParm};
        switch currentParm(1)
            case {'X','x'}
                % derivative with respect to center coordinate
                % dG/dx = (x-x0)/sx^2 * G(x)

                currentDim = str2double(currentParm(2:end));

                GaussGrad(:,varCt) = (coordList(:,currentDim) -...
                    parameters(currentDim))...
                    / sigma(currentDim)^2 .* tmpGauss;

            case {'S','s'}
                % derivative with respect to sigma
                % dG/dsx = ((x-x0)^2 - sx^2)/sx^3 * G(sx)
                % if not normed: no -sx^2
                % if sxy: +(y-y0)^2 in the numerator, -2sxy instead of
                %   -sx.

                % take care of the special case sxy
                if strcmpi(currentParm(2:end),'xy')
                    % use sx
                    % currentDim = 1;

                    GaussGrad(:,varCt) = ...
                        ((coordList(:,1) - parameters(1)).^2 + ...
                        (coordList(:,2) - parameters(2)).^2 -...
                        sigma(1)^2*isNormed)...
                        / sigma(1)^3 .* tmpGauss;
                else
                    % get dimension
                    currentDim = str2double(currentParm(2:end));

                    GaussGrad(:,varCt) = ...
                        ((coordList(:,currentDim) - parameters(currentDim)).^2 -...
                        sigma(currentDim)^2*isNormed)...
                        / sigma(currentDim)^3 .* tmpGauss;
                end

            case {'A','a'}
                % derivative with respect to amplitude.
                % dG/da = 1/a * G(a)
                GaussGrad(:,varCt) = tmpGauss / ampAndBg(1);

            case {'B','b'}
                % derivative with respect to amplitude.
                % dG/da = 1 (irrespective of how many kernels there are)
                if iKernel < nKernels
                    % don't assign
                else
                    GaussGrad(:,end) = 1;
                end

            otherwise
                error('fitParameter ''%s'' not recognized!',currentParm)
        end % switch fitParm


    end % loop fitParameters

    % sum the Gaussians
    Gauss = Gauss + tmpGauss;

end % loop kernels

% add background only once
Gauss = Gauss + ampAndBg(2);

% the end