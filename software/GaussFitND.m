function [parameters,sigmaParameters,Q,chiSquared,degreesOfFreedom,...
    residualImage,resAndGauss] = GaussFitND(intensities, coordList, ...
    fitParameters, parameters, isNormed, fitOptions)
%GAUSSFITND fits one or several N-D Gaussians to an input intensity distribution
%
% SYNOPSIS: [parameters,sigmaParameters,Q,chiSquared,degreesOfFreedom,...
%    residualImage,resAndGauss] = GaussFitND(intensities, coordList, ...
%    fitParameters, parameters, isNormed, fitOptions)
%
% INPUT intensities: Vector with the image intensities. Alternatively, you
%               can supply the n-d image directly.
%		coordList: nCoord-by-nDims list of coordinates corresponding to
%               every intensity in the image list. If you supplied the
%               image directly, you can pass empty here and GaussFitND will
%               calculate the coordinateList for you.
%		fitParameters: Cell array of strings that lists the parameters that
%               should be fitted
%                -'X#' position in dimension #
%                -'A' amplitude
%                -'S#' sigma in dimension # (Sxy if S1 and S2 are identical)
%                -'B' background (has to be the last parameter!)
%		 parameters: All the parameters (including initial guesses)
%               necessary to describe the N-D Gaussian:
%               [X1...n,A,S1...n,B]. If S1==S2, you specify SXY instead of
%               S1 and S2 (will result in one parameter less).If you have
%               no idea about the value of a specific parameter, supply
%               NaN, and the code will estimate the parameter.
%               If you want to jointly fit multiple kernels, supply
%               multiple rows of parameters (the background is assumed to
%               be constant, though) If you do so, you cannot enter
%               NaNs as parameters, because the estimation will be
%               meaningless.
%		 isNormed: (opt) Whether the Gaussian should be normed to an
%               integral of 1. [{0}/1]. Only relevant for sigma-fitting,
%               when the input is normed already, and thus the amplitude is
%               not a fit-parameter.
%        fitOptions: (opt) Options for lsqnonlin. It is not recommended
%               that you set 'Jacobian' to 'off'
%               Default options
%                   TolFun 1e-20
%                   TolX   1e-3
%                   Display 'off'
%
% OUTPUT parameters: Descriptors of the Gaussian(s) in the same order as
%               the initialGuess
%        sigmaParameters: uncertainty in the parameters
%                           sqrt(diag(Q)*chiSquared)
%        Q : covariance matrix of the unknowns
%        chiSquared : variance of the residuals
%        degreesOfFreedom : degrees of freedom
%        residualImage    : n-d matrix of residuals (NaN where there is no
%                           intensity given)
%
%
% created with MATLAB ver.: 7.1.0.246 (R14) Service Pack 3 on Windows_NT
%
% created by: jdorn
% DATE: 04-Mar-2006
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



%============================
% TEST INPUT
%============================

% defaults
def_isNormed = 0;
def_fitOptions = optimset('Display','off','Jacobian','on','TolFun',1e-20,'TolX',1e-3);
verbose = 0; % displays parameters
noFit = 0; % if 1, only estimation is done, but not actual fitting


% nargin
% if nargin < 4 || isempty(fitParameters) || isempty(parameters) || isempty(coordList) || isempty(intensities)
if nargin < 4 || isempty(fitParameters) || isempty(parameters) || isempty(intensities)
    error('not enough input arguments for GaussFitND!')
end

% check whether we have to create the coordList ourselves
if isempty(coordList)
    sizeInt = size(intensities);
    nd = length(sizeInt);
    % check for vector
    if nd == 2 && any(sizeInt == 1)
        % for a vector, coordList is just 1:n
        coordList = (1:length(intensities))';
    else
        % make coordList with eval to have n-dims
        inputString = sprintf('1:%i,',sizeInt);
        inputString = inputString(1:end-1);
        eval(['[',sprintf('X%i ',1:nd),']=ndgrid(' inputString ');']);
        eval(['coordList=[',sprintf('X%i(:) ',1:nd),'];']);
    end
end

% make sure image is a list
intensities = intensities(:);

% get dimensionality
nDims = size(coordList,2);
nCoords = size(coordList,1);

if nCoords ~= length(intensities)
    error('GaussFitND needs the coordinate of every pixel/voxel')
end

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
if nargin < 5 || isempty(isNormed);
    isNormed = def_isNormed;
end

if nargin < 6
    fitOptions = def_fitOptions;
else
    fitOptions = optimset(def_fitOptions,fitOptions);
end

% check nargout
if nargout > 1
    doStatistics = 1;
else
    doStatistics = 0;
end

%==========================


%==========================
% ESTIMATE PARAMETERS
%==========================
if verbose
    for i=1:nKernels
        disp(sprintf(['orig parms ',repmat('%1.4f ',1,nParameters)],parameters(i,:)))
    end
end

if any(isnan(parameters))


    % test whether we can estimate at all
    if nKernels > 1 || nDims > 3
        error('Parameter estimation has not been implemented for multiple kernels or more than 3 dimensions')
    end

    % make full image from coordList. Put NaN wherever there is nothing
    % transform to 1:n
    tmpCoordList = coordList;
    minList = min(tmpCoordList,[],1);
    tmpCoordList = tmpCoordList - repmat(minList,nCoords,1) + 1;
    % minimum image size is the maximum coordinates
    maxList = max(tmpCoordList,[],1);
    fullImage = repmat(NaN,maxList);
    % to assign values, coordinates must be indices, and sub2Ind doesn't accept
    % a coordinate matrix
    tmpCoordList = mat2cell(tmpCoordList,nCoords,ones(1,nDims));
    tmpIdxList = sub2ind(maxList,tmpCoordList{:});
    % assign intensities, have NaN otherwise
    fullImage(tmpIdxList) = intensities;

    % find parameters to estimate. Get NaN-indices so that we can
    % (theoretically) have only some of the centers or sigmas estimated
    centerIdx = find(isnan(parameters(1:nDims)));
    backgroundIdx = find(isnan(parameters(end)));
    amplitudeIdx = find(isnan(parameters(nDims+1)));
    sigmaIdx = find(isnan(parameters(nDims+2:end-1)));


    % do estimation only for up to 3D images, as I'm not going to write a
    % n-d centroid code right now

    % estimate center
    if any(centerIdx)
        % find centroid. Use exponent 10 to become a bit less sensitive to
        % noise. - Unfortunately, this means being more sensitive to
        % asymmetries, but that's why we do fitting
        center = centroid3D(fullImage,10);

        % transform center
        center = center(1:nDims) + minList - 1;

        % assign parameters
        parameters(centerIdx) = center(centerIdx);
    end
    % read new center. Transform to fullImage coordinates
    center = parameters(1:nDims) - minList + 1;



    % estimate background
    if any(backgroundIdx)
        % background is taken as the median of all the border pixels. If all is
        % NaN-masked, we might get only a few pixels to estimate the background
        % from, but at least they will be fairly definite background pixels.
        bgIdx = zeros(nCoords,nDims);
        % get border pixels by finding min and max coordinates along everz
        % dimension
        for i=1:nDims
            bgIdx(:,i) = tmpCoordList{i} == 1 | tmpCoordList{i} == maxList(i);
        end
        bgIdx = any(bgIdx,2);
        parameters(end) = nanmedian(fullImage(tmpIdxList(bgIdx)));
    end

    % read background and subtract
    background = parameters(end);
    fullImage = fullImage - background;

    if any(amplitudeIdx) && ~isNormed
        % amplitude is the maximum of a 3^nDims pixel region around the center.
        % Stamp3D is another of those functions that would have to be extended
        % for N-D.
        % if isNormed, there is no need for estimating amplitude
        subImage = stamp3d(fullImage,repmat(3,[1,nDims]),floor(center));
        parameters(nDims+1) = max(subImage(:));
    end
    % read amplitude, make image into 0/1
    amplitude = parameters(nDims+1);
    fullImage = fullImage./amplitude;

    if any(sigmaIdx)
        % for sigma: Threshold image at 1/sqrt(e) For a Gaussian, the width of
        % the remaining image should be two sigma
        fullImage = fullImage>exp(-0.5);
        % find largest group in threshold
        fullImage = bwlabeln(fullImage);
        [number,entry] = getMultiplicity(fullImage(:));
        [dummy,idx] = max(number(2:end));
        % get indices of largest group
        idx = find(fullImage == entry(idx+1));
        % get coordinates of largest group. Limit to 3D
        [cx,cy,cz] = ind2sub(maxList,idx);
        coords = [cx,cy,cz];
        minCoords = min(coords(:,1:nDims),[],1);
        maxCoords = max(coords(:,1:nDims),[],1);
        % sigma is half the extent of the thresholded region
        sigma = (maxCoords-minCoords+1)/2;
        % take care of Sxy
        if isSxy
            sigmaXY = mean(sigma(1:2));
            sigma = [sigmaXY,sigma(3:end)];
        end
        parameters(nDims+1+sigmaIdx) = sigma(sigmaIdx);
    end
    % no need to read sigma here


    if verbose
        for i=1:nKernels
            disp(sprintf(['est parms  ',repmat('%1.4f ',1,nParameters)],parameters(i,:)))
        end
    end


end % if estimate

%==========================

if noFit
    if nargout > 1
        [Q,chiSquared,degreesOfFreedom]=deal([]);
    end
    return
end

%==========================
% FIT
%==========================

% Choose fitting scheme:
% a) if sigmas are fixed, do all at once
% b) if not, do iterative fitting

anyAmplitude = any(strncmpi('a',fitParameters,1));
bgIdx = find(strncmpi('b',fitParameters,1));
anyBackground = any(bgIdx);
if anyBackground && any(bgIdx < nFitParameters)
    error('for fitting, background has to be at the last position in the fitParameters')
end
% In case we fit the amplitude: make amplitude comparable in magnitude to coordinates to
% avoid fitting problems. The amplitudes (max of image-background)
% are set to be roughly the average of the coordinates. Accordingly, the
% intensities and the background will have to be scaled, too!
if anyAmplitude || anyBackground
    % yes, coordList is correct here!
    %intensityScaling = mean(abs(coordList(:))) / max(intensities - parameters(1,end));

    % when things didn't work as with the MMF, I tried the very simple 0..1
    % norm from fitTest. Surprisingly, it works. In other words, the
    % norming of the amplitude is something that should be investigated
    % further - but I don't have the time right now.
    intensityScaling = 1/max(intensities);
    parameters(:,[nDims+1,end]) = ...
        parameters(:,[nDims+1,end]) * intensityScaling;
    intensities = intensities * intensityScaling;
end



% mask image
mask = zeros(size(intensities));
for k = 1:nKernels
    dist = normList(coordList-repmat(parameters(k,1:nDims),nCoords,1));
    % %     mask = mask | dist < 3*parameters(k,nDims+2);
    mask = mask | dist < 4*parameters(k,nDims+2); %--KJ (more reliable fit with bigger area)
end
intensities = intensities(mask);
coordList = coordList(mask,:);
nCoords = size(coordList,1);



%if fitSigma
% check for whether sigmas are to be fitted.
%fitSigma = any(strncmpi('s',fitParameters,1));
if 0
    %     % do iterative fitting
    %     % check for what we have to fit at all
    %     anyCenter = any(strncmpi('x',fitParameters,1));
    %
    %     for i=1:3
    %         % fit center if necessary
    %         if anyCenter
    %             % create xIdx
    %             xIdx = (1:nKernels*nDims)';
    %             % define fitFcn
    %             fitFcn = @(x)(GaussFitND_lsqnonlinFitFcn(x,fitParameters,parameters,xIdx,intensities,coordList,isNormed));
    %             % fit
    %             x = lsqnonlin(...
    %                 fitFcn,parameters(xIdx),[],[],fitOptions);
    %             % update parameters
    %             parameters(xIdx) = x;
    %             if verbose
    %                 for i=1:nKernels
    %                     disp(sprintf(['fitC parms ',repmat('%1.4f ',1,nParameters)],parameters(i,:)))
    %                 end
    %             end
    %         end
    %
    %         if anyAmplitude || anyBackground
    %             % create xIdx
    %             xIdx = zeros(nKernels*anyAmplitude + anyBackground,1);
    %             if anyAmplitude
    %                 xIdx(1:nKernels) = nKernels*nDims+1 : nKernels*(nDims+1);
    %             end
    %             if anyBackground
    %                 xIdx(end) = nParameters * nKernels;
    %             end
    %
    %             % define fitFcn
    %             fitFcn = @(x)(GaussFitND_lsqnonlinFitFcn(...
    %                 x,fitParameters,parameters,xIdx,intensities,coordList,isNormed));
    %             % fit
    %             x = lsqnonlin(...
    %                 fitFcn,parameters(xIdx),[],[],fitOptions);
    %             % update parameters
    %             parameters(xIdx) = x;
    %             if verbose
    %                 for i=1:nKernels
    %                     disp(sprintf(['fitA parms ',repmat('%1.4f ',1,nParameters)],parameters(i,:)))
    %                 end
    %             end
    %
    %         end
    %
    %         % fit sigma.
    %         % create xIdx
    %         xIdx = (nKernels*(nDims+1)+1:nKernels*(nDims*2+1-isSxy))';
    %         % define fitFcn
    %         fitFcn = @(x)(GaussFitND_lsqnonlinFitFcn(...
    %             x,fitParameters,parameters,xIdx,intensities,coordList,isNormed));
    %         % fit
    %         x = lsqnonlin(...
    %             fitFcn,parameters(xIdx),[],[],fitOptions);
    %         % update parameters
    %         parameters(xIdx) = x;
    %         if verbose
    %             for i=1:nKernels
    %                 disp(sprintf(['fitS parms ',repmat('%1.4f ',1,nParameters)],parameters(i,:)))
    %             end
    %         end
    %     end % loop iterative fitting
    %
    %
    %     % undo intensity scaling
    %     if anyAmplitude || anyBackground
    %         parameters(:,[nDims+1,end]) = ...
    %             parameters(:,[nDims+1,end]) / intensityScaling;
    %         intensities = intensities / intensityScaling;
    %     end
    %
    %     if doStatistics
    %
    %         % get jacobian, residuals for statistics. Find xIdx again
    %         xIdx = [];
    %         if anyCenter
    %             xIdx = [xIdx;(1:nKernels*nDims)'];
    %         end
    %         if anyAmplitude || anyBackground
    %             % create xIdx
    %             tmp = zeros(nKernels*anyAmplitude + anyBackground,1);
    %             if anyAmplitude
    %                 tmp(1:nKernels) = nKernels*nDims+1 : nKernels*(nDims+1);
    %             end
    %             if anyBackground
    %                 tmp(end) = nParameters * nKernels;
    %             end
    %             xIdx = [xIdx;tmp];
    %         end
    %         % we'll always be fitting sigma
    %         xIdx = [xIdx;(nKernels*(nDims+1)+1:nKernels*(nDims*2+1-isSxy))'];
    %
    %         % get residuals and jacobian
    %         [residuals, jacobian] = ...
    %             GaussFitND_lsqnonlinFitFcn(parameters(xIdx),...
    %             fitParameters,parameters,xIdx,intensities,coordList,isNormed);
    %     end
else
    % do one-step fitting

    % find the indices of the parameters that have to be fitted, and add to
    % xIdx
    xIdx = GaussFitND_xIdx(fitParameters,nFitParameters,nParameters,nKernels,nDims,isSxy);

    % define fitFcn
    fitFcn = @(x)(GaussFitND_lsqnonlinFitFcn(...
        x,fitParameters,parameters,xIdx,intensities,coordList,isNormed,nDims,isSxy));
    % fit
    x = lsqnonlin(...
        fitFcn,parameters(xIdx),[],[],fitOptions);

    % update parameters
    parameters(xIdx) = x;
    parameters(:,end) = parameters(end);

    % undo intensity scaling
    if anyAmplitude || anyBackground
        parameters(:,[nDims+1,end]) = ...
            parameters(:,[nDims+1,end]) / intensityScaling;
        intensities = intensities / intensityScaling;
    end

    if doStatistics
        % get residuals and jacobian
        [residuals, jacobian, gaussian] = ...
            GaussFitND_lsqnonlinFitFcn(parameters(xIdx),...
            fitParameters,parameters,xIdx,intensities,coordList,isNormed,nDims,isSxy);
    end



end

if verbose
    for i=1:nKernels
        disp(sprintf(['end  parms ',repmat('%1.4f ',1,nParameters)],parameters(i,:)))
    end
end

%=======================


%=======================
% Calculate statistics
%=======================

if doStatistics
    % calculate chi squared of fit, covarianceMatrix
    degreesOfFreedom = (nCoords-nFitParameters);
    chiSquared= sum(residuals.^2)/degreesOfFreedom;
    Q=inv(jacobian'*jacobian);
    sp = sqrt(chiSquared*diag(Q));
    % reshape sigmaParameters so that it looks like parameters
    sigmaParameters = repmat(NaN,size(parameters));
    sigmaParameters(xIdx) = sp;
    sigmaParameters(:,end) = sigmaParameters(end);

    if nargout > 5
        % make residual image

        tmpCoordList = coordList;
        minList = min(tmpCoordList,[],1);
        tmpCoordList = tmpCoordList - repmat(minList,nCoords,1) + 1;
        % minimum image size is the maximum coordinates
        maxList = max(tmpCoordList,[],1);
        residualImage = repmat(NaN,maxList);
        % to assign values, coordinates must be indices, and sub2Ind doesn't accept
        % a coordinate matrix
        tmpCoordList = mat2cell(tmpCoordList,nCoords,ones(1,nDims));
        tmpIdxList = sub2ind(maxList,tmpCoordList{:});
        % assign intensities, have NaN otherwise
        residualImage(tmpIdxList) = residuals;
        if nargout > 6
            resAndGauss = [residuals,gaussian];
        end
    end

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [residuals, jacobian, gaussian] = GaussFitND_lsqnonlinFitFcn(x,fitParameters,parameters,xIdx,intensities,coordList,isNormed,nDims,isSxy)
% GaussFitND_lsqnonlinFitFcn is the fit-routine in GaussFitND
% input: x: vector of unknowns.
%        fitParameters: 'labels' of the unknowns.
%        parameters: vector with all the data necessary to describe the
%                    Gaussian.
%        xIdx: index to place the unknowns in the parameter array
%        intensities : data to be fitted
%        coordList : coordinate corresponding to the intensities
%        isNormed : whether Gaussian is normed or not


% fill in parameters
parameters(xIdx) = x;
parameters(:,end) = parameters(end);

% check for amplitude norm - if normed, amplitude is a function of sigma.
% Fill it in.
if isNormed
    % amplitude is 1/(2*pi)^(nDims/2)*1/(sx*sy*sz)
    sigmas = parameters(:,end-(nDims-isSxy):end-1);
    if isSxy
        sigmas(:,1) = sigmas(:,1).^2;
    end
    parameters(:,nDims+1) = 1/(2*pi)^(nDims/2)*1./prod(sigmas,2);
end
    

% calculate gradient, Gaussian
[jacobian,gaussian] = ...
    GaussFitND_gradient(fitParameters,parameters,coordList,isNormed);

% get residuals
residuals = gaussian - intensities;

% disp parameters, suggested step
% disp(sprintf(['parms ',repmat('%1.4f ',1,length(parameters))],parameters))
% step = jacobian \ residuals;
% disp(sprintf(['step ',repmat('%1.4f ',1,length(step))],step))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function xIdx = GaussFitND_xIdx(fitParameters,nFitParameters,nParameters,nKernels,nDims,isSxy)
% GaussFitND_xIdx finds the incides of the parameters that have to be
% fitted.

% find the indices of the parameters that have to be fitted, and add to
% xIdx
% If we also fit the background, there will only be one index (the very
% last)
% xIdx should be in rows, i.e. if there are two kernels, there should first
% be all the indices to the first kernel, then all the indices to the
% second kernel etc.
fitBg = any(strncmpi('b',fitParameters,1));

% xIdx = zeros((nFitParameters - fitBg) * nKernels + fitBg,1);
% if fitBg
%     xIdx(end) = nParameters * nKernels;
% end
xIdx = zeros(nFitParameters - fitBg,nKernels);

% we don't need to go over background (we force it to be at the end of the
% list earlier in the program)
for iParm = 1:nFitParameters - fitBg
    currentParm = fitParameters{iParm};
    switch currentParm(1)
        case {'X','x'}

            % read dimension
            currentDim = str2double(currentParm(2:end));

            % add index to xIdx
            xIdx(iParm,:) = (currentDim-1) * nKernels + 1 : currentDim * nKernels;



        case {'S','s'}

            % take care of the special case sxy
            if strcmpi(currentParm(2:end),'xy')
                % use sx. Also correct for the fact that sigmas start
                % after center and amplitude
                currentDim = 1  + nDims + 1;
            else
                % get dimension. If there is a sxy, we lack S2. Also
                % correct for the fact that sigmas start after center
                % and amplitude
                currentDim = str2double(currentParm(2:end)) - isSxy;
                currentDim = currentDim  + nDims + 1;
            end

            % add index to xIdx
            xIdx(iParm,:) =...
                (currentDim-1) * nKernels + 1 : currentDim * nKernels;

        case {'A','a'}
            % add index to xIdx
            xIdx(iParm,:) = ...
                nDims * nKernels + 1 : (nDims+1) * nKernels;

        case {'B','b'}
            % we've taken care of that already

        otherwise
            error('fitParameter ''%s'' not recognized!',currentParm)
    end % switch fitParm
end % loop fitParm

% make xIdx into a vector, add background
if fitBg
    xIdx = [xIdx(:);nKernels*nParameters];
else
    xIdx = xIdx(:);
end
