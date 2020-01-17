function [kalmanFilterInfo,errFlag] = kalmanInitLinearMotion(frameInfo,...
    probDim,costMatParam)
%KALMANINITLINEARMOTION initializes Kalman filter state vector and covariance matrix for features in a frame
%
%SYNOPSIS [kalmanFilterInfo,errFlag] = kalmanInitLinearMotion(frameInfo,...
%    probDim,costMatParam)
%
%INPUT  frameInfo       : Structure with fields (for 1 frame):
%             .allCoord     : Image coordinate system x,dx,y,dy,[z,dz] of 
%                             detected features and their uncertainties (in
%                             pixels).
%             .num          : Number of features in frame.
%       probDim         : Problem dimension. 2 (for 2D) or 3 (for 3D).
%       costMatParam    : Linking cost matrix parameters. In addition
%                         to the fields needed for costMatLinearMotionLink,
%                         if it has the field 
%             .kalmanInitParam, then sub-fields within will be used for
%                               initialization. 
%                   .convergePoint: Convergence point (x, y, [z]
%                                   coordinates) of tracks if motion is
%                                   radial, in image coordinate system.
%                                   Should be a row vector. If supplied,
%                                   radial form is assumed and value is
%                                   used to estimate initial velocities. If
%                                   not supplied, then initial velocity is
%                                   taken as zero.
%                   .searchRadiusFirstIteration: max number of pixels to
%                                   use in the first frame pair where the
%                                   search is essentially a nearest
%                                   neighbor. here we use this parameter to
%                                   get the initial noise variance. in
%                                   kathryn's alternate function
%                                   costMatLinearMotionLink_EB3.m, the
%                                   min/max search radii are not applied in
%                                   the first iteration so that this
%                                   initial radius can be sufficiently big
%                                   if features are moving quickly (ie
%                                   faster than the max search radius)
%                   .initVelocity : Initial velocity guess (vx, vy, [vz]),
%                                   in whatever units coordinates are in
%                                   per frame.
%
%OUTPUT kalmanFilterInfo: Structure with fields:
%             .stateVec     : State vector for each feature.
%             .stateCov     : State covariance matrix for each feature.
%             .noiseVar     : Variance of state noise for each feature.
%       errFlag         : 0 if function executes normally, 1 otherwise.
%
%
%Khuloud Jaqaman, March 2007
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

%% Output

kalmanFilterInfo = [];
errFlag = 0;

%% Input

%check whether correct number of input arguments was used
if nargin < 3
    disp('--kalmanInitLinearMotion: Incorrect number of input arguments!');
    errFlag  = 1;
    return
end

%extract min. and max. search radii and brownStdMult from costMatParam
minSearchRadius = costMatParam.minSearchRadius;
maxSearchRadius = costMatParam.maxSearchRadius;
brownStdMult = costMatParam.brownStdMult;

%check whether additional initialization parameters were input
if isfield(costMatParam,'kalmanInitParam')
    
    initParam = costMatParam.kalmanInitParam;
    
    if isfield(initParam,'convergePoint')
        convergePoint = initParam.convergePoint;
    else
        initParam.convergePoint = [];
        convergePoint = [];
    end
    
    if isfield(initParam,'initVelocity')
        initVelGuess = initParam.initVelocity;
    else
        initParam.initVelocity = [];
        initVelGuess = [];
    end
    
    % calculate noiseVarInit from max search radius (if given) or from
    % min/max search radii
    if isfield(initParam,'searchRadiusFirstIteration') && ~isempty(initParam.searchRadiusFirstIteration)
        searchRadiusFirstIteration = initParam.searchRadiusFirstIteration;
        noiseVarInit = -(searchRadiusFirstIteration/brownStdMult)^2 / probDim; %negative to flag as first appearance
    else
        noiseVarInit = (maxSearchRadius / brownStdMult) ^ 2 / probDim;
        %         noiseVarInit = (mean([minSearchRadius maxSearchRadius]) / brownStdMult) ^ 2 / probDim;
    end

else
    initParam.convergePoint = [];
    convergePoint = [];
    initVelGuess = [];
    noiseVarInit = (mean([minSearchRadius maxSearchRadius]) / brownStdMult) ^ 2 / probDim;
end

if ~isempty(convergePoint)
    [nRow,nCol] = size(convergePoint);
    if nRow ~= 1 && nCol == 1
        convergePoint = convergePoint';
    end
    nCol = size(convergePoint,2);
    if nCol ~= probDim
        disp('--kalmanInitLinearMotion: initParam.convergePoint of wrong dimension!');
        errFlag = 1;
    end
end
    
if ~isempty(initVelGuess)
    [nRow,nCol] = size(initVelGuess);
    if nRow ~= 1 && nCol == 1
        initVelGuess = initVelGuess';
    end
    nCol = size(initVelGuess,2);
    if nCol ~= probDim
        disp('--kalmanInitLinearMotion: initParam.initVelocity of wrong dimension!');
        errFlag = 1;
    end
end
    
%exit if there are problems with input
if errFlag
    disp('--kalmanInitLinearMotion: Please fix input parameters.');
    return
end

%% Initialization

%find number of features in frame
numFeatures = frameInfo.num;

%estimate initial velocity
if isempty(initVelGuess) %if there is no initial guess for velocity
    
    if isempty(convergePoint) %if there is no convergence point information
        
        %assume zero initial velocity
        velocityInit = zeros(numFeatures,probDim);
        
    else %if there is a convergence point
        
        %assign initial speed
        speedInit = 1;
        
        %get displacement and distance of features from convergence point
        displacement = repmat(convergePoint,numFeatures,1) - ...
            frameInfo.allCoord(1:2:end);
        distance = sqrt(sum(displacement.^2,2));
        
        %calculate initial velocity
        velocityInit = speedInit * displacement ./ repmat(distance,1,probDim);
        
    end
    
else %if there is an initial guess for velocity
    
    velocityInit = repmat(initVelGuess,numFeatures,1);
    
end

%initialize state vector
kalmanFilterInfo.stateVec = [frameInfo.allCoord(:,1:2:end) velocityInit];

%initialize state covariance matrix
for iFeature = 1 : numFeatures
    posVar = frameInfo.allCoord(iFeature,2:2:end).^2;
    posVar = max(eps,posVar);
    kalmanFilterInfo.stateCov(:,:,iFeature) = diag([posVar 4*ones(1,probDim)]);
end

%initialize state noise variance
for iFeature = 1 : numFeatures
    kalmanFilterInfo.noiseVar(:,:,iFeature) = diag(...
        ones(2*probDim,1)*noiseVarInit);
end


%% %%%% ~~ the end ~~ %%%%
