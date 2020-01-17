function [kalmanFilterInfoOut,errFlag] = plusTipKalmanGainLinearMotion(trackedFeatureIndx,...
    frameInfo,kalmanFilterInfoTmp,propagationScheme,kalmanFilterInfoIn,probDim,...
    filterInfoPrev,costMatParam,initFunctionName)
%KALMANGAINLINEARMOTION uses the Kalman gain from linking to get better estimates of the state vector, its covariance matrix, state noise and its variance
%
%SYNOPSIS [kalmanFilterInfoOut,errFlag] = plusTipKalmainGainLinearMotion(trackedFeatureIndx,...
%    frameInfo,kalmanFilterInfoTmp,propagationScheme,kalmanFilterInfoIn,probDim,...
%    filterInfoPrev,costMatParam,initFunctionName)
%
%INPUT  trackedFeatureIndx : Matrix showing connectivity between features
%                            in current frame (listed in last column of matrix) 
%                            and features in previous frames. A zero in 
%                            columns before last indicates that feature is
%                            not connected to any previous features.
%       frameInfo          : Structure with fields (for current frame):
%             .allCoord        : x,dx,y,dy,[z,dz] of features collected in one
%                                matrix.kalmanFunctions.
%             .amp             : Amplitudes of PSFs fitting detected features. 
%                                1st column for values and 2nd column 
%                                for standard deviations.
%       kalmanFilterInfoTmp: Structure with fields (for current frame):
%             .stateVec        : Kalman filter prediction of state
%                                vector of all features in current frame 
%                                based on all 3 motion models.
%             .stateCov        : Kalman filter prediction of state
%                                covariance matrix of all features in
%                                current frame based on all 3 motion models.
%             .obsVec          : Kalman filter prediction of the
%                                observed variables for all features in 
%                                current frame based on all 3 motion models.
%       propagationScheme  : Matrix indicating the propagation scheme that
%                            yielded the lowest cost for a link between two
%                            features.
%       kalmanFilterInfoIn : Structure with fields (for all previous frames):
%             .stateVec        : Kalman filter state vector for all features.
%             .stateCov        : Kalman filter state covariance matrix
%                                for all features.
%             .noiseVar        : Variance of state noise for all feature.
%             .stateNoise      : Estimated state noise for each feature in
%                                frame.
%             .scheme          : 1st column: propagation scheme connecting
%                                feature to previous feature. 2nd column:
%                                propagation scheme connecting feature to
%                                next feature.
%       probDim            : Problem dimension. 2 (for 2D) or 3 (for 3D).
%       filterInfoPrev     : Structure with fields (for current frame):
%             .stateVec        : Kalman filter state vector for each
%                                feature in frame.
%             .stateCov        : Kalman filter state covariance matrix
%                                for each feature in frame.
%             .noiseVar        : Variance of state noise for each
%                                feature in frame.
%                            Enter [] if there is no previous information.
%       costMatParam       : Linking cost matrix parameters.
%       initFunctionName   : Name of function for Kalman filter
%                            initialization.
%
%OUTPUT kalmanFilterInfoOut: Structure with fields (for all frames upto current):
%             .stateVec        : Kalman filter state vector for all features.
%             .stateCov        : Kalman filter state covariance matrix
%                                for all features.
%             .noiseVar        : Variance of state noise for all features.
%             .stateNoise      : Estimated state noise for each feature in
%                                frame.
%             .scheme          : 1st column: propagation scheme connecting
%                                feature to previous feature. 2nd column:
%                                propagation scheme connecting feature to
%                                next feature.
%       errFlag            : 0 if function executes normally, 1 otherwise.
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

kalmanFilterInfoOut = [];
errFlag = [];

%% Input

%check whether correct number of input arguments was used
if nargin < 9
    disp('--plusTipKalmainGainLinearMotion: Incorrect number of input arguments!');
    errFlag  = 1;
    return
end

if isempty(filterInfoPrev)
    usePriorInfo = 0;
else
    usePriorInfo = 1;
end

%% Gain calculation and update

%take absolute value of all noise variances - we do this because we make
%noiseVar negative in the initialization for first appearances...
for iFrame = 1 : length(kalmanFilterInfoIn)
    kalmanFilterInfoIn(iFrame).noiseVar = abs(kalmanFilterInfoIn(iFrame).noiseVar);
end

%copy kalmanFilterInfoIn into kalmanFilterInfoOut
kalmanFilterInfoOut = kalmanFilterInfoIn;

%get number of features in current frame and frame number
[numFeatures,iFrame] = size(trackedFeatureIndx);

%construct Kalman filter observation matrix
observationMat = [diag(ones(probDim,1)) zeros(probDim)];

%go over all features in current frame
for iFeature = 1 : numFeatures

    %find index of feature in previous frame that this feature is connected to
    iFeaturePrev = trackedFeatureIndx(iFeature,end-1);

    %if this feature is connected to a feature in previous frame
    if iFeaturePrev ~= 0

        %find propagation scheme leading to this link and save in kalmanFilterInfo
        iScheme = propagationScheme(iFeaturePrev,iFeature);
        kalmanFilterInfoOut(iFrame).scheme(iFeature,1) = iScheme; %to current feature
        kalmanFilterInfoOut(iFrame-1).scheme(iFeaturePrev,2) = iScheme; %from previous feature

        %get the corresponding state vector, covariance and observation vector
        stateVecOld = kalmanFilterInfoTmp.stateVec(iFeaturePrev,:,iScheme)';
        stateCovOld = kalmanFilterInfoTmp.stateCov(:,:,iFeaturePrev,iScheme);
        obsVecOld = kalmanFilterInfoTmp.obsVec(iFeaturePrev,:,iScheme)';
        
        %calculate Kalman gain
        kalmanGain = stateCovOld*observationMat'*...
            inv(observationMat*stateCovOld*observationMat'+...
            diag(frameInfo.allCoord(iFeature,2:2:end).^2));

        %estimate state noise in previous frame and save in kalmanFilterInfo
        stateNoise = kalmanGain * (frameInfo.allCoord(iFeature,1:2:end)' - obsVecOld);
        kalmanFilterInfoOut(iFrame-1).stateNoise(iFeaturePrev,:) = stateNoise';
        
        %update estimate of state vector in current frame
        stateVec = stateVecOld + stateNoise;
        
        %update estimate of state covariance matrix in current frame
        stateCov = stateCovOld - kalmanGain*observationMat*stateCovOld;

        %find indices of previous features connected to this feature
        indx = trackedFeatureIndx(iFeature,1:end-1); 
        
        %determine where the track starts
        indxLength = length(find(indx));
        
        %collect all of the error terms
        stateNoiseAll = zeros(indxLength,2*probDim);
        j = 0;
        for i = iFrame-indxLength : iFrame-1
            j = j + 1;
            stateNoiseAll(j,:) = kalmanFilterInfoOut(i).stateNoise(indx(i),:);
        end

        %impose isotropy (all directions are equivalent)
        stateNoisePos = stateNoiseAll(:,1:probDim);
        stateNoiseVel = stateNoiseAll(:,probDim+1:2*probDim);

        %estimate positional and speed noise variances in current frame
        noiseVar = zeros(1,2*probDim);
        noiseVar(1:probDim) = var(stateNoisePos(:));
        noiseVar(probDim+1:2*probDim) = var(stateNoiseVel(:));

        %save this information in kalmanFilterInfo
        kalmanFilterInfoOut(iFrame).stateVec(iFeature,:) = stateVec';
        kalmanFilterInfoOut(iFrame).stateCov(:,:,iFeature) = stateCov;
        kalmanFilterInfoOut(iFrame).noiseVar(:,:,iFeature) = diag(noiseVar);

    else %if this feature is not connected to anything in previous frame

        %initialize Kalman filter for this feature
        if usePriorInfo %use prior information if supplied
            kalmanFilterInfoOut(iFrame).stateVec(iFeature,:) = filterInfoPrev.stateVec(iFeature,:);
            kalmanFilterInfoOut(iFrame).stateCov(:,:,iFeature) = filterInfoPrev.stateCov(:,:,iFeature);
            kalmanFilterInfoOut(iFrame).noiseVar(:,:,iFeature) = filterInfoPrev.noiseVar(:,:,iFeature);
        else
            featureInfo.allCoord = frameInfo.allCoord(iFeature,:);
            featureInfo.num = 1;
            eval(['[filterTmp,errFlag] = ' initFunctionName '(featureInfo,probDim,'...
                'costMatParam);']);
            kalmanFilterInfoOut(iFrame).stateVec(iFeature,:) = filterTmp.stateVec;
            kalmanFilterInfoOut(iFrame).stateCov(:,:,iFeature) = filterTmp.stateCov;
            kalmanFilterInfoOut(iFrame).noiseVar(:,:,iFeature) = filterTmp.noiseVar;
        end

    end %(if iFeaturePrev ~= 0 ... else ...)

end %(for iFeature = 1 : numFeatures)


%% %%%% ~~ the end ~~ %%%%
