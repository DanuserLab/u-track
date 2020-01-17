function [trackType,xyzVelS,xyzVelE,noiseStd,trackCenter,trackMeanDispMag,...
    errFlag] = estimTrackTypeParamRDS(trackedFeatIndx,trackedFeatInfo,...
    kalmanFilterInfo,lenForClassify,probDim)
%estimTrackTypeParamRDS ...
%
%SYNOPSIS [trackType,xyzVelS,xyzVelE,noiseStd,trackCenter,trackMeanDispMag,...
%    errFlag] = estimTrackTypeParamRDS(trackedFeatIndx,trackedFeatInfo,...
%    kalmanFilterInfo,lenForClassify,probDim)
%
%INPUT  trackedFeatIndx : Connectivity matrix of features between time
%                         points from initial linking. Rows indicate tracks, while columns
%                         indicate frames. A track that ends before the
%                         last frame is followed by zeros, and a track
%                         that starts at a time after the first frame
%                         is preceded by zeros. 
%       trackedFeatInfo : The positions and amplitudes of the tracked
%                         features from initial linking. Number of rows = number of tracks,
%                         while number of columns = 8*number of time points.
%                         Each row consists of
%                         [x1 y1 z1 a1 dx1 dy1 dz1 da1 x2 y2 z2 a2 dx2 dy2 dz2 da2 ...]
%                         in image coordinate system (coordinates in
%                         pixels). NaN is used to indicate time points
%                         where the track does not exist.
%       kalmanFilterInfo: Kalman filter information as calculated in
%                         linkFeaturesKalman for the linear motion model.
%       lenForClassify  : Minimum length of a track to classify it as
%                         directed or Brownian.
%       probDim        : Problem dimensionality. 2 (for 2D) or 3 (for 3D).
%
%OUTPUT trackType
%       xyzVelS
%       xyzVelE
%       noiseStd
%       trackCenter
%       trackMeanDispMag
%       errFlag         : 0 if function executes normally, 1 otherwise.
%
%Khuloud Jaqaman, April 2007
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

trackType = [];
xyzVelS = [];
xyzVelE = [];
noiseStd = [];
trackCenter = [];
errFlag = 0;

%% Input

%check whether correct number of input arguments was used
if nargin ~= nargin('estimTrackTypeParamRDS')
    disp('--estimTrackTypeParamRDS: Incorrect number of input arguments!');
    return
end

%get number of tracks from initial linking and number of frames
[numTracksLink,numFrames] = size(trackedFeatIndx);

%% Estimation of type and motion parameters

%reserve memory for output variables
trackType = NaN(numTracksLink,1);
xyzVelS = zeros(numTracksLink,probDim);
xyzVelE = zeros(numTracksLink,probDim);
noiseStd = zeros(numTracksLink,1);
trackCenter = zeros(numTracksLink,probDim);
trackMeanDispMag = NaN(numTracksLink,1);

%get the start times, end times and lifetimes of all tracks
trackSEL = getTrackSEL(trackedFeatInfo);
trackStartTime = trackSEL(:,1);
trackEndTime   = trackSEL(:,2);
trackLifeTime  = trackSEL(:,3);
clear trackSEL

%assign the asymmetry parameter thresholds that indicate directed motion
%for different track lengths
if probDim == 2
    % % %     %50th percentile:
    % % %     asymThresh = [[NaN NaN 1.6 0.95 0.75 0.66 0.64 0.6 0.57 0.55 0.54 ...
    % % %         0.53 0.52 0.5 0.5]'; 0.5*ones(numFrames-15,1)];
    % % %     %80th percentile:
    % % %     asymThresh = [[NaN NaN 3.5 2 1.5 1.4 1.3 1.2 1.2 1.1 1.1 1.1 1.1 1.1 1.1 ...
    % % %         1.05 1.05 1.05 1.05 1.05]'; ones(numFrames-20,1)];
    %90th percentile:
    asymThresh = [[NaN NaN 5 2.7 2.1 1.8 1.7 1.6 1.5 1.45 1.45 1.4 1.4 ...
        1.4 1.4 1.4 1.4 1.35 1.35 1.35]'; 1.3*ones(numFrames-20,1)];
    % % %     %99th percentile:
    % % %     asymThresh = [[NaN NaN 9 5 3.5 3 2.7 2.5 2.4 2.4 2.3 2.2 ...
    % % %         2.2 2.2 2.2 2.2 2.2 2.2 2.2 2.2]'; 2.1*ones(numFrames-20,1)];
else
    %90th percentile:
    asymThresh = [[NaN NaN 2.9 1.9 1.5 1.4 1.3 1.3 1.2 1.2 1.2 1.2 ...
        1.2 1.2 1.2]'; 1.1*ones(numFrames-15,1)];
    % % %     %99th percentile:
    % % %     asymThresh = [[NaN NaN 5.5 3 2.5 2.3 2 2 1.9 1.9 1.9 1.8 1.8 ...
    % % %         1.8 1.7]'; 1.7*ones(numFrames-15,1)];
end

%go over all tracks
for iTrack = 1 : numTracksLink
    
    %get current track's coordinates
    currentTrack = (reshape(trackedFeatInfo(iTrack,:)',8,[]))';
    currentTrack = currentTrack(:,1:probDim);
    currentTrack = full(currentTrack(trackStartTime(iTrack):trackEndTime(iTrack),:));

    %calculate the track's center of mass
    trackCenter(iTrack,:) = mean(currentTrack);
    
    %calculate the track's mean displacement
    if size(currentTrack,1) > 1
        trackMeanDispMag(iTrack) = mean(sqrt(sum(diff(currentTrack,1,1).^2,2)));
    end
    
    %assign default value of track type (NaN implies track is too short to
    %determine its type)
    overallType = NaN;
    
    %if track's lifetime is at least "lenForClassify" frames
    %(shorter tracks are not reliable) ...
    if trackLifeTime(iTrack) >= lenForClassify

        %evaluate the asymmetry parameter as defined in the Huet et al. BJ 2006 paper
        asymmetry = asymDeterm2D3D(currentTrack);

        %if the asymmetry is larger than threshold ...
        overallType = asymmetry > asymThresh(trackLifeTime(iTrack));
        
    end %(if trackLifeTime(iTrack) >= lenForClassify)
    
    %store track type in vector
    trackType(iTrack) = overallType;
    
    %assign motion parameters for track based on overallType
    switch overallType

        case 1 %if track is directed
            
            %assign velocity
            xyzVelS(iTrack,:) = kalmanFilterInfo(trackStartTime(...
                iTrack)).stateVec(trackedFeatIndx(iTrack,...
                trackStartTime(iTrack)),probDim+1:2*probDim);
            xyzVelE(iTrack,:) = kalmanFilterInfo(trackEndTime(...
                iTrack)).stateVec(trackedFeatIndx(iTrack,...
                trackEndTime(iTrack)),probDim+1:2*probDim);
            
            %assign noise std
            noiseStd(iTrack) = sqrt( abs( kalmanFilterInfo(trackEndTime(...
                iTrack)).noiseVar(1,1,trackedFeatIndx(iTrack,...
                trackEndTime(iTrack))) ) );
            
        case 0 %if track is Brownian

            %give track a velocity of zero - which it
            %already has from initialization

            %assign noise std
            noiseStd(iTrack) = sqrt( abs( kalmanFilterInfo(trackEndTime(...
                iTrack)).noiseVar(1,1,trackedFeatIndx(iTrack,...
                trackEndTime(iTrack))) ) );

        otherwise %if track is undetermined
            
            %give track a velocity of zero - which it
            %already has from initialization

            %give track a noise std of 1 (this value will be
            %overwritten in the actual calculation of average displacement
            %and search radius
            noiseStd(iTrack) = 1;

    end %(switch overallType)

end %(for iTrack = 1 : numTracksLink)


%% ~~~ the end ~~~