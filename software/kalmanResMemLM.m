function kalmanFilterInfo = kalmanResMemLM(numFrames,numFeatures,probDim)
%KALMANRESMEMLM reserves memory for Kalman filter structure for linear motion model
%
%SYNPOSIS kalmanFilterInfo = kalmanResMemLM(numFrames,numFeatures,probDim)
%
%INPUT   numFrames  : Number of frames in movie.
%        numFeatures: An array with number of feaures in each frame.
%        probDim    : Problem dimensionality.
%
%OUTPUT   kalmanFilterInfo: Structure array with number of entries equal to 
%                           number of frames in movie. Contains the fields:
%             .stateVec        : Kalman filter state vector for each
%                                feature in frame.
%             .stateCov        : Kalman filter state covariance matrix
%                                for each feature in frame.
%             .noiseVar        : Variance of state noise for each
%                                feature in frame.
%             .stateNoise      : Estimated state noise for each feature in
%                                frame.
%             .scheme          : 1st column: propagation scheme connecting
%                                feature to previous feature. 2nd column:
%                                propagation scheme connecting feature to
%                                next feature.
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

%calculate vector size
vecSize = 2 * probDim;

%go over all frames and reserve memory
for iFrame = numFrames : -1 : 1

    kalmanFilterInfo(iFrame) = struct('stateVec',zeros(numFeatures(iFrame),vecSize),...
        'stateCov',zeros(vecSize,vecSize,numFeatures(iFrame)),...
        'noiseVar',zeros(vecSize,vecSize,numFeatures(iFrame)),...
        'stateNoise',zeros(numFeatures(iFrame),vecSize),...
        'scheme',zeros(numFeatures(iFrame),2));

end
