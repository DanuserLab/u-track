function kalmanFilterInfo = kalmanReverseLinearMotion(kalmanFilterInfo,probDim)
%KALMANREVERSELINEARMOTION revese Kalman filter information in time
%
%SYNPOSIS kalmanFilterInfo = kalmanResMemLM(numFrames,numFeatures,probDim)
%
%INPUT  kalmanFilterInfo: Kalman filter information from a previous round
%                         of linking (as initialized by kalmanResMemLM).
%
%OUTPUT kalmanFilterInfo: Kalman filter information reversed in time.
%
%Khuloud Jaqaman, September 2008
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

%reverse time
kalmanFilterInfo = kalmanFilterInfo(end:-1:1);


%go over all frames and reverse velocity
for iFrame = length(kalmanFilterInfo) : -1 : 1
    kalmanFilterInfo(iFrame).stateVec(:,probDim+1:end) = ...
        -kalmanFilterInfo(iFrame).stateVec(:,probDim+1:end);
end
