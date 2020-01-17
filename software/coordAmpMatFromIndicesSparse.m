function trackedFeatureInfo = coordAmpMatFromIndicesSparse(trackedFeatureIndx,...
    movieInfo,numFrames,probDim)
%COORDAMPMATFROMINDICESSPARSE converts a matrix of particle index connectivity to a matrix of coordinate and amplitude connectivity
%
%SYNOPSIS trackedFeatureInfo = coordAmpMatFromIndicesSparse(trackedFeatureIndx,...
%    movieInfo,numFrames,probDim)
%
%INPUT  
%
%OUTPUT
%
%REMARKS This is a quick hack, I took these lines out of the function
%"linkFeaturesKalmanSparse".
%
%Khuloud Jaqaman, June 2010
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

%% Conversion

%reserve space for matrix
trackedFeatureInfo = sparse(size(trackedFeatureIndx,1),8*numFrames);
% trackedFeatureInfo = NaN(size(trackedFeatureIndx,1),8*numFrames);

%for now, the matrix always has space for z, hence I have to treat 2D and
%3D differently.
%at some point I can change the postprocessing functions to handle 2D and
%3D differently, and so I can save them differently, in which case I won't
%need the if-else here.
if probDim ==2

    %go over all frames
    for iFrame = 1 : numFrames

        %find rows that have a feature index
        indx1 = find(trackedFeatureIndx(:,iFrame)~=0);

        %if there are detected features in this frame
        if ~isempty(indx1)

            %find the feature index
            indx2 = trackedFeatureIndx(indx1,iFrame);

            %store its information
            trackedFeatureInfo(indx1,8*(iFrame-1)+1:8*iFrame) = ...
                [movieInfo(iFrame).allCoord(indx2,1:2:2*probDim) ...
                zeros(length(indx2),1) movieInfo(iFrame).amp(indx2,1) ...
                movieInfo(iFrame).allCoord(indx2,2:2:2*probDim) ...
                zeros(length(indx2),1) movieInfo(iFrame).amp(indx2,2)];

        end

    end

else

    %go over all frames
    for iFrame = 1 : numFrames

        %find rows that have a feature index
        indx1 = find(trackedFeatureIndx(:,iFrame)~=0);

        %if there are detected features in this frame
        if ~isempty(indx1)

            %find the feature index
            indx2 = trackedFeatureIndx(indx1,iFrame);

            %store its information
            trackedFeatureInfo(indx1,8*(iFrame-1)+1:8*iFrame) = ...
                [movieInfo(iFrame).allCoord(indx2,1:2:2*probDim) ...
                movieInfo(iFrame).amp(indx2,1) ...
                movieInfo(iFrame).allCoord(indx2,2:2:2*probDim) ...
                movieInfo(iFrame).amp(indx2,2)];

        end

    end

end

%% %%%%% ~~ the end ~~ %%%%%
