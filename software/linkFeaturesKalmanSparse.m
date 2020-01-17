function [trackedFeatureIndx,trackedFeatureInfo,kalmanFilterInfo,...
    nnDistFeatures,prevCost,errFlag] = linkFeaturesKalmanSparse(movieInfo,...
    costMatFunc,costMatParam,kalmanFunctions,probDim,filterInfoPrev,...
    prevCost,verbose)
%LINKFEATURESKALMAN links features between consecutive frames using LAP and possibly motion propagation using the Kalman filter
%
%SYNOPSIS [trackedFeatureIndx,trackedFeatureInfo,kalmanFilterInfo,...
%    nnDistFeatures,prevCost,errFlag] = linkFeaturesKalmanSparse(movieInfo,...
%    costMatFunc,costMatParam,kalmanFunctions,probDim,filterInfoPrev,...
%    prevCost,verbose)
%
%INPUT  movieInfo      : Array of size equal to the number of frames
%                        in a movie, containing the fields:
%             .xCoord      : x-coordinates of detected features. 
%                            1st column: value, 2nd column: standard
%                            deviation (zeros if not available).
%             .yCoord      : y-coordinates of detected features.
%                            1st column: value, 2nd column: standard
%                            deviation (zeros if not available).
%             .zCoord      : z-coordinates of detected features.
%                            1st column: value, 2nd column: standard
%                            deviation (zeros if not available).
%                            Optional. Skipped if problem is 2D. Default: zeros.
%             .amp         : "Intensities" of detected features.
%                            1st column: values (ones if not available),
%                            2nd column: standard deviation (zeros if not
%                            available).
%             .num         : Number of features. 
%                            Optional. Calculated if not supplied.  
%             .allCoord    : x,dx,y,dy,[z,dz] of features collected in one
%                            matrix. Optional. Calculated if not supplied.
%             .nnDist      : Distance from each feature to its nearest
%                            neighbor. Optional. Calculated if not supplied.
%    %%%costMatName    : Name of cost matrix function used for linking. % Replaced!!!!!
%       costMatFunc    : Function handle pointing to the cost matrix function used for linking. % Updated by Carmen Klein Herenbrink and Brian Devree
%       costMatParam   : Parameters needed for cost matrix calculation. 
%                        Structure with fields specified by particular
%                        cost matrix used (costMatName).
%       kalmanFunctions: Names of Kalman filter functions for self-adaptive
%                        tracking. Structure with fields:
%             .reserveMem   : Reserves memory for kalmanFilterInfo.
%             .initialize   : Initializes the Kalman filter for an appearing
%                             feature.
%             .calcGain     : Calculates the Kalman gain after linking.
%                        For non-self-adaptive tracking, enter [].
%                        Optional. Default: [].
%       probDim        : Problem dimensionality. 2 (for 2D) or 3 (for 3D).
%                        Optional. If not input, dimensionality will be
%                        derived from movieInfo.
%       filterInfoPrev : Structure array with number of entries equal to
%                        number of frames in movie. Contains at least the
%                        fields:
%             .stateVec    : Kalman filter state vector for each feature in frame.
%             .stateCov    : Kalman filter state covariance matrix for each feature in frame.
%             .noiseVar    : Variance of state noise for each feature in frame.
%                            Optional. Enter [] or nothing if not to be used.
%       prevCost       : Matrix of costs of actual assignments in previous
%                        round of linking. Optional. Default: Empty.
%       verbose        : 1 to show calculation progress, 0 otherwise.
%                        Optional. Default: 1.
%
%OUTPUT trackedFeatureIndx: Connectivity matrix of features between frames.
%                           Rows indicate continuous tracks, while columns 
%                           indicate frames. A track that ends before the
%                           last frame is followed by zeros, and a track
%                           that starts in a frame after the first frame
%                           is preceded by zeros. 
%       trackedFeatureInfo: The positions and "intensities" of the tracked
%                           features, in the same units as input. 
%                           Number of rows = number of tracks.
%                           Number of columns = 8*number of frames.
%                           Each row consists of 
%                           [x1 y1 z1 a1 dx1 dy1 dz1 da1 x2 y2 z2 a2 dx2 dy2 dz2 da2 ...]
%                           NaN is used to indicate frames where the tracks
%                           do not exist.
%       kalmanFilterInfo  : Structure array with number of entries equal to 
%                           number of frames in movie. Contains the fields
%                           defined in kalmanFunctions.reserveMem (at
%                           least stateVec, stateCov and noiseVar).
%       nnDistFeatures    : Matrix indicating the nearest neighbor
%                           distances of features linked together within
%                           tracks.
%       prevCost          : Matrix of costs of actual assignments.
%       errFlag           : 0 if function executes normally, 1 otherwise.
%
%REMARKS Algorithm can handle cases where some frames do not have any
%        features at all. However, the very first frame must not be empty.
%
%Khuloud Jaqaman, March 2007
%
% Updated in Jan 2020 to incorporate the changes made by Carmen Klein Herenbrink 
% and Brian Devree from Copenhagen University to reduce the tracking time.
% Changes made in this function are modified the code to pass the function named 
% "costMatRandomDirectedSwitchingMotionLink" as a handle, instead of calling it through 
% an eval statement.
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

trackedFeatureIndx = [];
trackedFeatureInfo = [];
kalmanFilterInfo = [];
nnDistFeatures = [];
errFlag = [];

%% Input

%check whether correct number of input arguments was used
if nargin < 3
    disp('--linkFeaturesKalmanSparse: Incorrect number of input arguments!');
    errFlag  = 1;
    return
end

%check whether tracking is self-adaptive 
if nargin < 4 || isempty(kalmanFunctions)
    kalmanFunctions = [];
    selfAdaptive = 0;
else
    selfAdaptive = 1;
end

%check whether z-coordinates were input, making problem potentially 3D
if isfield(movieInfo,'zCoord')
    probDimT = 3;
else    
    probDimT = 2;
end

%assign problem dimensionality if not input
if nargin < 5 || isempty(probDim)
    probDim = probDimT;
else
    if probDim == 3 && probDimT == 2
        disp('--linkFeaturesKalmanSparse: Inconsistency in input. Problem 3D but no z-coordinates.');
        errFlag = 1;
    end
end

%check whether a priori Kalman filter information is given
if nargin < 6 || isempty(filterInfoPrev)
    filterInfoPrev = [];
    usePriorInfo = 0;
else
    usePriorInfo = 1;
end

%check whether previous costs have been input
if nargin < 7 || isempty(prevCost)
    prevCost = [];
end

%check whether verbose
if nargin < 8 || isempty(verbose)
    verbose = 1;
end

%exit if there are problems with input
if errFlag
    disp('--linkFeaturesKalmanSparse: Please fix input parameters.');
    return
end

%% preamble

%get number of frames in movie
numFrames = length(movieInfo);

%get number of features in each frame
if ~isfield(movieInfo,'num')
    for iFrame = 1 : numFrames
        movieInfo(iFrame).num = size(movieInfo(iFrame).xCoord,1);
    end
end

%collect coordinates and their std in one matrix in each frame
if ~isfield(movieInfo,'allCoord')
    if probDim == 2
        for iFrame = 1 : numFrames
            movieInfo(iFrame).allCoord = [movieInfo(iFrame).xCoord ...
                movieInfo(iFrame).yCoord];
        end
    else
        for iFrame = 1 : numFrames
            movieInfo(iFrame).allCoord = [movieInfo(iFrame).xCoord ...
                movieInfo(iFrame).yCoord movieInfo(iFrame).zCoord];
        end
    end
end

%calculate nearest neighbor distance for each feature in each frame
if ~isfield(movieInfo,'nnDist')

    for iFrame = 1 : numFrames
        
        switch movieInfo(iFrame).num

            case 0 %if there are no features

                %there are no nearest neighbor distances
                nnDist = zeros(0,1);

            case 1 %if there is only 1 feature

                %assign nearest neighbor distance as 1000 pixels (a very big
                %number)
                nnDist = 1000;

            otherwise %if there is more than 1 feature

                %compute distance matrix
                nnDist = createDistanceMatrix(movieInfo(iFrame).allCoord(:,1:2:end),...
                    movieInfo(iFrame).allCoord(:,1:2:end));

                %sort distance matrix and find nearest neighbor distance
                nnDist = sort(nnDist,2);
                nnDist = nnDist(:,2);

        end

        %store nearest neighbor distance
        movieInfo(iFrame).nnDist = nnDist;

    end

end

%% Linking

%make an array of the number of features per frame
numFeatures = zeros(numFrames,1);
for iFrame = 1 : numFrames
    numFeatures(iFrame) = movieInfo(iFrame).num;
end

%reserve memory for kalmanFilterInfo
if selfAdaptive

    % -- USER DEFINED FUNCTION -- %
    eval(['kalmanFilterInfo = ' kalmanFunctions.reserveMem ...
        '(numFrames,numFeatures,probDim);']);
    
else
    
    kalmanFilterInfo = zeros(numFrames,1);

end

%fill the feature indices in 1st frame in the connectivity matrix
trackedFeatureIndx = (1:movieInfo(1).num)';

%fill the nearest neighbor distances of features in first frame
nnDistFeatures = movieInfo(1).nnDist;

%initialize Kalman filter for features in 1st frame
if selfAdaptive

    % -- USER DEFINED FUNCTION -- %
    if usePriorInfo %use a priori information if available
        kalmanFilterInfo(1).stateVec = filterInfoPrev(1).stateVec; %state vector
        kalmanFilterInfo(1).stateCov = filterInfoPrev(1).stateCov; %state covariance
        kalmanFilterInfo(1).noiseVar = filterInfoPrev(1).noiseVar; %noise variance
    else
        eval(['[filterInit,errFlag] = ' kalmanFunctions.initialize ...
            '(movieInfo(1),probDim,costMatParam);'])
        kalmanFilterInfo(1).stateVec = filterInit.stateVec;
        kalmanFilterInfo(1).stateCov = filterInit.stateCov;
        kalmanFilterInfo(1).noiseVar = filterInit.noiseVar;
    end
    
end

%get number of particles in whole movie and calculate a worst-case scenario
%number of tracks
%it can be that the final number of tracks is even larger than this worst
%case scenario. Every time the auxiliary matrices (defined below) run out
%of rows, another "numTracksWorstCase" rows are added to them.
numTracksWorstCase = max(1,round(sum(numFeatures)/10));

%initialize auxiliary matrices for storing information related to tracks
%that end in the middle of the movie
trackedFeatureIndxAux = zeros(numTracksWorstCase,numFrames);
nnDistFeaturesAux = NaN(numTracksWorstCase,numFrames);
prevCostAux = NaN(numTracksWorstCase,numFrames);
rowEnd = numTracksWorstCase;

%store the costs of previous links for features in first frame
%if no previous costs have been input
%in this case, store NaN since there are no previous links
if isempty(prevCost)
    prevCost = NaN(movieInfo(1).num,1);
else
    prevCost = max(prevCost(:))*ones(movieInfo(1).num,1);
end
prevCostStruct.all = prevCost;
prevCostStruct.max = max(prevCost(:));
% prevCostStruct.allAux = [];
prevCostStruct.allAux = HandleObject(prevCostAux); % Updated by Carmen Klein Herenbrink and Brian Devree

%assign the lifetime of features in first frame
featLifetime = ones(movieInfo(1).num,1);

% % % %for paper - get number of potential link per feature
% % % numPotLinksPerFeature = [];

%initialize progress display
if verbose
    progressText(0,'Linking frame-to-frame');
end

%go over all frames
for iFrame = 1 : numFrames-1

    %get number of features
    numFeaturesFrame1 = movieInfo(iFrame).num; %in 1st frame
    numFeaturesFrame2 = movieInfo(iFrame+1).num; %in 2nd frame

    if numFeaturesFrame1 ~= 0 %if there are features in 1st frame
        
        if numFeaturesFrame2 ~= 0 %if there are features in 2nd frame
            
            %calculate cost matrix
            % -- USER DEFINED FUNCTION -- %
            % eval(['[costMat,propagationScheme,kalmanFilterInfoTmp,nonlinkMarker]'...
            %     ' = ' costMatName '(movieInfo,kalmanFilterInfo(iFrame),'...
            %     'costMatParam,nnDistFeatures(1:numFeaturesFrame1,:),'...
            %     'probDim,prevCostStruct,featLifetime,trackedFeatureIndx,iFrame);'])
            [costMat,propagationScheme,kalmanFilterInfoTmp,nonlinkMarker] = costMatFunc(...
                movieInfo,kalmanFilterInfo(iFrame),costMatParam,nnDistFeatures(1:numFeaturesFrame1,:),...
                probDim,prevCostStruct,featLifetime,trackedFeatureIndx,iFrame); % Updated by Carmen Klein Herenbrink and Brian Devree

            % % %             %for paper - get number of potential links per feature
            % % %             numPotLinksPerFeature = [numPotLinksPerFeature; sum(...
            % % %                 costMat(1:numFeaturesFrame1,1:numFeaturesFrame2)...
            % % %                 ~=nonlinkMarker,2)];

            if any(costMat(:)~=nonlinkMarker) %if there are potential links

                %link features based on cost matrix, allowing for birth and death
                [link12,link21] = lap(costMat,nonlinkMarker,0);

                %get indices of features in 2nd frame that are connected to features in 1st frame
                indx2C = find(link21(1:numFeaturesFrame2)<=numFeaturesFrame1);

                %get indices of corresponding features in 1st frame
                indx1C = link21(indx2C);

                %find existing tracks that are not connected to features in 2nd frame
                numExistTracks = size(trackedFeatureIndx,1);
                indx1U = setdiff(1:numExistTracks,indx1C);
                numRows = length(indx1U);
                
                %determine where to store these tracks in auxiliary matrix
                %extend auxiliary matrices if necessary
                rowStart = rowEnd - numRows + 1;
                while(rowStart <= 1)
                    trackedFeatureIndxAux = [zeros(numTracksWorstCase,numFrames); ...
                        trackedFeatureIndxAux];
                    nnDistFeaturesAux = [NaN(numTracksWorstCase,numFrames); ...
                        nnDistFeaturesAux];
                    prevCostAux = [NaN(numTracksWorstCase,numFrames); ...
                        prevCostAux];
                    rowEnd = rowEnd + numTracksWorstCase;
                    rowStart = rowStart + numTracksWorstCase;
                end
                
                %move rows of tracks that are not connected to points in
                %2nd frame to auxilary matrix
                trackedFeatureIndxAux(rowStart:rowEnd,1:iFrame) = trackedFeatureIndx(indx1U,:);
                
                %assign space for new connectivity matrix
                tmp = zeros(numFeaturesFrame2,iFrame+1);

                %fill in the feature numbers in 2nd frame
                tmp(1:numFeaturesFrame2,iFrame+1) = (1:numFeaturesFrame2)';

                %shuffle existing tracks to get the correct connectivity with 2nd frame
                tmp(indx2C,1:iFrame) = trackedFeatureIndx(indx1C,:);
                
                %update the connectivity matrix "trackedFeatureIndx"
                trackedFeatureIndx = tmp;

                %repeat for the matrix of nearest neighbor distances
                nnDistFeaturesAux(rowStart:rowEnd,1:iFrame) = nnDistFeatures(indx1U,:);
                tmp = NaN(numFeaturesFrame2,iFrame+1);
                tmp(1:numFeaturesFrame2,iFrame+1) = movieInfo(iFrame+1).nnDist;
                tmp(indx2C,1:iFrame) = nnDistFeatures(indx1C,:);
                nnDistFeatures = tmp;

                %repeat for the matrix of linking costs
                prevCostAux(rowStart:rowEnd,1:iFrame) = prevCost(indx1U,:);
                tmp = NaN(numFeaturesFrame2,iFrame+1);
                for i = 1 : length(indx2C)
                    tmp(indx2C(i),iFrame+1) = costMat(indx1C(i),indx2C(i));
                end
                tmp(indx2C,1:iFrame) = prevCost(indx1C,:);
                prevCost = tmp;
                
                %update rowEnd to indicate until which row the auxiliary
                %matrices are ampty
                rowEnd = rowStart - 1;

                %get track lifetimes for features in 2nd frame
                featLifetime = ones(numFeaturesFrame2,1);
                for iFeat = 1 : numFeaturesFrame2
                    featLifetime(iFeat) = length(find(trackedFeatureIndx(iFeat,:)~=0));
                end

                %use the Kalman gain from linking to get better estimates
                %of the state vector and its covariance matrix in 2nd frame
                %as well as state noise and its variance
                if selfAdaptive

                    % -- USER DEFINED FUNCTION -- %
                    if usePriorInfo %if prior information is supplied
                        eval(['[kalmanFilterInfo,errFlag] = ' kalmanFunctions.calcGain ...
                            '(trackedFeatureIndx(1:numFeaturesFrame2,:),'...
                            'movieInfo(iFrame+1),kalmanFilterInfoTmp,'...
                            'propagationScheme,kalmanFilterInfo,probDim,'...
                            'filterInfoPrev(iFrame+1),costMatParam,kalmanFunctions.initialize);'])
                    else %if no prior information is supplied
                        eval(['[kalmanFilterInfo,errFlag] = ' kalmanFunctions.calcGain ...
                            '(trackedFeatureIndx(1:numFeaturesFrame2,:),'...
                            'movieInfo(iFrame+1),kalmanFilterInfoTmp,'...
                            'propagationScheme,kalmanFilterInfo,probDim,'...
                            '[],costMatParam,kalmanFunctions.initialize);'])
                    end
                    
                end
                
            else %if there are no potential links
                
                %determine where to store the tracks up to 1st frame in
                %auxiliary matrix
                %extend auxiliary matrices if necessary
                numRows = size(trackedFeatureIndx,1);
                rowStart = rowEnd - numRows + 1;
                while( rowStart <= 1)
                    trackedFeatureIndxAux = [zeros(numTracksWorstCase,numFrames); ...
                        trackedFeatureIndxAux];
                    nnDistFeaturesAux = [NaN(numTracksWorstCase,numFrames); ...
                        nnDistFeaturesAux];
                    prevCostAux = [NaN(numTracksWorstCase,numFrames); ...
                        prevCostAux];
                    rowEnd = rowEnd + numTracksWorstCase;
                    rowStart = rowStart + numTracksWorstCase;
                end
                
                %move tracks upto 1st frame to auxiliary matrix
                trackedFeatureIndxAux(rowStart:rowEnd,1:iFrame) = trackedFeatureIndx;

                %assign space for new connectivity matrix
                trackedFeatureIndx = zeros(numFeaturesFrame2,iFrame+1);

                %fill in the feature numbers in 2nd frame
                trackedFeatureIndx(1:numFeaturesFrame2,iFrame+1) = (1:numFeaturesFrame2)';
                
                %repeat for the matrix of nearest neighbor distances
                nnDistFeaturesAux(rowStart:rowEnd,1:iFrame) = nnDistFeatures;
                nnDistFeatures = NaN(numFeaturesFrame2,iFrame+1);
                nnDistFeatures(1:numFeaturesFrame2,iFrame+1) = movieInfo(iFrame+1).nnDist;

                %repeat for the matrix of linking costs
                prevCostAux(rowStart:rowEnd,1:iFrame) = prevCost;
                prevCost = NaN(numFeaturesFrame2,iFrame+1);
                
                %update rowEnd to indicate until which row the auxiliary
                %matrices are ampty
                rowEnd = rowStart - 1;

                %assign track lifetimes for features in 2nd frame
                featLifetime = ones(numFeaturesFrame2,1);

                %initialize Kalman filter for features in 2nd frame
                if selfAdaptive

                    % -- USER DEFINED FUNCTION -- %
                    if usePriorInfo %use a priori information if available
                        kalmanFilterInfo(iFrame+1).stateVec = filterInfoPrev(iFrame+1).stateVec; %state vector
                        kalmanFilterInfo(iFrame+1).stateCov = filterInfoPrev(iFrame+1).stateCov; %state covariance
                        kalmanFilterInfo(iFrame+1).noiseVar = filterInfoPrev(iFrame+1).noiseVar; %noise variance
                    else
                        eval(['[filterInit,errFlag] = ' kalmanFunctions.initialize ...
                            '(movieInfo(iFrame+1),probDim,costMatParam);'])
                        kalmanFilterInfo(iFrame+1).stateVec = filterInit.stateVec;
                        kalmanFilterInfo(iFrame+1).stateCov = filterInit.stateCov;
                        kalmanFilterInfo(iFrame+1).noiseVar = filterInit.noiseVar;
                    end
                    
                end

            end %(if any(costMat(:)~=nonlinkMarker))
            
        else %if there are no features in 2nd frame
            
            %determine where to store the tracks up to 1st frame in
            %auxiliary matrix
            %extend auxiliary matrices if necessary
            numRows = size(trackedFeatureIndx,1);
            rowStart = rowEnd - numRows + 1;
            if rowStart <= 1
                trackedFeatureIndxAux = [zeros(numTracksWorstCase,numFrames); ...
                    trackedFeatureIndxAux];
                nnDistFeaturesAux = [NaN(numTracksWorstCase,numFrames); ...
                    nnDistFeaturesAux];
                prevCostAux = [NaN(numTracksWorstCase,numFrames); ...
                    prevCostAux];
                rowEnd = rowEnd + numTracksWorstCase;
                rowStart = rowStart + numTracksWorstCase;
            end
            
            %move tracks upto 1st frame to auxiliary matrix
            trackedFeatureIndxAux(rowStart:rowEnd,1:iFrame) = trackedFeatureIndx;
            
            %update the connectivity matrix "trackedFeatureIndx"
            trackedFeatureIndx = zeros(numFeaturesFrame2,iFrame+1);
            
            %repeat for the matrix of nearest neighbor distances
            nnDistFeaturesAux(rowStart:rowEnd,1:iFrame) = nnDistFeatures;
            nnDistFeatures = NaN(numFeaturesFrame2,iFrame+1);
            
            %repeat for the matrix of linking costs
            prevCostAux(rowStart:rowEnd,1:iFrame) = prevCost;
            prevCost = NaN(numFeaturesFrame2,iFrame+1);
            
            %update rowEnd to indicate until which row the auxiliary
            %matrices are ampty
            rowEnd = rowStart - 1;
            
            %assign track lifetimes for features in 2nd frame
            featLifetime = [];
            
        end %(if numFeaturesFrame2 ~= 0 ... else ...)
        
    else %if there are no features in 1st frame
        
        if numFeaturesFrame2 ~= 0 %if there are features in 2nd frame
            
            %assign space for new connectivity matrix
            trackedFeatureIndx = zeros(numFeaturesFrame2,iFrame+1);
            
            %fill in the feature numbers in 2nd frame
            trackedFeatureIndx(1:numFeaturesFrame2,iFrame+1) = (1:numFeaturesFrame2)';
            
            %repeat for the matrix of nearest neighbor distances
            nnDistFeatures = NaN(numFeaturesFrame2,iFrame+1);
            nnDistFeatures(1:numFeaturesFrame2,iFrame+1) = movieInfo(iFrame+1).nnDist;
            
            %repeat for the matrix of linking costs
            prevCost = NaN(numFeaturesFrame2,iFrame+1);
            
            %assign track lifetimes for features in 2nd frame
            featLifetime = ones(numFeaturesFrame2,1);
            
            %initialize Kalman filter for features in 2nd frame
            if selfAdaptive
                
                % -- USER DEFINED FUNCTION -- %
                if usePriorInfo %use a priori information if available
                    kalmanFilterInfo(iFrame+1).stateVec = filterInfoPrev(iFrame+1).stateVec; %state vector
                    kalmanFilterInfo(iFrame+1).stateCov = filterInfoPrev(iFrame+1).stateCov; %state covariance
                    kalmanFilterInfo(iFrame+1).noiseVar = filterInfoPrev(iFrame+1).noiseVar; %noise variance
                else
                    eval(['[filterInit,errFlag] = ' kalmanFunctions.initialize ...
                        '(movieInfo(iFrame+1),probDim,costMatParam);'])
                    kalmanFilterInfo(iFrame+1).stateVec = filterInit.stateVec;
                    kalmanFilterInfo(iFrame+1).stateCov = filterInit.stateCov;
                    kalmanFilterInfo(iFrame+1).noiseVar = filterInit.noiseVar;
                end
                
            end

        else %if there are no features in 2nd frame

            %assign space for new connectivity matrix
            trackedFeatureIndx = zeros(numFeaturesFrame2,iFrame+1);
            
            %repeat for the matrix of nearest neighbor distances
            nnDistFeatures = NaN(numFeaturesFrame2,iFrame+1);
            
            %repeat for the matrix of linking costs
            prevCost = NaN(numFeaturesFrame2,iFrame+1);
            
            %assign track lifetimes for features in 2nd frame
            featLifetime = [];

        end %(if numFeaturesFrame2 ~= 0 ... else ...)

    end %(if numFeaturesFrame1 ~= 0 ... else ...)

    %update structure of previous costs
    prevCostStruct.all = prevCost;
    prevCostStruct.max = max([prevCostStruct.max; prevCost(:,end)]);
    % prevCostStruct.allAux = prevCostAux;
    % prevCostStruct.allAux = HandleObject(prevCostAux); % Updated by Carmen Klein Herenbrink and Brian Devree
    
    %display progress
    if verbose
        progressText(iFrame/(numFrames-1),'Linking frame-to-frame');
    end
    
end %(for iFrame=1:numFrames-1)

%add information from last frame to auxiliary matrices
numRows = size(trackedFeatureIndx,1);
rowStart = rowEnd - numRows + 1;
if rowStart <= 1
    trackedFeatureIndxAux = [zeros(numRows,numFrames); ...
        trackedFeatureIndxAux];
    nnDistFeaturesAux = [NaN(numRows,numFrames); ...
        nnDistFeaturesAux];
    prevCostAux = [NaN(numRows,numFrames); ...
        prevCostAux];
    rowEnd = rowEnd + numRows;
    rowStart = rowStart + numRows;
end
trackedFeatureIndxAux(rowStart:rowEnd,:) = trackedFeatureIndx;
nnDistFeaturesAux(rowStart:rowEnd,:) = nnDistFeatures;
prevCostAux(rowStart:rowEnd,:) = prevCost;

%remove all empty rows
trackedFeatureIndx = trackedFeatureIndxAux(rowStart:end,:);
clear trackedFeatureIndxAux
nnDistFeatures = nnDistFeaturesAux(rowStart:end,:);
clear nnDistFeaturesAux
prevCost = prevCostAux(rowStart:end,:);
clear prevCostAux

%get total number of tracks
numTracks = size(trackedFeatureIndx,1);

%find the frame where each track begins and then sort the vector
frameStart = zeros(numTracks,1);
for i=1:numTracks
    frameStart(i) = find((trackedFeatureIndx(i,:)~=0),1,'first');
end
[frameStart,indx] = sort(frameStart);

%rearrange "trackedFeatureIndx" such that tracks are sorted in ascending order by their
%starting point. Note that this ends up also arranging tracks starting at the 
%same time in descending order from longest to shortest.
trackedFeatureIndx = trackedFeatureIndx(indx,:);

%also re-arrange the matrix indicating nearest neighbor distances and
%previouc costs
nnDistFeatures = nnDistFeatures(indx,:);
prevCost = prevCost(indx,:);

%clear some memory
clear costMat tmp

%store feature positions and amplitudes in a matrix that also shows their connectivities
%information is stored as [x y z a dx dy dz da] in image coordinate system
%trackedFeatureInfo is in sparse format
trackedFeatureInfo = coordAmpMatFromIndicesSparse(trackedFeatureIndx,movieInfo,...
    numFrames,probDim);

%take absolute value of all noise variances - this takes care of the
%negative variances used to indicate first appearances
if selfAdaptive && ~usePriorInfo
    for iFrame = 1 : numFrames
        kalmanFilterInfo(iFrame).noiseVar = abs(kalmanFilterInfo(iFrame).noiseVar);
    end
end


%% %%%%% ~~ the end ~~ %%%%%

