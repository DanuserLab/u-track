function [costMat,propagationScheme,kalmanFilterInfoFrame2,nonlinkMarker,...
    errFlag] = plusTipCostMatLinearMotionLink(movieInfo,kalmanFilterInfoFrame1,...
    costMatParam,nnDistFeatures,probDim,prevCost,featLifetime,...
    trackedFeatureIndx,currentFrame)
%COSTMATLINEARMOTIONLINK provides a cost matrix for linking features based on competing linear motion models
%
%SYNOPSIS [costMat,propagationScheme,kalmanFilterInfoFrame2,nonlinkMarker,...
%     errFlag] = costMatLinearMotionLink(movieInfo,kalmanFilterInfoFrame1,...
%     costMatParam,nnDistFeatures,probDim,prevCost,featLifetime,...
%     trackedFeatureIndx,currentFrame)
%
%INPUT  movieInfo             : An nx1 array (n = number of frames in
%                               movie) containing the fields:
%             .allCoord           : x,dx,y,dy,[z,dz] of features collected in one
%                                   matrix.
%             .amp                : Amplitudes of PSFs fitting detected features.
%                                   1st column for values and 2nd column
%                                   for standard deviations.
%             .num                : Number of features in each frame.
%             .nnDist             : Distance from each feature to its nearest
%                                   neighbor. Not needed at the moment.
%      kalmanFilterInfoFrame1 : Structure with at least the following fields:
%             .stateVec           : Kalman filter state vector for each
%                                   feature in 1st frame.
%             .stateCov           : Kalman filter state covariance matrix
%                                   for each feature in 1st frame.
%             .noiseVar           : Variance of state noise for each
%                                   feature in 1st frame.
%      costMatParam           : Structure with fields:
%             .linearMotion       : 1 to propagate using a linear motion
%                                   model, 0 otherwise.
%             .minSearchRadius    : Minimum allowed search radius.
%             .maxSearchRadius    : Maximum allowed search radius.
%             .brownStdMult       : Factor multiplying Brownian
%                                   displacement std to get search radius.
%             .lftCdf             : Lifetime cumulative density function.
%                                   Column vector, specifying cdf for
%                                   lifetime = 0 to movie length.
%                                   Enter [] if cdf is not to be used.
%                                   Optional. Default: [].
%             .useLocalDensity    : Logical variable indicating whether to use
%                                   local density in search radius estimation.
%             .nnWindow           : Number of past frames for calculating
%                                   nearest neighbor distance.
%      nnDistFeatures         : Matrix of nearest neighbor distances of
%                               features in first frame as well as of
%                               features in previous frames that they are
%                               connected to.
%      probDim                : Problem dimensionality. 2 (for 2D) or 3 (for 3D).
%      prevCost               : Matrix of previous linking costs.
%      featLifetime           : Lengths of tracks that features in
%                               first frame belong to.
%      trackedFeatureIndx     : The matrix of feature index connectivity up
%                               to current frame.
%                               Currently not used in this cost function.
%      currentFrame           : Current frame that is being linked to the
%                               next frame.
%
%OUTPUT costMat               : Cost matrix.
%       propagationScheme     : Propagation scheme corresponding to each
%                               cost in the cost matrix.
%       kalmanFilterInfoFrame2: Structure with at least the following fields:
%             .stateVec           : Kalman filter prediction of state
%                                   vector in 2nd frame based on all 3
%                                   motion models.
%             .stateCov           : Kalman filter prediction of state
%                                   covariance matrix in 2nd frame based on
%                                   all 3 motion models.
%             .obsVec             : Kalman filter prediction of the
%                                   observed variables in 2nd frame based
%                                   on all 3 motion models.
%       nonlinkMarker         : Value indicating that a link is not allowed.
%       errFlag               : 0 if function executes normally, 1 otherwise.
%
%REMARKS Three competing linear motion models: 1, 2 and 3.
%        1: forward drift, 2: backward drift, 3: zero drift (Brownian).
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

costMat = [];
propagationScheme = [];
kalmanFilterInfoFrame2 = [];
nonlinkMarker = [];
errFlag = [];

%% Input

%check whether correct number of input arguments was used
if nargin ~= nargin('plusTipCostMatLinearMotionLink')
    disp('--plusTipCostMatLinearMotionLink: Incorrect number of input arguments!');
    errFlag  = 1;
    return
end

%get cost function parameters
linearMotion = costMatParam.linearMotion;
minSearchRadius = costMatParam.minSearchRadius;
maxSearchRadius = costMatParam.maxSearchRadius;
brownStdMult    = costMatParam.brownStdMult;
useLocalDensity = costMatParam.useLocalDensity;
nnWindow = costMatParam.nnWindow;
if useLocalDensity
    closestDistScale = 2;
    maxStdMult = 100;
end
if isfield('costMatParam','lftCdf')
    lftCdf = costMatParam.lftCdf;
else
    lftCdf = [];
end
if isfield(costMatParam,'diagnostics')
    diagnostics = costMatParam.diagnostics;
else
    diagnostics = 0;
end

%calculate nearest neighbor distance given feature history
frameNum = size(nnDistFeatures,2);
tmpNN = max(1,frameNum-nnWindow);
nnDistTracks = min(nnDistFeatures(:,tmpNN:end),[],2);

%extract the two frames of interest from movieInfo
movieInfo = movieInfo(currentFrame:currentFrame+1);

%find features that appear and remove negative signs
firstAppearanceFlag = sign(squeeze(kalmanFilterInfoFrame1.noiseVar(1,1,:)));
kalmanFilterInfoFrame1.noiseVar = abs(kalmanFilterInfoFrame1.noiseVar);

%% Motion propagation

%specify number of propagation schemes used
numSchemes = 3;

%calculate vector sizes
vecSize = 2 * probDim;

%construct transition matrices
if linearMotion
    transMat(:,:,1) = eye(vecSize) + diag(ones(probDim,1),probDim); %forward drift transition matrix
    transMat(:,:,2) = eye(vecSize) + diag(ones(probDim,1),probDim); %forward drift transition matrix
    transMat(:,:,3) = eye(vecSize) + diag(ones(probDim,1),probDim); %forward drift transition matrix
    %transMat(:,:,2) = eye(vecSize) + diag(-ones(probDim,1),probDim); %backward drift transition matrix
    %transMat(:,:,3) = eye(vecSize); %zero drift transition matrix
else
    transMat(:,:,3) = eye(vecSize) + diag(ones(probDim,1),probDim); %forward drift transition matrix
    transMat(:,:,2) = eye(vecSize) + diag(-ones(probDim,1),probDim); %backward drift transition matrix
    transMat(:,:,1) = eye(vecSize); %zero drift transition matrix
end

%construct observation matrix
observationMat = [eye(probDim) zeros(probDim)]; %observation matrix

%get number of features in the 2 frames
numFeaturesFrame1 = movieInfo(1).num;
numFeaturesFrame2 = movieInfo(2).num;

%reserve memory for "kalmanFilterInfoframe2"
kalmanFilterInfoFrame2 = struct('stateVec',zeros(numFeaturesFrame1,vecSize,numSchemes),...
    'stateCov',zeros(vecSize,vecSize,numFeaturesFrame1,numSchemes),...
    'obsVec',zeros(numFeaturesFrame1,probDim,numSchemes));

%apply Kalman filters to each feature in 1st frame
for iFeature = 1 : numFeaturesFrame1
    
    %get state vector and its covariance matrix of feature in 1st frame
    stateOld = kalmanFilterInfoFrame1.stateVec(iFeature,:)';
    stateCovOld = kalmanFilterInfoFrame1.stateCov(:,:,iFeature);
    noiseVar = kalmanFilterInfoFrame1.noiseVar(:,:,iFeature);
    
    %go over all possible propagation schemes
    for iScheme = 1 : numSchemes
        
        %predict state vector of feature in 2nd frame
        stateVec = transMat(:,:,iScheme)*stateOld;
        
        %predict state covariance matrix of feature in 2nd frame
        stateCov = transMat(:,:,iScheme)*stateCovOld*transMat(:,:,iScheme)' ...
            + noiseVar;
        
        %determine observation vector of feature in 2nd frame (i.e. the
        %propagated position of the feature)
        obsVec = observationMat*stateVec;
        
        %save information in kalmanFilterInfoFrame2
        kalmanFilterInfoFrame2.stateVec(iFeature,:,iScheme) = stateVec';
        kalmanFilterInfoFrame2.stateCov(:,:,iFeature,iScheme) = stateCov;
        kalmanFilterInfoFrame2.obsVec(iFeature,:,iScheme) = obsVec';
        
    end
    
end

%get the propagated positions of features in 1st frame based on the three propagation schemes
propagatedPos = kalmanFilterInfoFrame2.obsVec;

%put the coordinates of features in the 2nd frame in one matrix
coord2 = movieInfo(2).allCoord(:,1:2:end);

%calculate the cost matrices for all three propagation schemes
for iScheme = 1 : numSchemes
    
    %put the propagated x and y coordinates of features from 1st frame in
    %one matrix
    coord1 = propagatedPos(:,:,iScheme);
    
    %calculate the distances between features
    costMatTmp(:,:,iScheme) = createDistanceMatrix(coord1,coord2);
    
end

%find the minimum cost for the link between every pair, which also
%determines the best propagation scheme to perform that link
[costMat,propagationScheme] = min(costMatTmp,[],3);

%% Search radius

%get the Kalman standard deviation of all features in frame 1
kalmanStd = sqrt(probDim * squeeze(kalmanFilterInfoFrame1.noiseVar(1,1,:)));

%copy brownStdMult into vector
stdMultInd = repmat(brownStdMult,numFeaturesFrame1,1);

%if local density information is used to expand search radius ...
if useLocalDensity
    
    %divide each feature's nearest neighbor distance/closestDistScale by kalmanStd
    ratioDist2Std = nnDistTracks./kalmanStd/closestDistScale;
    
    %make ratios larger than maxStdMult equal to maxStdMult
    ratioDist2Std(ratioDist2Std > maxStdMult) = maxStdMult;
    
    %expand search radius multiplication factor if possible
    stdMultInd = max([stdMultInd ratioDist2Std],[],2);
    
end

%get the search radius of each feature in frame 1 and make sure it falls
%within reasonable limits
searchRadius = stdMultInd .* kalmanStd;

% if it is the first frame, we don't bound the searchRadius with min/max
% values because it essentially does a nearest-neighbor search.  in the
% subsequent iterations, we can bound to min/max because we want to look
% around the propagated positions for features to link

oldIdx = firstAppearanceFlag==1; % these are not new features - apply min/max
tooBigIdx = searchRadius>maxSearchRadius;
tooSmallIdx = searchRadius<minSearchRadius;

searchRadius(oldIdx & tooBigIdx)   = maxSearchRadius;
searchRadius(oldIdx & tooSmallIdx) = minSearchRadius;

%replicate the search radius to compare to cost matrix
searchRadius = repmat(searchRadius,1,numFeaturesFrame2);

%assign NaN to costs corresponding to distance > searchRadius
costMat(costMat>searchRadius) = NaN;

%square the cost matrix to make the cost = distance squared
costMat = costMat.^2;

%% Lifetime penalty

if ~isempty(lftCdf)
    
    %specify 1 - lifetime cumulative probability
    oneMinusLftCdf = 1 - lftCdf;
    
    %calculate 1 / (lifetime penalty), which is 1 / (1-cumulative probability
    %of lifetime of feature in first frame)
    oneOverLftPen = oneMinusLftCdf(featLifetime+1);
    
    %multiple each cost by the lifetime penalty
    costMat = costMat ./ repmat(oneOverLftPen,1,numFeaturesFrame2);
    
    %replace infinite costs by NaN
    costMat(isinf(costMat)) = NaN;
    
end

%% Birth and death

%append matrix to allow birth and death
if isstruct(prevCost)
    prevCostMax = prevCost.max;
else
    prevCostMax = max(prevCost(:));
end

if ~isnan(prevCostMax) && prevCostMax ~= 0
    maxCost = 1.05*prevCostMax;
else
    maxCost = max(prctile(costMat(:),80),eps);
end
deathCost = maxCost * ones(numFeaturesFrame1,1);
birthCost = maxCost * ones(numFeaturesFrame2,1);

%generate upper right and lower left block
deathBlock = diag(deathCost); %upper right
deathBlock(deathBlock==0) = NaN;
birthBlock = diag(birthCost); %lower left
birthBlock(birthBlock==0) = NaN;

%get the cost for the lower right block
costLR = min(min(min(costMat))-1,-1);
lrBlock = costMat';
lrBlock(~isnan(lrBlock)) = costLR;

%append cost matrix
costMat = [costMat deathBlock; birthBlock lrBlock];

%% nonLinkMarker

%determine the nonlinkMarker
nonlinkMarker = min(floor(min(min(costMat)))-5,-5);

%replace NaN, indicating pairs that cannot be linked, with nonlinkMarker
costMat(isnan(costMat)) = nonlinkMarker;



%% Histogram of linking distances

%get current frame
if isstruct(prevCost)
    currentFrame = size(prevCost.all,2);
else
    currentFrame = size(prevCost,2);
end
%check whether current frame matches any of the diagnostics frames
if currentFrame ~= 1 && any(diagnostics == currentFrame)
    
    %get linking distances
    if isstruct(prevCost)
        prevCostNoCol1 = prevCost.all(:,2:end);
    else
        prevCostNoCol1 = prevCost(:,2:end);
    end
    linkingDistances = sqrt(prevCostNoCol1(~isnan(prevCostNoCol1)));
    
    %plot histogram
    figure('Name',['frame # ' num2str(currentFrame)],'NumberTitle','off');
    try
        %         % create bins
        %         n=1:20;
        %         % bin the samples
        %         [x1,nbins1] = histc(linkingDistances,n); % forward
        %         % make the plot
        %         figure
        %         bar(n,x1,'stack')
        optimalHistogram(linkingDistances,[],0);
        xlabel('Linking distance');
        ylabel('Counts');
    catch
        disp('histogram plot failed');
    end
    
end
%% ~~~ the end ~~~


%% old snippets of code

%             .cost2probExpo      : Exponent coefficient in formula
%                                   converting linking cost to scheme
%                                   probability.
% cost2probExpo    = costMatParam.cost2probExpo;
% %reserve memory
% costMat = NaN*ones(numFeaturesFrame1,numFeaturesFrame2);
% propagationScheme = ones(numFeaturesFrame1,numFeaturesFrame2);
%
% %pick a set of uniformly distributed random numbers
% randomNumber = rand(numFeaturesFrame1,numFeaturesFrame2);
%
% %go over all features in frame 1
% for i1 = 1 : numFeaturesFrame1
%
%     %get the square search radius of feature i1 in frame 1
%     %     searchRadius2 = min(noiseStdMult2*2*...
%     %         kalmanFilterInfoFrame1.noiseVar(1,1,i1),maxSearchRadius2);
%     searchRadius2 = min(noiseStdMult2*2*...
%         max(kalmanFilterInfoFrame1.noiseVar(1,1,:)),maxSearchRadius2);
%
%     %go over all features in frame 2
%     for i2 = 1 : numFeaturesFrame2
%
%         %get the costs for this pair
%         costPair = squeeze(costMatTmp(i1,i2,:));
%
%         %assign NaN to schemes whose cost = square distance > searchRadius2
%         costPair(costPair>searchRadius2) = NaN;
%
%         %find the propagation scheme with the minimum cost
%         [costTmp,propSchemeTmp] = min(costPair);
%
%         %if there is a possible link between the two features
%         if ~isnan(costTmp)
%
%             %if Brownian motion has the minumum cost, accept it
%             if propSchemeTmp == 3
%
%                 propagationScheme(i1,i2) = propSchemeTmp;
%                 costMat(i1,i2) = costTmp;
%
%             else %if forward or backward propagation has minimum cost
%
%                 %make the other propagation scheme impossible
%                 if propSchemeTmp == 1
%                     costPair(2) = NaN;
%                 else
%                     costPair(1) = NaN;
%                 end
%
%                 %calculate the probability of each scheme using its cost
%                 schemeProb = exp(-cost2probExpo*costPair);
%
%                 %assign a probability of zero to impossible schemes
%                 schemeProb(isnan(schemeProb)) = 0;
%
%                 %normalize the probabilities
%                 schemeProb = schemeProb/sum(schemeProb);
%
%                 %calculate the cumulative probability
%                 schemeProb(2) = schemeProb(2) + schemeProb(1);
%                 schemeProb(3) = 1;
%
%                 %choose propagation scheme based on probabilities and random number
%                 propagationScheme(i1,i2) = find((randomNumber(i1,i2)<=schemeProb),1,'first');
%
%                 %assign cost based on chosen propagation scheme
%                 costMat(i1,i2) = costMatTmp(i1,i2,propagationScheme(i1,i2));
%
%             end %(if propSchemeTmp == 3)
%
%         end %(if any(~isnan(costPair)))
%
%     end %(for i2 = 1 : numFeaturesFrame2)
% end %(for i1 = 1 : numFeaturesFrame1)
