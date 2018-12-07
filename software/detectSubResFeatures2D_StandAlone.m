function [movieInfo,exceptions,localMaxima,background,psfSigma] = ...
    detectSubResFeatures2D_StandAlone(movieParam,detectionParam,saveResults,verbose)
%DETECTSUBRESFEATURES2D_STANDALONE detects subresolution features in a series of images
%
%SYNOPSIS [movieInfo,exceptions,localMaxima,background,psfSigma] = ...
%    detectSubResFeatures2D_StandAlone(movieParam,detectionParam,saveResults)
%
%INPUT  movieParam    : Structure with fields
%           .imageDir     : Directory where images are stored.
%           .filenameBase : Filename base.
%           .firstImageNum: Numerical index of first image in movie.
%           .lastImageNum : Numerical index of last image in movie.
%           .digits4Enum  : Number of digits used to enumerate frames.
%       detectionParam: Structure with fields
%           .psfSigma     : Initial guess for standard deviation of point
%                           spread function (in pixels).
%           .calcMethod   : 'g' for Gaussian fitting, 'c' for centroid
%                           calculation.
%                           Optional. Default: 'g'.
%           .testAlpha    : Alpha-values for statistical tests in
%                           detectSubResFeatures2D.
%                           (See detectSubResFeatures2D for details).
%                           Optional. Default values 0.05.s
%           .visual       : 1 if user wants to view results; 0 otherwise.
%                           Optional. Default: 0.
%           .doMMF        : 1 if user wants to do mixture-model fitting, 0
%                           otherwise.
%                           Optional. Default: 1.
%           .bitDepth     : Camera bit depth.
%                           Optional. Default: 16.
%           .alphaLocMax  : Alpha-value for statistical test comparing the
%                           amplitudes of local maxima to the local
%                           background.
%                           Optional. default: 0.05.
%                           --- alphaLocMax must be a row vector if
%                           integWindow is a row vector. See description of
%                           integWindow below.
%           .numSigmaIter : Maximum number of iterations to perform when
%                           trying to estimate PSF sigma. Input 0 for no
%                           estimation.
%                           Optional. Default: 10.
%           .integWindow  : Number of frames on each side of a frame
%                           used for time integration.
%                           Optional. Default: 0.
%                           --- integWindow can be a row vector, in which
%                           case alphaLocMax should be a row vector of the
%                           same length. When integWindow is a row vector,
%                           the initial local maxima detection is done by
%                           using all specified integration windows.
%           .background   : Structure with fields
%               .imageDir      : Directory where background images are
%                                stored.
%               .filenameBase  : Background images filename base.
%               .alphaLocMaxAbs: Alpha-value for statistical test comparing
%                                the amplitudes of local maxima to the
%                                absolute background. Same dimensions as
%                                alphaLocMax above.
%                           If this structure is supplied, the code will
%                           calculate a uniform background mean and std
%                           using the specified background images. The code
%                           assumes the same firstImageNum, lastImageNum
%                           and digits4Enum as the images to be analyzed.
%                           Optional. Default: Absolute background is not
%                           calculated and the significance of local maxima
%                           is assessed based on local background only.
%                           Omit field or assign as [] to use default.
%           .maskLoc      : Name (including full path) of mask file
%                           specifying ROI for detection.
%           .roiMask      : binary mask specifying ROI for detection.
%       saveResults   : 0 if no saving is requested.
%                       If saving is requested, structure with fields:
%           .dir          : Directory where results should be saved.
%                           Optional. Default: current directory.
%           .filename     : Name of file where results should be saved.
%                           Optional. Default: detectedFeatures.
%                       Whole structure optional.
%       verbose       : 1 to show progress while running, 0 otherwise.
%                       Optional. Default: 1.
%
%       All optional variables can be entered as [] to use default values.
%
%OUTPUT movieInfo     : Structure array of length = number of frames in
%                       movie, containing the fields:
%             .xCoord    : Image coordinate system x-coordinate of detected
%                          features [x dx] (in pixels).
%             .yCoord    : Image coordinate system y-coordinate of detected
%                          features [y dy] (in pixels).
%             .amp       : Amplitudes of PSFs fitting detected features [a da].
%       exceptions    : Structure with fields:
%             .emptyFrames: Array indicating frames where no features were
%                           detected.
%             .framesFailedMMF: Array indicating frames where mixture-model
%                               fitting failed.
%             .framesFailedLocMax: Array indicating frames where initial
%                                  detection of local maxima failed.
%       localMaxima   : Structure array of length = number of frames in
%                       movie, containing the field "cands", which is a
%                       structure array of length = number of local maxima
%                       in each frame, containing the fields:
%             .IBkg       : Mean background intensity around local maximum.
%             .Lmax       : Position of local maximum (in matrix coordinates).
%             .amp        : Amplitude of local maximum.
%             .pValue     : P-value of local maximum in statistical test
%                           determining its significance.
%       background    : Structure with fields:
%             .meanRawLast5: Mean background intensity in raw movie as
%                            calculated from the last 5 frames.
%             .stdRawLast5 : Standard deviation of background intensity in
%                            raw movie as calculated from the 5 frames.
%             .meanIntegFLast1: Mean background intensity in last frame of
%                               integrated movie.
%             .stdIntegFLast1 : Standard deviation of background intensity
%                               in last frame of integrated movie.
%             .meanIntegFFirst1: Mean background intensity in first frame of
%                                integrated movie.
%             .stdIntegFFirst1 : Standard deviation of background intensity
%                                in first frame of integrated movie.
%       psfSigma      : Standard deviation of point spread function as
%                       estimated from fitting to local maxima in the movie.
%       signal2noiseRatio: Number of features - by - number of frames
%                       array showing signal to noise ratio of all
%                       features in all frames (SNR = signal amplitude
%                       above background / local background std). - WILL
%                       IMPLEMENT SOON.
%       errFlag       : 0 if function executes normally, 1 otherwise.
%
%Khuloud Jaqaman, September 2007
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

%% Output

movieInfo = [];
exceptions = [];
localMaxima = [];
background = [];
psfSigma = [];

%% Input + initialization

%check whether correct number of input arguments was used
if nargin < 2
    disp('--detectSubResFeatures2D_StandAlone: Incorrect number of input arguments!');
    return
end

%get movie parameters
hasImageDir = isfield(movieParam,'imageDir');
if hasImageDir
    imageDir = movieParam.imageDir;
    filenameBase = movieParam.filenameBase;
    digits4Enum = movieParam.digits4Enum;
else
    channel = movieParam.channel;
end
firstImageNum = movieParam.firstImageNum;
lastImageNum = movieParam.lastImageNum;

%get initial guess of PSF sigma
psfSigma = detectionParam.psfSigma;

%get position calculation method
if ~isfield(detectionParam,'calcMethod')
    calcMethod = 'g';
else
    calcMethod = detectionParam.calcMethod;
end

%get statistical test alpha values
if ~isfield(detectionParam,'testAlpha') || isempty(detectionParam.testAlpha)
    testAlpha = struct('alphaR',0.05,'alphaA',0.05,'alphaD',0.05,'alphaF',0);
else
    testAlpha = detectionParam.testAlpha;
end

%get visualization option
if ~isfield(detectionParam,'visual') || isempty(detectionParam.visual)
    visual = 0;
else
    visual = detectionParam.visual;
end

%check whether to do MMF
if ~isfield(detectionParam,'doMMF') || isempty(detectionParam.doMMF)
    doMMF = 1;
else
    doMMF = detectionParam.doMMF;
end

%get camera bit depth
if ~isfield(detectionParam,'bitDepth') || isempty(detectionParam.bitDepth)
    bitDepth = 16;
else
    bitDepth = detectionParam.bitDepth;
end

%get alpha-value for local maxima detection
if ~isfield(detectionParam,'alphaLocMax') || isempty(detectionParam.alphaLocMax)
    alphaLocMax = 0.05;
else
    alphaLocMax = detectionParam.alphaLocMax;
end
numAlphaLocMax = length(alphaLocMax);

%check whether to estimate PSF sigma from the data
if ~isfield(detectionParam,'numSigmaIter') || isempty(detectionParam.numSigmaIter)
    numSigmaIter = 10;
else
    numSigmaIter = detectionParam.numSigmaIter;
end

%get integration time window
if ~isfield(detectionParam,'integWindow')
    integWindow = 0;
else
    integWindow = detectionParam.integWindow;
end
numIntegWindow = length(integWindow);

%make sure that alphaLocMax is the same size as integWindow
if numIntegWindow > numAlphaLocMax
    alphaLocMax = [alphaLocMax ...
        alphaLocMax(1)*ones(1,numIntegWindow-numAlphaLocMax)];
end

%get background information if supplied
if ~isfield(detectionParam,'background') || isempty(detectionParam.background)
    absBG = 0;
else
    absBG = 1;
    bgImageDir = detectionParam.background.imageDir;
    bgImageBase = detectionParam.background.filenameBase;
    alphaLocMaxAbs = detectionParam.background.alphaLocMaxAbs;
end

%get mask information if any
if isfield(detectionParam,'roiMask') && ~isempty(detectionParam.roiMask)
    %get ROImask directly
    maskFlag = 1;
    maskImage = detectionParam.roiMask;
    if ndims(maskImage) == 3
        % For now assume rois remains constant for entire time series.
        maskImage = maskImage(:,:,1);
    end
elseif ~isfield(detectionParam,'maskLoc') || isempty(detectionParam.maskLoc)
    maskFlag = 0;
else 
    % load binary mark from path
    maskFlag = 1;
    maskImage = double(imread(detectionParam.maskLoc));
end

%determine where to save results
if nargin < 3 || isempty(saveResults) %if nothing was input
    saveResDir = pwd;
    saveResFile = 'detectedFeatures.mat';
    saveResults.dir = pwd;
else
    if isstruct(saveResults)
        if ~isfield(saveResults,'dir') || isempty(saveResults.dir)
            saveResDir = pwd;
        else
            saveResDir = saveResults.dir;
        end
        if ~isfield(saveResults,'filename') || isempty(saveResults.filename)
            saveResFile = 'detectedFeatures.mat';
        else
            saveResFile = saveResults.filename;
        end
    else
        saveResults = 0;
    end
end

%check verbose state
if nargin < 4 || isempty(verbose)
    verbose = 1;
end

if hasImageDir
    %store the string version of the numerical index of each image
    enumString = repmat('0',lastImageNum,digits4Enum);
    formatString = ['%0' num2str(digits4Enum) 'i'];
    for i=1:lastImageNum
        enumString(i,:) = num2str(i,formatString);
    end
end

%initialize some variables
emptyFrames = [];
framesFailedLocMax = [];
framesFailedMMF = [];

%turn warnings off
warningState = warning('off','all');

%% General image information

%get image indices and number of images
imageIndx = firstImageNum : lastImageNum;
numImagesRaw = lastImageNum - firstImageNum + 1; %raw images
numImagesInteg = repmat(numImagesRaw,1,numIntegWindow) - 2 * integWindow; %integrated images

if hasImageDir
    %read first image and get image size
    if exist([imageDir filenameBase enumString(imageIndx(1),:) '.tif'],'file')
        imageTmp = imread([imageDir filenameBase enumString(imageIndx(1),:) '.tif']);
    else
        disp('First image does not exist! Exiting ...');
        return
    end
    [imageSizeX,imageSizeY] = size(imageTmp);
    clear imageTmp
else
    imageSizeX=channel.owner_.imSize_(1);
    imageSizeY=channel.owner_.imSize_(2);
end

%make mask of one if no mask is supplied
if ~maskFlag
    maskImage = ones(imageSizeX,imageSizeY);
end

%check which images exist and which don't
imageExists = true(numImagesRaw,1);
if hasImageDir
    for iImage = 1 : numImagesRaw
        if ~exist([imageDir filenameBase enumString(imageIndx(iImage),:) '.tif'],'file')
            imageExists(iImage) = 0;
        end
    end
end

%calculate background properties at movie end
last5start = max(numImagesRaw-4,1);
i = 0;
imageLast5 = NaN(imageSizeX,imageSizeY,5);
for iImage = last5start : numImagesRaw
    i = i + 1;
    if imageExists(iImage) && hasImageDir
        imageLast5(:,:,i) = double(imread([imageDir filenameBase ...
            enumString(imageIndx(iImage),:) '.tif']));
    elseif ~hasImageDir
        imageLast5(:,:,i) = double(channel.loadImage(imageIndx(iImage)));
    end
end
%apply mask
imageLast5 = imageLast5 .* repmat(maskImage,[1 1 size(imageLast5,3)]);

imageLast5 = double(imageLast5) / (2^bitDepth-1);
imageLast5(imageLast5==0) = NaN;
[bgMeanRaw,bgStdRaw] = spatialMovAveBG(imageLast5,imageSizeX,imageSizeY);

%get size of absolute background images if supplied
if absBG
    imageTmp = imread([bgImageDir bgImageBase enumString(imageIndx(1),:) '.tif']);
    [bgSizeX,bgSizeY] = size(imageTmp);
end

%% Local maxima detection

%initialize output structure
localMaxima = repmat(struct('cands',[]),numImagesRaw,1);

for iWindow = 1 : numIntegWindow
    
    %initialize progress text
    if verbose
        progressText(0,['Detecting local maxima with integration window = ' num2str(integWindow(iWindow))]);
    end
    
    for iImage = 1 : numImagesInteg(iWindow)
        
        %store raw images in array
        imageRaw = NaN(imageSizeX,imageSizeY,1+2*integWindow(iWindow));
        for jImage = 1 : 1 + 2*integWindow(iWindow)
            if hasImageDir && imageExists(jImage+iImage-1)
                imageRaw(:,:,jImage) = double(imread([imageDir filenameBase ...
                    enumString(imageIndx(jImage+iImage-1),:) '.tif']));
            elseif ~hasImageDir
                imageRaw(:,:,jImage) = double(channel.loadImage(imageIndx(jImage+iImage-1)));
            end
        end
        %apply mask
        imageRaw = imageRaw .* repmat(maskImage,[1 1 1+2*integWindow(iWindow)]);
        
        %replace zeros with NaNs
        %zeros result from cropping that leads to curved boundaries
        imageRaw(imageRaw==0) = NaN;
        
        %normalize images
        imageRaw = imageRaw / (2^bitDepth-1);
        
        %integrate images
        imageInteg = nanmean(imageRaw,3);
        
        %filter integrated image
        %         imageIntegF = filterGauss2D(imageInteg,min(1,psfSigma));
        imageIntegF = filterGauss2D(imageInteg,psfSigma);

        %use robustMean to get mean and std of background intensities
        %in this method, the intensities of actual features will look like
        %outliers, so we are effectively getting the mean and std of the background
        %account for possible spatial heterogeneity by taking a spatial moving
        %average
        
        %get integrated image background noise statistics
        [bgMeanInteg,bgStdInteg] = ...
            spatialMovAveBG(imageInteg,imageSizeX,imageSizeY);
        
        %calculate absolute background mean and std if supplied
        if absBG
            bgRaw = NaN(bgSizeX,bgSizeY,1+2*integWindow(iWindow));
            for jImage = 1 : 1 + 2*integWindow(iWindow)
                if imageExists(jImage+iImage-1)
                    bgRaw(:,:,jImage) = double(imread([bgImageDir bgImageBase ...
                        enumString(imageIndx(jImage+iImage-1),:) '.tif']));
                end
            end
            bgRaw(bgRaw==0) = NaN;
            bgRaw = bgRaw / (2^bitDepth-1);
            bgInteg = nanmean(bgRaw,3);
            bgAbsMeanInteg = nanmean(bgInteg(:))*ones(imageSizeX,imageSizeY);
            bgAbsStdInteg = nanstd(bgInteg(:))*ones(imageSizeX,imageSizeY);
        end
        
        background = [];
        
        %clear some variables
        clear imageInteg
        
        try
            
            %call locmax2d to get local maxima in filtered image
            fImg = locmax2d(imageIntegF,[3 3],1);
            
            %get positions and amplitudes of local maxima
            [localMaxPosX,localMaxPosY,localMaxAmp] = find(fImg);
            localMax1DIndx = find(fImg(:));
            
            %get background values corresponding to local maxima
            bgMeanInteg1 = bgMeanInteg;
            bgMeanMaxF = bgMeanInteg1(localMax1DIndx);
            bgStdInteg1 = bgStdInteg;
            bgStdMaxF = bgStdInteg1(localMax1DIndx);
            bgMeanMax = bgMeanRaw(localMax1DIndx);
            
            %calculate the p-value corresponding to the local maxima's amplitudes
            %assume that background intensity in integrated image is normally
            %distributed with mean bgMeanMaxF and standard deviation bgStdMaxF
            pValue = 1 - normcdf(localMaxAmp,bgMeanMaxF,bgStdMaxF);
            
            if absBG
                bgAbsMeanMaxF = bgAbsMeanInteg(localMax1DIndx);
                bgAbsStdMaxF = bgAbsStdInteg(localMax1DIndx);
                pValueAbs = 1 - normcdf(localMaxAmp,bgAbsMeanMaxF,bgAbsStdMaxF);
            end
            
            %retain only those maxima with significant amplitude
            if absBG
                keepMax = find((pValue < alphaLocMax(iWindow)) & ...
                    (pValueAbs < alphaLocMaxAbs));
            else
                keepMax = find(pValue < alphaLocMax(iWindow));
            end
            localMaxPosX = localMaxPosX(keepMax);
            localMaxPosY = localMaxPosY(keepMax);
            localMaxAmp = localMaxAmp(keepMax);
            bgMeanMax = bgMeanMax(keepMax);
            pValue = pValue(keepMax);
            numLocalMax = length(keepMax);
            
            %construct cands structure
            if numLocalMax == 0 %if there are no local maxima
                
                cands = [];
                %                 emptyFrames = [emptyFrames; iImage+integWindow]; %#ok<AGROW>
                
            else %if there are local maxima
                
                %define background mean and status
                cands = repmat(struct('status',1,'IBkg',[],...
                    'Lmax',[],'amp',[],'pValue',[]),numLocalMax,1);
                
                %store maxima positions, amplitudes and p-values
                for iMax = 1 : numLocalMax
                    cands(iMax).IBkg = bgMeanMax(iMax);
                    cands(iMax).Lmax = [localMaxPosX(iMax) localMaxPosY(iMax)];
                    cands(iMax).amp = localMaxAmp(iMax);
                    cands(iMax).pValue = pValue(iMax);
                end
                
            end
            
            %add the cands of the current image to the rest - this is done
            %for the raw images, not the integrated ones
            localMaxima(iImage+integWindow(iWindow)).cands = ...
                [localMaxima(iImage+integWindow(iWindow)).cands; cands];
            
        catch %#ok<CTCH>
            
            %             %if local maxima detection fails, make cands empty
            %             localMaxima(iImage+integWindow).cands = [];
            %
            %             %add this frame to the array of frames with failed local maxima
            %             %detection and to the array of empty frames
            %             framesFailedLocMax = [framesFailedLocMax; iImage+integWindow]; %#ok<AGROW>
            %             emptyFrames = [emptyFrames; iImage+integWindow]; %#ok<AGROW>
            
        end
        
        %display progress
        if verbose
            progressText(iImage/numImagesInteg(iWindow),['Detecting local maxima with integration window = ' num2str(integWindow(iWindow))]);
        end
        
    end %(for iImage = 1 : numImagesInteg(iWindow))
    
    %assign local maxima for frames left out due to time integration
    for iImage = 1 : integWindow(iWindow)
        localMaxima(iImage).cands = [localMaxima(iImage).cands; ...
            localMaxima(integWindow(iWindow)+1).cands];
    end
    for iImage = numImagesRaw-integWindow(iWindow)+1 : numImagesRaw
        localMaxima(iImage).cands = [localMaxima(iImage).cands; ...
            localMaxima(end-integWindow(iWindow)).cands];
    end
    
    % if any(emptyFrames==integWindow(iWindow)+1)
    %     emptyFrames = [emptyFrames; (1:integWindow)'];
    % end
    % if any(emptyFrames==numImagesRaw-integWindow)
    %     emptyFrames = [emptyFrames; (numImagesRaw-integWindow+1:numImagesRaw)'];
    % end
    % if any(framesFailedLocMax==integWindow+1)
    %     framesFailedLocMax = [framesFailedLocMax; (1:integWindow)'];
    % end
    % if any(framesFailedLocMax==numImagesRaw-integWindow)
    %     framesFailedLocMax = [framesFailedLocMax; (numImagesRaw-integWindow+1:numImagesRaw)'];
    % end
    
end %(for iWindow = 1 : numIntegWindow)

%delete whatever local maxima were found in the frames that don't exist,
%because they are clearly an artifact of time-integration
for iFrame = (find(imageExists==0))'
    localMaxima(iFrame).cands = [];
end

%go over all frames, remove redundant cands, and register empty frames
if verbose
    progressText(0,'Removing redundant local maxima');
end
for iImage = 1 : numImagesRaw
    
    %get the cands of this frame
    candsCurrent = localMaxima(iImage).cands;
    
    %if there are no cands, register that this is an empty frame
    if isempty(candsCurrent)
        
        emptyFrames = [emptyFrames; iImage]; %#ok<AGROW>
        
    else
        
        %get the local maxima positions in this frame
        maxPos = vertcat(candsCurrent.Lmax);
        
        %find the unique local maxima positions
        [~,indxUnique] = unique(maxPos,'rows');
        
        %keep only these unique cands
        candsCurrent = candsCurrent(indxUnique);
%         maxPos = vertcat(candsCurrent.Lmax);
%         
%         %if there is more than one surviving cand
%         if size(maxPos,1) > 1
%             
%             %remove cands that are closer than 2*psfSigma to each other ...
%             
%             %first do that by clustering the cands ...
%             
%             %calculate the distances between cands
%             y = pdist(maxPos);
%             
%             %get the linkage between cands using maximum distance
%             Z = linkage(y,'complete');
%             
%             %cluster the cands and keep only 1 cand from each cluster
%             T = cluster(Z,'cutoff',2*psfSigma,'criterion','distance');
%             [~,cands2keep] = unique(T);
%             
%             %update list of cands
%             candsCurrent = candsCurrent(cands2keep);
%             maxPos = vertcat(candsCurrent.Lmax);
%             
%             if size(maxPos,1) > 1
%                 
%                 %then refine that by removing cands one by one ...
%                 
%                 %calculate the distances between surviving cands
%                 distBetweenCands = createDistanceMatrix(maxPos,maxPos);
%                 
%                 %find the minimum distance for each cand
%                 distBetweenCandsSort = sort(distBetweenCands,2);
%                 distBetweenCandsSort = distBetweenCandsSort(:,2:end);
%                 minDistBetweenCands = distBetweenCandsSort(:,1);
%                 
%                 %find the minimum minimum distance
%                 minMinDistBetweenCands = min(minDistBetweenCands);
%                 
%                 %if this distance is smaller than 2*psfSigma, remove the
%                 % maximum with smallest average distance to its neighbors
%                 while minMinDistBetweenCands <= (2 * psfSigma)
%                     
%                     %find the cands involved
%                     candsInvolved = find(distBetweenCandsSort(:,1) == minMinDistBetweenCands);
%                     
%                     %determine which one of them has the smallest average distance
%                     %to the other cands
%                     aveDistCands = mean(distBetweenCandsSort(candsInvolved,:),2);
%                     cand2remove = candsInvolved(aveDistCands==min(aveDistCands));
%                     cands2keep = setdiff((1:size(maxPos,1))',cand2remove(1));
%                     
%                     %remove it from the list of cands
%                     candsCurrent = candsCurrent(cands2keep);
%                     maxPos = vertcat(candsCurrent.Lmax);
%                     
%                     %repeat the minimum distance calculation
%                     if size(maxPos,1) > 1
%                         distBetweenCands = createDistanceMatrix(maxPos,maxPos);
%                         distBetweenCandsSort = sort(distBetweenCands,2);
%                         distBetweenCandsSort = distBetweenCandsSort(:,2:end);
%                         minDistBetweenCands = distBetweenCandsSort(:,1);
%                         minMinDistBetweenCands = min(minDistBetweenCands);
%                     else
%                         minMinDistBetweenCands = 3 * psfSigma;
%                     end
%                     
%                 end %(while minMinDistBetweenCands <= (2 * psfSigma))
%                 
%             end %(if size(maxPos,1) > 1)
%             
%         end %(if size(maxPos,1) > 1)
        
        localMaxima(iImage).cands = candsCurrent;
        
    end
    
    %display progress
    if verbose
        progressText(iImage/numImagesRaw,'Removing redundant local maxima');
    end
    
end

%make a list of images that have local maxima
goodImages = setxor(1:numImagesRaw,emptyFrames);
numGoodImages = length(goodImages);

%clear some variables
clear ImageIntegF

%% PSF sigma estimation

if numSigmaIter
    
    %specify which parameters to fit for sigma estimation
    fitParameters = [{'X1'} {'X2'} {'A'} {'Sxy'} {'B'}];
    
    %store original input sigma
    psfSigmaIn = psfSigma;
    
    %give a dummy value for psfSigma0 and acceptCalc to start while loop
    psfSigma0 = 0;
    acceptCalc = 1;
    
    %initialize variable counting number of iterations
    numIter = 0;
    
    %iterate as long as estimated sigma is larger than initial sigma
    while numIter <= numSigmaIter && acceptCalc && ((psfSigma-psfSigma0)/psfSigma0 > 0.05)
        
        %add one to number of iterations
        numIter = numIter + 1;
        
        %save input PSF sigma in new variable and empty psfSigma for estimation
        psfSigma0 = psfSigma;
        psfSigma = [];
        
        %calculate some numbers that get repeated many times
        psfSigma5 = ceil(5*psfSigma0);
        
        %initialize progress display
        if verbose
            switch numIter
                case 1
                    progressText(0,'Estimating PSF sigma');
                otherwise
                    progressText(0,'Repeating PSF sigma estimation');
            end
        end
        
        %go over the first 50 good images and find isolated features
        images2use = goodImages(1:min(50,numGoodImages));
        images2use = setdiff(images2use,1:integWindow);
        for iImage = images2use(:)'
            
            %read raw image
            if hasImageDir
                imageRaw = imread([imageDir filenameBase enumString(imageIndx(iImage),:) '.tif']);                
            else
                imageRaw= channel.loadImage(imageIndx(iImage));
            end
            imageRaw = double(imageRaw) / (2^bitDepth-1);
            %apply mask
            imageRaw = imageRaw .* maskImage;
            
            %get feature positions and amplitudes and average background
            featPos = vertcat(localMaxima(iImage).cands.Lmax);
            featAmp = vertcat(localMaxima(iImage).cands.amp);
            featBG  = vertcat(localMaxima(iImage).cands.IBkg);
            featPV  = vertcat(localMaxima(iImage).cands.pValue);
            
            %retain only features that are more than 5*psfSigma0 away from boundaries
            feat2use = find(featPos(:,1) > psfSigma5 & ...
                featPos(:,1) < imageSizeX - psfSigma5 & ...
                featPos(:,2) > psfSigma5 & featPos(:,2) < imageSizeY - psfSigma5);
            featPos = featPos(feat2use,:);
            featAmp = featAmp(feat2use);
            featBG = featBG(feat2use);
            featPV = featPV(feat2use);
            
            %if there is more than one feature ...
            if length(feat2use) > 1
                
                %find nearest neighbor distances
                nnDist = createDistanceMatrix(featPos,featPos);
                nnDist = sort(nnDist,2);
                nnDist = nnDist(:,2);
                
                %retain only features whose nearest neighbor is more than 10*psfSigma0
                %away
                feat2use = find(nnDist > ceil(10*psfSigma0));
                featPos = featPos(feat2use,:);
                featAmp = featAmp(feat2use);
                featBG = featBG(feat2use);
                featPV = featPV(feat2use);
                
                %retain only features with pValue between the 25th and 75th
                %percentiles
                percentile25 = prctile(featPV,25);
                percentile75 = prctile(featPV,75);
                feat2use = find(featPV > percentile25 & featPV < percentile75);
                featPos = featPos(feat2use,:);
                featAmp = featAmp(feat2use);
                featBG = featBG(feat2use);
                
            end
            
            %go over the selected features and estimate psfSigma
            numFeats = length(featAmp);
            parameters = zeros(numFeats,5);
            if numFeats >= 1
                
                for iFeat = 1 : numFeats
                    
                    %crop image around selected feature
                    lowerBound = featPos(iFeat,:) - psfSigma5;
                    upperBound = featPos(iFeat,:) + psfSigma5;
                    imageCropped = imageRaw(lowerBound(1):upperBound(1),...
                        lowerBound(2):upperBound(2),1);
                    
                    %estimate sigma if image region contains no NaNs
                    %NaNs appear due to cropping
                    if all(~isnan(imageCropped(:)))
                        
                        %make initial guess for fit (in the order given in fitParameters)
                        initGuess = [psfSigma5+1 psfSigma5+1 featAmp(iFeat) ...
                            psfSigma0 featBG(iFeat)];
                        
                        %fit image and estimate sigma of Gaussian
                        parameters(iFeat,:) = GaussFitND(imageCropped,[],...
                            fitParameters,initGuess);
                        
                        %                         %just to test whether GaussFitND works either way
                        %                         %the answer is yes
                        %                         [imSizeY,imSizeX] = size(imageCropped);
                        %                         coord4fitX = repmat((1:imSizeX)',imSizeY,1);
                        %                         coord4fitY = repmat((1:imSizeY),imSizeX,1);
                        %                         coord4fitY = coord4fitY(:);
                        %                         parameters2(iFeat,:) = GaussFitND(imageCropped(:),[coord4fitX coord4fitY],...
                        %                             fitParameters,initGuess);
                        
                    else %otherwise assign NaN
                        
                        parameters(iFeat,:) = NaN;
                        
                    end
                    
                end
                
                %add to array of sigmas
                psfSigma = [psfSigma; parameters(:,4)]; %#ok<AGROW>
                
            end %(if numFeats >= 1)
            
            %display progress
            if verbose
                switch numIter
                    case 1
                        progressText(iImage/max(images2use),'Estimating PSF sigma');
                    otherwise
                        progressText(iImage/max(images2use),'Repeating PSF sigma estimation');
                end
            end
            
        end %(for iImage = images2use)
        
        %estimate psfSigma as the robust mean of all the sigmas from the fits
        psfSigma = psfSigma(~isnan(psfSigma)); %get rid of NaNs from cropped regions
        numCalcs = length(psfSigma);
        if numCalcs > 0
            
            [psfSigma,~,inlierIndx] = robustMean(psfSigma,[],3,0,true);
            
            numInlierIndx = sum(inlierIndx);
            
            %accept new sigma if there are enough observations and inliers
            acceptCalc = (numCalcs >= 100 && numInlierIndx >= 0.7*numCalcs) || ...
                (numCalcs >= 50 && numInlierIndx >= 0.9*numCalcs) || ...
                (numCalcs >= 10 && numInlierIndx == numCalcs);
            
        else
            
            acceptCalc = 0;
            
        end
        
        %show new sigma if estimation is accepted
        if acceptCalc
            fprintf('PSF sigma = %1.3f (%d inliers out of %d observations)',...
                psfSigma,numInlierIndx,numCalcs);
        else %otherwise alert user that input sigma was retained
            psfSigma = psfSigmaIn;
            disp('Not enough observations to change PSF sigma, using input PSF sigma');
        end
        
    end %(while numIter <= numSigmaIter && acceptCalc && ((psfSigma-psfSigma0)/psfSigma0 > 0.05))
    
    %if maximum number of iterations has been performed but sigma value is not converging
    if numIter == numSigmaIter+1 && acceptCalc && ((psfSigma-psfSigma0)/psfSigma0 > 0.05)
        psfSigma = psfSigmaIn;
        disp('Estimation terminated (no convergence), using input PSF sigma');
    end
    
end %(if numSigmaIter)

%% Mixture-model fitting

%initialize movieInfo
clear movieInfo
movieInfo = repmat(struct('xCoord',[],'yCoord',[],'amp',[]),numImagesRaw,1);

%initialize progress display
if verbose
    if strcmp(calcMethod,'g')
        progressText(0,'Mixture-model fitting');
    else
        progressText(0,'Centroid calculation');
    end
end

%go over all non-empty images ...
for iImage = goodImages(:)'
    
    %read raw image
    if hasImageDir
        imageRaw = imread([imageDir filenameBase enumString(imageIndx(iImage),:) '.tif']);
    else
        imageRaw = channel.loadImage(imageIndx(iImage));
    end
    imageRaw = double(imageRaw) / (2^bitDepth-1);
    %apply mask
    imageRaw = imageRaw .* maskImage;
    
    try %try to detect features in this frame
        
        %fit with mixture-models
        if strcmp(calcMethod,'g')
            featuresInfo = detectSubResFeatures2D_V2(imageRaw,...
                localMaxima(iImage).cands,psfSigma,testAlpha,visual,...
                doMMF,1,0,mean(bgStdRaw(:)));
        else
            featuresInfo = centroidSubResFeatures2D(imageRaw,...
                localMaxima(iImage).cands,psfSigma,visual,1,0);
        end
        
        %save results
        movieInfo(iImage) = featuresInfo;
        
        %check whether frame is empty
        if isempty(featuresInfo.xCoord)
            emptyFrames = [emptyFrames; iImage]; %#ok<AGROW>
        end
        
    catch %#ok<CTCH> %if detection fails
        
        %label frame as empty
        emptyFrames = [emptyFrames; iImage]; %#ok<AGROW>
        
        %add this frame to the array of frames with failed mixture-model
        %fitting
        framesFailedMMF = [framesFailedMMF; iImage]; %#ok<AGROW>
        
    end
    
    %display progress
    if verbose
        if strcmp(calcMethod,'g')
            progressText(iImage/numImagesRaw,'Mixture-model fitting');
        else
            progressText(iImage/numImagesRaw,'Centroid calculation');
        end
    end
    
end

%% Post-processing

%sort list of empty frames
emptyFrames = sort(emptyFrames);

%store empty frames and frames where detection failed in structure
%exceptions
exceptions = struct('emptyFrames',emptyFrames,'framesFailedLocMax',...
    framesFailedLocMax,'framesFailedMMF',framesFailedMMF');

%indicate correct frames in movieInfo
tmptmp = movieInfo;
clear movieInfo
movieInfo(firstImageNum:lastImageNum,1) = tmptmp;

%save results
if isstruct(saveResults)
    save([saveResDir filesep saveResFile],'movieParam','detectionParam',...
        'movieInfo','exceptions','localMaxima','background','psfSigma');
end

%go back to original warnings state
warning(warningState);


%% ~~~ the end ~~~

