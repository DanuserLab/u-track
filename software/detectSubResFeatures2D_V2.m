function [detectedFeatures,clustersMMF,imageN3,errFlag] = ...
    detectSubResFeatures2D_V2(image,cands,psfSigma,testAlpha,visual,...
    doMMF,bitDepth,saveResults,bgNoiseSigma)
%DETECTSUBRESFEATURES2D_V2 determines the positions and intensity amplitudes of sub-resolution features using mixture model fitting
%
%SYNOPSIS [detectedFeatures,clustersMMF,imageN3,errFlag] = ...
%    detectSubResFeatures2D_V2(image,cands,psfSigma,testAlpha,visual,...
%    doMMF,bitDepth,saveResults,bgNoiseSigma)
%
%INPUT  image      : Image being analyzed.
%       cands      : Cands structure as output from fsmCenter.
%       psfSigma   : Standard deviation of point spread function (in pixels).
%       testAlpha  : Alpha-values for statistical tests. Structure with fields:
%             .alphaR: For the residuals test, comparing N+1-kernel fit to
%                      N-kernal fit. Optional. Default: 0.05.
%             .alphaA: For amplitude test. Optional. Default: 0.05.
%             .alphaD: For distance test. Optional. Default: 0.05.
%             .alphaF: Final residuals test, comparing residuals from final
%                      fit to estimated background noise.
%                      Optional. Default: 0.
%       visual     : 1 if user wants to view results; 0 otherwise.
%                    Optional. Default: 0.
%       doMMF      : 1 if user wants to do mixture-model fitting, 0
%                    otherwise. Optional. Default: 1.
%       bitDepth   : Camera bit depth. Optional. Default: 14.
%       saveResults: 1 if results are to be saved (in file 'detectedFeatures.mat'),
%                    0 otherwise. Optional. Default: 0.
%       bgNoiseSigma:Standard deviation of background noise. Optional. If
%                    not input, the code will estimate it from the image.
%
%       All optional variables can be entered as [] to use default values.
%
%OUTPUT detectedFeatures: Structure with fields:
%             .xCoord    : Image coordinate system x-coordinate of detected
%                          features [x dx] (in pixels).
%             .yCoord    : Image coorsinate system y-coordinate of detected
%                          features [y dy] (in pixels).
%             .amp       : Amplitudes of PSFs fitting detected features [a da].
%       clustersMMF: Array of clusters of sub-resolution features.
%                    Structure with fields:
%             .position  : Position of each feature in image coordinate
%                          system (in pixels): [x y dx dy] = [y x dy dx]
%                          in matrix coordinate system.
%             .amplitude : Intensity of each feature [A dA].
%             .bgAmp     : Background intensity [Abg dAbg].
%             .varCovMat : Variance-covariance matrix of estimated parameters.
%       imageN3    : Image with labeled features. Blue: those from cands;
%                    Red: those from mixture-model fitting; Magenta: those
%                    from MMF which coincide with those from cands.
%                    Will be output only if visual = 1.
%       errFlag    : 0 if function executes normally, 1 otherwise.
%
%Khuloud Jaqaman, August 2005
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

detectedFeatures = [];
clustersMMF = [];
imageN3 = [];
errFlag = 0;

%% Input

%check whether correct number of input arguments was used
if nargin < 3
    disp('--detectSubResFeatures2D_V2: Incorrect number of input arguments!');
    errFlag  = 1;
    return
end

%check whether statistical test alpha values were inputted
if nargin < 4 || isempty(testAlpha) %if not, assign defaults
    
    testAlpha = struct('alphaR',0.05,'alphaA',0.05,'alphaD',0.05,'alphaF',0.05);
    
else %if some were, check their values and assign default for the rest
    
    if ~isfield(testAlpha,'alphaR')
        testAlpha.alphaR = 0.05;
    else
        if testAlpha.alphaR < 0 || testAlpha.alphaR > 1
            disp('--detectSubResFeatures2D_V2: testAlpha.alphaR should be between 0 and 1!');
            errFlag = 1;
        end
    end
    if ~isfield(testAlpha,'alphaA')
        testAlpha.alphaA = 0.05;
    else
        if testAlpha.alphaA < 0 || testAlpha.alphaA > 1
            disp('--detectSubResFeatures2D_V2: testAlpha.alphaA should be between 0 and 1!');
            errFlag = 1;
        end
    end
    if ~isfield(testAlpha,'alphaD')
        testAlpha.alphaD = 0.05;
    else
        if testAlpha.alphaD < 0 || testAlpha.alphaD > 1
            disp('--detectSubResFeatures2D_V2: testAlpha.alphaD should be between 0 and 1!');
            errFlag = 1;
        end
    end
    if ~isfield(testAlpha,'alphaF')
        testAlpha.alphaF = 0;
    else
        if testAlpha.alphaF < 0 || testAlpha.alphaF > 1
            disp('--detectSubResFeatures2D_V2: testAlpha.alphaF should be between 0 and 1!');
            errFlag = 1;
        end
    end
    
end

%check visualization option
if nargin < 5 || isempty(visual)
    visual = 0;
else
    if visual ~= 0 && visual ~= 1
        disp('--detectSubResFeatures2D_V2: Variable "visual" should be 0 or 1!');
        errFlag = 1;
    end
end

%check whether to do MMF
if nargin < 6 || isempty(doMMF)
    doMMF = 1;
else
    if doMMF ~= 0 && doMMF ~= 1
        disp('--detectSubResFeatures2D_V2: Variable "doMMF" should be 0 or 1!');
        errFlag = 1;
    end
end

%check the bit depth
if nargin < 7 || isempty(bitDepth)
    bitDepth = 14;
else
    if bitDepth <= 0 || bitDepth-floor(bitDepth) ~= 0
        disp('--detectSubResFeatures2D_V2: Variable "bitDepth" should be a positive integer!');
    end
end

%check whether results are to be saved
if nargin < 8 || isempty(saveResults)
    saveResults = 0;
else
    if saveResults ~= 0 && saveResults ~= 1
        disp('--detectSubResFeatures2D_V2: Variable "saveResults" should be 0 or 1!');
    end
end

%check whether back ground noise sigma is input or whether it should be
%calculated on the fly
if nargin < 9 || isempty(bgNoiseSigma)
    bgNoiseSigma = 0;
    estimateBgNoise = 1;
else
    estimateBgNoise = 0;
end

%exit if there are problems with input data
if errFlag
    disp('--detectSubResFeatures2D_V2: Please fix input data!');
    return
end

%get number of pixels in each direction (in image coordinate system)
[numPixelsY,numPixelsX] = size(image);

%extract test alpha values from input
alphaR = testAlpha.alphaR;
alphaA = testAlpha.alphaA;
alphaD = testAlpha.alphaD;
alphaF = testAlpha.alphaF;

%Divide image by bit depth, to normalize it between 0 and 1
image = double(image)/(2^bitDepth-1);

%get background intensity information from cands
bgAmp = vertcat(cands.IBkg);
status = vertcat(cands.status);
bgAmp = bgAmp(status==1);
bgAmpMax = max(bgAmp);
bgAmpAve = mean(bgAmp);

%% Determine signal overlap

%determine which signals are overlapping, in which case they must
%be fitted together later on
[clusters,errFlag] = findOverlapPSFs2D(cands,numPixelsX,numPixelsY,psfSigma);
if errFlag
    disp('--detectSubResFeatures2D_V2: Could not place signals in clusters!');
    return
end

%% Mixture Model Fitting

%initialize vector indicating whether clusters should be retained
keepCluster = ones(length(clusters),1);

%set optimization options
options = optimset('Jacobian','on',...
    'MaxFunEvals', 1e4, ...
    'MaxIter', 1e4, ...
    'Display', 'off', ...
    'TolX', 1e-6, ...
    'Tolfun', 1e-6);

%reserve memory for clustersMMF
numClusters = length(clusters);
clustersMMF = repmat(struct('position',[],'amplitude',[],'bgAmp',[],...
    'numDegFree',[],'residuals',[]),numClusters,1);

%go over all clusters
for i = 1 : numClusters
    
    %get initial guess of positions and amplitudes
    numMaximaT = clusters(i).numMaxima;
    maximaPosT = clusters(i).maximaPos(:,1:2);
    maximaAmpT = clusters(i).maximaAmp;
    bgAmpT = bgAmpAve;
    clusterPixels = clusters(i).pixels(:,1:2);
    
    %% remove superfluous maxima
    if numMaximaT > 1
        
        %calculate the distances between all features
        distBetweenMax = createDistanceMatrix(maximaPosT,maximaPosT);
        
        %find the minimum distance for each maximum
        distBetweenMaxSort = sort(distBetweenMax,2);
        distBetweenMaxSort = distBetweenMaxSort(:,2:end);
        minDistBetweenMax = distBetweenMaxSort(:,1);
        
        %find the minimum minimum distance
        minMinDistBetweenMax = min(minDistBetweenMax);
        
        %if this distance is smaller than 2*psfSigma, remove the maximum with
        %smallest average distance to its neighbors
        while minMinDistBetweenMax <= (2 * psfSigma)
            
            %find the two maxima involved
            ijMax = find(distBetweenMaxSort(:,1) == minMinDistBetweenMax);
            
            %determine which one of them has the smaller average distance to
            %the other maxima
            aveDistIJ = mean(distBetweenMaxSort(ijMax,:),2);
            max2remove = ijMax(aveDistIJ==min(aveDistIJ));
            max2keep = setdiff((1:numMaximaT)',max2remove(1));
            
            %remove it from the cluster
            numMaximaT = numMaximaT - 1;
            maximaPosT = maximaPosT(max2keep,:);
            maximaAmpT = maximaAmpT(max2keep,:);
            
            %repeat the minimum distance calculation
            if numMaximaT > 1
                distBetweenMax = createDistanceMatrix(maximaPosT,maximaPosT);
                distBetweenMaxSort = sort(distBetweenMax,2);
                distBetweenMaxSort = distBetweenMaxSort(:,2:end);
                minDistBetweenMax = distBetweenMaxSort(:,1);
                minMinDistBetweenMax = min(minDistBetweenMax);
            else
                minMinDistBetweenMax = 3 * psfSigma;
            end
            
        end
        
    end
    
    %% attempt first fit and iterative mixture-model fitting if requested

    %crop part of image that is relevant to this fitting
    imageC = image(clusters(i).pixels(:,3));
    
    firstFit = 1; %logical variable indicating first fit attempt
    fit = 1; %logical variable indicating whether to attempt to fit
    
    %fit if flag is 1
    while fit
        
        %collect initial guesses and lower and upper bounds for fit
        [x0,lb,ub] = mmfInitGuessLowerUpperBounds(maximaPosT,maximaAmpT,...
            bgAmpT,psfSigma,clusterPixels,firstFit);
        
        %calculate number of degrees of freedom in system
        numDegFreeT = size(clusterPixels,1)-3*numMaximaT-1;
        
        %determine feature positions and amplitudes and estimate background
        %intensity, using nonlinear least squares data fitting
        %also get residuals and Jacobian matrix to calculate variance-covariance
        %of estimated parameters
        [solutionT,dummy,residualsT,dummy,dummy,dummy,jacMatT] = ...
            lsqnonlin(@fitNGaussians2D,x0,lb,ub,options,imageC,...
            clusterPixels,psfSigma);
        jacMatT = full(jacMatT);
        residualsT = -residualsT; %minus sign so that residuals = real image - model image
        residVarT = (sum(residualsT.^2)/numDegFreeT);
        
        %check whether addition of 1 PSF has significantly improved the fit
        if firstFit %if this is the first fit
            
            firstFit = 0; %next one won't be
            
        else %if this is not the first fit
            
            %get test statistic, which is F-distributed
            testStat = residVarT/residVar;
            
            %get p-value of test statistic
            pValue = fcdf(testStat,numDegFree,numDegFreeT);
            
            %compare p-value to alpha
            %1-sided F-test: H0: F=1, H1: F<1
            if pValue < alphaR %if p-value is smaller, accept this fit
                fit = 1; %and attempt another one with an additional kernel
            else %if p-value is larger, do not accept this fit and exit
                fit = 0;
            end
            
        end
        
        if fit %if this fit is accepted (which is the default if it's the first fit)
            
            %update variables
            numMaxima = numMaximaT;
            numDegFree = numDegFreeT;
            solution = solutionT;
            residuals = residualsT;
            residVar = residVarT;
            jacMat = jacMatT;
            
            %calculate the parameters' variance-covariance matrix and get their
            %uncertainties
            varCovMat = residVar * inv(jacMat'*jacMat);
            standDevVec = sqrt(diag(varCovMat));
            
            %if nothing weird happened in the fit...
            if all(isreal(standDevVec))
                
                %extract estimate and std of background intensity and
                %remove from vectors
                bgAmp = [solution(end) standDevVec(end)];
                solution = solution(1:end-1);
                standDevVec = standDevVec(1:end-1);
                
                %reshape 3nx1 vectors "solution" and "standDevVec" into nx3 matrices
                solution = reshape(solution,3,numMaxima);
                solution = solution';
                standDevVec = reshape(standDevVec,3,numMaxima);
                standDevVec = standDevVec';
                
                %extract feature positions and amplitudes and their uncertainties
                maximaPos = [solution(:,1:2) standDevVec(:,1:2)];
                maximaAmp = [solution(:,3) standDevVec(:,3)];
                
                %if attempting iterative mixture-model fitting
                if doMMF
                    
                    %add a kernel positioned at the pixel with maximum residual
                    %for fit with numMaxima + 1 Gaussians
                    numMaximaT = numMaxima + 1; %update number of maxima
                    maximaAmpT = [maximaAmp(:,1); mean(maximaAmp(:,1))]; %signal amplitude
                    indxMaxRes = find(residuals==max(residuals)); %index of pixel with maximum residuals
                    coord = clusterPixels(indxMaxRes(1),:); %position of new kernel (using indxMaxRes(1) avoids crashing if it so happens that two pixels have the exact same value which happens to be the maximum value - it has happened!!!)
                    maximaPosT = [maximaPos(:,1:2); coord]; %signal positions
                    bgAmpT = bgAmp(1); %background amplitude
                    
                else %if no iterative mixture model fitting
                    
                    fit = 0;
                    
                end
                
            else %if something weird did happen in the fit
                
                fit = 0; %stop fitting this cluster
                keepCluster(i) = 0; %mark it to be discarded
                
            end %(if all(isreal(standDevVec)))
            
        end %(if fit)
        
    end %(while fit)
    
    %if nothing weird happened in the fit and this cluster is to be
    %retained
    if keepCluster(i)
        
        %% check amplitudes and make sure that they are all significant
        
        %1-sided t-test: H0: T=0, H1: T>0
        %calculate test statistic (t-distributed)
        testStat = maximaAmp(:,1)./sqrt((maximaAmp(:,2).^2+residVar));
        
        %get p-value
        pValue = 1-tcdf(testStat,numDegFree);
        
        %find largest p-value and decide whether to remove one kernel,
        %repeat the fit and test again for amplitude
        pValueMax = max(pValue);
        testAmp = pValueMax > alphaA;
        
        %if any of the amplitudes is not significant
        while testAmp
            
            %remove a kernel and refit only if there is originally more
            %than one kernel
            %otherwise, discard the whole cluster
            if numMaxima > 1
                
                %determine which kernels to keep
                indxBad = find(pValue == pValueMax);
                indxBad = indxBad(1);
                indx = setdiff(1:numMaxima,indxBad);
                
                %remove the information of the kernel to be discarded
                maximaPos = maximaPos(indx,:);
                maximaAmp = maximaAmp(indx,:);
                numMaxima = numMaxima - 1;
                
                %collect initial guesses and lower and upper bounds for fit
                [x0,lb,ub] = mmfInitGuessLowerUpperBounds(maximaPos(:,1:2),...
                    maximaAmp(:,1),bgAmp(1),psfSigma,clusterPixels,1);
                
                %calculate number of degrees of freedom in system
                numDegFree = size(clusterPixels,1) - 3*numMaxima - 1;
                
                %determine feature positions and amplitudes and estimate background
                %intensity, using nonlinear least squares data fitting
                %also get residuals and Jacobian matrix to calculate variance-covariance
                %of estimated parameters
                [solution,dummy,residuals,dummy,dummy,dummy,jacMat] = ...
                    lsqnonlin(@fitNGaussians2D,x0,lb,ub,options,imageC,...
                    clusterPixels,psfSigma);
                jacMat = full(jacMat);
                residuals = -residuals; %minus sign so that residuals = real image - model image
                residVar = sum(residuals.^2)/numDegFree;
                
                %calculate the parameters' variance-covariance matrix and get their
                %uncertainties
                varCovMat = residVar * inv(jacMat'*jacMat);
                standDevVec = sqrt(diag(varCovMat));
                
                %if nothing weird happened in the fit...
                if all(isreal(standDevVec))
                    
                    %extract estimate and std of background intensity and
                    %remove from vectors
                    bgAmp = [solution(end) standDevVec(end)];
                    solution = solution(1:end-1);
                    standDevVec = standDevVec(1:end-1);
                    
                    %reshape 3nx1 vectors "solution" and "standDevVec" into nx3 matrices
                    solution = reshape(solution,3,numMaxima);
                    solution = solution';
                    standDevVec = reshape(standDevVec,3,numMaxima);
                    standDevVec = standDevVec';
                    
                    %extract feature positions and amplitudes and their uncertainties
                    maximaPos = [solution(:,1:2) standDevVec(:,1:2)];
                    maximaAmp = [solution(:,3) standDevVec(:,3)];
                    
                    %check amplitudes and make sure that they are all significant
                    
                    %1-sided t-test: H0: T=0, H1: T>0
                    %calculate test statistic (t-distributed)
                    testStat = maximaAmp(:,1)./sqrt((maximaAmp(:,2).^2+residVar));
                    
                    %get p-value
                    pValue = 1-tcdf(testStat,numDegFree);
                    
                    %find largest p-value and decide whether to remove one kernel,
                    %repeat the fit and test again for amplitude
                    pValueMax = max(pValue);
                    testAmp = pValueMax > alphaA;
                    
                else %(if all(isreal(standDevVec)))
                    
                    testAmp = 0; %stop fitting this cluster
                    keepCluster(i) = 0; %mark it to be discarded
                    
                end %(if all(isreal(standDevVec)))
                
            else %(if numMaxima > 1)
                
                testAmp = 0; %stop fitting this cluster
                keepCluster(i) = 0; %mark it to be discarded
                
            end %(if numMaxima > 1)
            
        end %(while testAmp)
        
        %% if there is more than one maximum, test for sigificance of distances between maxima
        
        if numMaxima > 1
            
            %get p-values of distances between maxima
            pValue = mmfDistPV(maximaPos,varCovMat,numMaxima,numDegFree);
            
            %find largest p-value and decide whether to remove one kernel,
            %repeat the fit and test again for distances
            pValueMax = max(pValue(:));
            testDist = pValueMax > alphaD;
            
            %if any of the distances is not significant
            while testDist
                
                %find pair with maximum p-value
                [indx1,indx2] = find(pValue==pValueMax);
                indx1 = indx1(1);
                indx2 = indx2(1);
                
                %out of this pair, identify maximum with smaller amplitude
                ampBoth = maximaAmp([indx1 indx2],1);
                if ampBoth(1) < ampBoth(2)
                    indxBad = indx1;
                else
                    indxBad = indx2;
                end
                
                %determine which kernels to keep
                indx = setdiff(1:numMaxima,indxBad);
                
                %remove the information of the kernel to be discarded
                maximaPos = maximaPos(indx,:);
                maximaAmp = maximaAmp(indx,:);
                numMaxima = numMaxima - 1;
                
                %collect initial guesses and lower and upper bounds for fit
                [x0,lb,ub] = mmfInitGuessLowerUpperBounds(maximaPos(:,1:2),...
                    maximaAmp(:,1),bgAmp(1),psfSigma,clusterPixels,1);
                
                %calculate number of degrees of freedom in system
                numDegFree = size(clusterPixels,1) - 3*numMaxima - 1;
                
                %determine feature positions and amplitudes and estimate background
                %intensity, using nonlinear least squares data fitting
                %also get residuals and Jacobian matrix to calculate variance-covariance
                %of estimated parameters
                [solution,dummy,residuals,dummy,dummy,dummy,jacMat] = ...
                    lsqnonlin(@fitNGaussians2D,x0,lb,ub,options,imageC,...
                    clusterPixels,psfSigma);
                jacMat = full(jacMat);
                residuals = -residuals; %minus sign so that residuals = real image - model image
                residVar = sum(residuals.^2)/numDegFree;

                %calculate the parameters' variance-covariance matrix and get their
                %uncertainties
                varCovMat = residVar * inv(jacMat'*jacMat);
                standDevVec = sqrt(diag(varCovMat));
                
                %if nothing weird happened in the fit...
                if all(isreal(standDevVec))
                    
                    %extract estimate and std of background intensity and
                    %remove from vectors
                    bgAmp = [solution(end) standDevVec(end)];
                    solution = solution(1:end-1);
                    standDevVec = standDevVec(1:end-1);
                    
                    %reshape 3nx1 vectors "solution" and "standDevVec" into nx3 matrices
                    solution = reshape(solution,3,numMaxima);
                    solution = solution';
                    standDevVec = reshape(standDevVec,3,numMaxima);
                    standDevVec = standDevVec';
                    
                    %extract feature positions and amplitudes and their uncertainties
                    maximaPos = [solution(:,1:2) standDevVec(:,1:2)];
                    maximaAmp = [solution(:,3) standDevVec(:,3)];
                    
                    %if there is more than one maximum, test for
                    %sigificance of distances between maxima
                    if numMaxima > 1
                        
                        %get p-values of distances between maxima
                        pValue = mmfDistPV(maximaPos,varCovMat,numMaxima,numDegFree);
                        
                        %find largest p-value and decide whether to remove one kernel,
                        %repeat the fit and test again for distances
                        pValueMax = max(pValue);
                        testDist = pValueMax > alphaD;
                        
                    else
                        
                        testDist = 0;
                        
                    end %(if numMaxima > 1)
                    
                else %(if all(isreal(standDevVec)))
                    
                    testDist = 0; %stop fitting this cluster
                    keepCluster(i) = 0; %mark it to be discarded
                    
                end %(if all(isreal(standDevVec)))
                
            end %(while testDist)
            
        end %(if numMaxima > 1)
        
    else %(if all(isreal(standDevVec)))
        
        keepCluster(i) = 0; %mark this cluster to be discarded
        
    end %(if all(isreal(standDevVec)))
    
    %% repeat the fit with the final number of kernels to get a final
    %% estimate of the parameters and the residuals from the fit
    if keepCluster(i)
        
        %collect initial guesses and lower and upper bounds for fit
        [x0,lb,ub] = mmfInitGuessLowerUpperBounds(maximaPos(:,1:2),...
            maximaAmp(:,1),bgAmp(1),psfSigma,clusterPixels,1);
        
        %calculate number of degrees of freedom in system
        numDegFree = size(clusterPixels,1)-3*numMaxima-1;
        
        %determine feature positions and amplitudes and estimate background
        %intensity, using nonlinear least squares data fitting
        %also get residuals and Jacobian matrix to calculate variance-covariance
        %of estimated parameters
        [solution,dummy,residuals,dummy,dummy,dummy,jacMat] = ...
            lsqnonlin(@fitNGaussians2D,x0,lb,ub,options,imageC,...
            clusterPixels,psfSigma);
        jacMat = full(jacMat);
        residuals = -residuals; %minus sign so that residuals = real image - model image
        residVar = sum(residuals.^2)/numDegFree;
        
        %calculate the parameters' variance-covariance matrix and get their
        %uncertainties
        varCovMat = residVar * inv(jacMat'*jacMat);
        standDevVec = sqrt(diag(varCovMat));
        
        %extract estimate and std of background intensity and
        %remove from vectors
        bgAmp = [solution(end) standDevVec(end)];
        solution = solution(1:end-1);
        standDevVec = standDevVec(1:end-1);
        
        %reshape 3nx1 vectors "solution" and "standDevVec" into nx3 matrices
        solution = reshape(solution,3,numMaxima);
        solution = solution';
        standDevVec = reshape(standDevVec,3,numMaxima);
        standDevVec = standDevVec';
        
        %extract feature positions and amplitudes and their uncertainties
        maximaPos = [solution(:,1:2) standDevVec(:,1:2)];
        maximaAmp = [solution(:,3) standDevVec(:,3)];
        
        %store solution in clustersMMF
        clustersMMF(i).position = maximaPos;
        clustersMMF(i).amplitude = maximaAmp;
        clustersMMF(i).bgAmp = bgAmp;
        clustersMMF(i).numDegFree = numDegFree;
        clustersMMF(i).residuals = residuals;
        
    end %(if keepCluster(i))
    
end %(for i=length(clusters):-1:1)

%find the clusters to retain
indx = find(keepCluster);

%% final test on fit

%if there are clusters with significant signal ...
if ~isempty(indx)
    
    %retain only those clusters which have significant signals in them
    clustersMMF = clustersMMF(indx);
    
    %re-initialize keepClusters to use in next test
    keepCluster = ones(length(indx),1);
    
    %estimate background noise variance if not supplied
    if estimateBgNoise
    else
        bgNoiseVar = bgNoiseSigma^2;
    end
    
    %go over all clusters and check that their residuals are
    %comparable to the background noise
    for iCluster = 1 : length(clustersMMF)
        
        %get test statistic, which is F-distributed
        testStat = (sum(clustersMMF(iCluster).residuals.^2)/numDegFree)/...
            bgNoiseVar;
        
        %get p-value of test statistic
        %THIS NEEDS CHECKING IF I EVER WANT TO USE IT AGAIN
        %I THINK THE DEGREES OF FREEDOM ARE SWITCHED AROUND
        pValue = 1 - fcdf(testStat,numDegFree,numDegFree+3*size(clustersMMF(iCluster).bgAmp,1));
        
        %compare p-value to alpha
        %1-sided F-test: H0: F=1, H1: F>1
        if pValue < alphaF %if p-value is smaller than alpha, reject overall fit
            keepCluster(iCluster) = 0;
        end
        
    end %(for iCluster = 1 : length(clustersMMF))
    
    %determine clusters that pass the residuals test
    indx = find(keepCluster);
    
    if ~isempty(indx)
        
        %keep only clusters that pass the residuals test
        clustersMMF = clustersMMF(indx);
        
        %store information in structure "detectedFeatures"
        tmp = vertcat(clustersMMF.position);
        detectedFeatures.xCoord = [tmp(:,1) tmp(:,3)];
        detectedFeatures.yCoord = [tmp(:,2) tmp(:,4)];
        detectedFeatures.amp = vertcat(clustersMMF.amplitude);
        
    end
    
end %(if ~isempty(indx))

if isempty(indx)
    
    %store empty structures
    clustersMMF = struct('position',zeros(0,4),'amplitude',zeros(0,2),'bgAmp',zeros(0,2));
    detectedFeatures = struct('xCoord',zeros(0,2),'yCoord',zeros(0,2),'amp',zeros(0,2));
    
end

%save output if requested
if saveResults
    save('detectedFeatures','detectedFeatures','clustersMMF');
end


%% Visualization

if visual
    
    %make 3 layers out of original image (normalized)
    %     imageNorm = image/max(image(:));
    imageNorm = image/prctile(image(:),99.9);
    imageNorm(imageNorm>1) = 1;
    imageN3 = repmat(imageNorm,[1 1 3]);
    
    %place zeros in pixels of maxima from cands
    for i=1:length(clusters)
        for j=1:3
            imageN3(clusters(i).maximaPos(:,3)+(j-1)*numPixelsX*numPixelsY)=0;
        end
    end
    
    %place zeros in pixels of maxima from mixture-model fitting
    for i=1:length(clustersMMF)
        pos = (round(clustersMMF(i).position(:,1)-1))*numPixelsY ...
            + round(clustersMMF(i).position(:,2));
        for j=1:3
            imageN3(pos+(j-1)*numPixelsX*numPixelsY)=0;
        end
    end
    
    %label maxima from cands in blue
    for i=1:length(clusters)
        imageN3(clusters(i).maximaPos(:,3)+2*numPixelsX*numPixelsY)=1;
    end
    
    %label maxima from mixture-model fitting in red
    %a maximum from mixture-model fitting that falls in the same pixel
    %as that from cands will appear in magenta
    for i=1:length(clustersMMF)
        pos = (round(clustersMMF(i).position(:,1)-1))*numPixelsY ...
            + round(clustersMMF(i).position(:,2));
        imageN3(pos)=1;
    end
    
    %plot image
    imtool(imageN3);
    
end

%%%%% ~~ the end ~~ %%%%%

