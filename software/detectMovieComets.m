function detectMovieComets(movieData,varargin)
% detectMovieComets detect comets in a movie using successive thresholds
%
% SYNOPSIS detectMovieComets(movieData,paramsIn)
%
% INPUT   
%   movieData - A MovieData object describing the movie to be processed
%
%   paramsIn - Structure with inputs for optional parameters. The
%   parameters should be stored as fields in the structure, with the field
%   names and possible values as described below:
%
%   Possible Parameter Structure Field Names:
%       ('FieldName' -> possible values)
%
%       ('OutputDirectory' -> character string)
%       Optional. A character string specifying the directory to save the
%       detected comets to.
%
%       ('ChannelIndex' -> Positive integer scalar or vector)
%       Optional. The integer index of the plusTip channel(s). If not
%       input, all channels will be segmented. If multiple
%       channels are selected, comets are detected for each independently.
%
%       ('sigma1' -> Positive scalar) Optional.
%       Smallest standard deviation to be used when filtering the image 
%       with the difference of Gaussian. Default is 1.
% 
%       ('sigma2' -> Positive scalar) Optional.
%       Lagest standard deviation to be used when filtering the image 
%       with the difference of Gaussian. Default is 5.
% 
%       ('multFactThresh' -> Positive scalar) 
%       Multiplicative factor to apply to the standard deviation of the
%       filtered images to determine the threshold for detecting comets.
%       Optional. Default is 3
%
%       ('multFactorStepSize' -> Positive scalar) 
%       Multiplicative factor to apply to the standard deviation of the
%       filtered images to determine the difference of intensity between 
%       successive detection thresholds.
%       Optional. Default is 1
%
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


% Sebastien Besson, Sep 2011 (last modified Sep 2011)

%% ----------- Input ----------- %%

%Check input
ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('movieData', @(x) isa(x,'MovieData'));
ip.addOptional('paramsIn',[], @isstruct);
ip.parse(movieData,varargin{:});
paramsIn=ip.Results.paramsIn;

%Get the indices of any previous speckle detection processes                                                                     
iProc = movieData.getProcessIndex('CometDetectionProcess',1,0);

%If the process doesn't exist, create it
if isempty(iProc)
    iProc = numel(movieData.processes_)+1;
    movieData.addProcess(CometDetectionProcess(movieData,...
        movieData.outputDirectory_));                                                                                                 
end
cometDetProc = movieData.processes_{iProc};
%Parse input, store in parameter structure
p = parseProcessParams(cometDetProc,paramsIn);

%% --------------- Initialization ---------------%%
if feature('ShowFigureWindows')
    wtBar = waitbar(0,'Initializing...','Name',cometDetProc.getName());
else 
    wtBar = -1;
end

% Reading various constants
imDirs  = movieData.getChannelPaths;
bitDepth = movieData.camBitdepth_;
nFrames = movieData.nFrames_;
maxIntensity =(2^bitDepth-1);
nChan =  numel(movieData.channels_); 

% Set up the input directories
inFilePaths = cell(1,nChan);
for j = p.ChannelIndex,  inFilePaths{1,j} = imDirs{j}; end
cometDetProc.setInFilePaths(inFilePaths);
    
% Set up the output directories
outFilePaths = cell(1,nChan);

for i = p.ChannelIndex;    
    %Create string for current directory
    outFilePaths{1,i} = [p.OutputDirectory filesep 'channel_' num2str(i) '.mat'];
end
filteredImagesDirectory = [p.OutputDirectory filesep 'filtered_images'];
if isdir(filteredImagesDirectory), rmdir(filteredImagesDirectory, 's'); end
mkClrDir(p.OutputDirectory);
cometDetProc.setOutFilePaths(outFilePaths);

%% --------------- Comet detection ---------------%%% 
disp('Starting detecting comets...')

% Get region of interest
roiMask = movieData.getROIMask;
allMovieInfo(nFrames,nChan) = struct('xCoord', [], 'yCoord', [], 'amp', [],...
        'int',[],'ecc',[]);
    
for iChan= p.ChannelIndex

    % Create temporary directory for storing filtered images
    if ~isdir(filteredImagesDirectory), mkdir(filteredImagesDirectory); end
    if ishandle(wtBar), waitbar(0,wtBar,'Loading image stack'); end    
    
    % get difference of Gaussians image for each frame and standard deviation
    % of the cell background, stored in stdList
    stdList=NaN(nFrames,1);
    
    logMsg=['Filtering images for comet detection for channel ' num2str(iChan)];
    progressText(0,logMsg);
    if ishandle(wtBar), waitbar(0,wtBar,logMsg); end  
    for i = p.firstFrame:p.lastFrame
        % Loading image
        im = double(movieData.channels_(iChan).loadImage(i))/maxIntensity;
        
        % Combine region of interest with preselected segmentation output
        if ~isempty(p.MaskProcessIndex)
            mask= roiMask(:,:,i) & movieData.processes_{p.MaskProcessIndex}.loadChannelOutput(iChan,i);
        else
            mask=roiMask(:,:,i);
        end
        
        % SB: copied and pasted from plusTipCometDetector
%       
%         % create kernels for gauss filtering
%         blurKernelLow  = fspecial('gaussian', 21, p.sigma1);
%         blurKernelHigh = fspecial('gaussian', 21, p.sigma2);
%         
%         Wlow = imfilter(double(mask), blurKernelLow);
%         Whigh = imfilter(double(mask), blurKernelHigh);
%         
%         % use subfunction that calls imfilter to take care of edge effects
%         % for now don't apply roiMask
%         im(~mask)=0;
%         lowPass = imfilter(im, blurKernelLow)./Wlow;
%         lowPass(~mask)=NaN;
%         highPass = imfilter(im, blurKernelHigh)./Whigh;
%         highPass(~mask)=NaN;
%         filterDiff=lowPass-highPass;

        % Difference of Gaussians
        filterDiff=filterGauss2D(im,p.sigma1)-filterGauss2D(im,p.sigma2);
        filterDiff(~mask)=NaN;
        
        % Save filtered images on disk (avoid memory errors)
        save(fullfile(filteredImagesDirectory,['filterDiff_' num2str(i) '.mat']),'filterDiff');

        % Perform maximum filter and mask out significant pixels
        thFilterDiff = ordfilt2(filterDiff,9,ones(3,3));
        threshold = thresholdOtsu(thFilterDiff)/3 + ...
            thresholdRosin(thFilterDiff)*2/3;
        stdList(i)=nanstd(filterDiff(thFilterDiff<threshold));
%         stdList(i)=nanstd(filterDiff(:));

        % Update progress status
        frac = (i-p.firstFrame+1)/(p.lastFrame-p.firstFrame+1);
        progressText(frac,'Filtering images for comet detection');
        if ishandle(wtBar) && mod(i,5)==0, waitbar(frac,wtBar,logMsg);  end  
    end
    
    meanStd = arrayfun(@(x) nanmean(stdList(max(1,x-1):min(nFrames,x+1))),1:nFrames);
    % meanStd = smooth(stdList,3); % Discrepancy for endpoint
    
    % loop thru frames and detect
    logMsg='Detecting comets';
    progressText(0,logMsg);
    if ishandle(wtBar), waitbar(0,wtBar,logMsg); end  
    for i = p.firstFrame:p.lastFrame
        % Reload band-pass filtered images
        s=load(fullfile(filteredImagesDirectory, ['filterDiff_' num2str(i) '.mat']));
        filterDiff=s.filterDiff;
        stepSize=p.multFactorStepSize*meanStd(i);
        thresh= p.multFactorThresh*meanStd(i);
        
        % Detect comets using watershed detection
        allMovieInfo(i,iChan) = detectComets(filterDiff,stepSize,thresh);
             
        % Update progress status
        frac = (i-p.firstFrame+1)/(p.lastFrame-p.firstFrame+1);
        progressText(frac,'Detecting comets');
        if ishandle(wtBar) && mod(i,5)==0, waitbar(frac,wtBar,logMsg); end
    end
    
    rmdir(filteredImagesDirectory,'s')
    movieInfo=allMovieInfo(:,iChan);
    save(outFilePaths{1,iChan} ,'movieInfo','stdList');
    
end

if ishandle(wtBar), close(wtBar); end
disp('Finished detecting comets');