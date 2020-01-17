function detectMoviePointSources3D(movieDataOrProcess,varargin)
% detectMoviePointSource detect diffraction-limited objects in a movie
%
% detectMoviePointSources 
%
% SYNOPSIS detectMoviePointSources(movieData,paramsIn)
%
% INPUT   
%   movieData - A MovieData object describing the movie to be processed
%
%   paramsIn - Structure with inputs for optional parameters. The
%   parameters should be stored as fields in the structure, with the field
%   names and possible values as described below
%
% OUTPUT   
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

% Joy Xu / Sebastien Besson, July 2014
% Andrew R. Jamieson, March 2017

%% ----------- Input ----------- %%

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('movieDataOrProcess', @isProcessOrMovieData);
ip.addOptional('paramsIn',[], @isstruct);
ip.addParameter('UseIntersection',false, @islogical);
ip.addParamValue('ROI',[], @isnumeric);
ip.parse(movieDataOrProcess,varargin{:});
ip.KeepUnmatched = true;
paramsIn = ip.Results.paramsIn;
ROI = ip.Results.ROI;

% Get MovieData object and Process
[movieData, pointSourceDetProc3D] = getOwnerAndProcess(movieDataOrProcess,'PointSourceDetectionProcess3D',true);
assert(movieData.is3D, 'detectMoviePointSources3D is specifically designed for 3D images. Please use detectMoviePointSources for 2D images.');

%Parse input, store in parameter structure
p = parseProcessParams(pointSourceDetProc3D, paramsIn);

%% --------------- Initialization ---------------%%
if feature('ShowFigureWindows')
    wtBar = waitbar(0,'Initializing...','Name',pointSourceDetProc3D.getName());
else
    wtBar=-1;
end

% Reading various constants
nFrames = movieData.nFrames_;
nChanDet = numel(p.ChannelIndex);
imSize = movieData.imSize_;

%Find the  the segmentation process.
if isempty(p.MaskProcessIndex) && ~isempty(p.MaskChannelIndex);
    p.MaskProcessIndex =movieData.getProcessIndex('MaskProcess',1,1);
end    

if ~isempty(p.MaskProcessIndex) && ~isempty(p.MaskChannelIndex)
    maskProc = movieData.processes_{p.MaskProcessIndex};
    if ~all(maskProc.checkChannelOutput(p.MaskChannelIndex))
        error('All channels must be segmented!')
    end
    
    %Create mask directory if several masks need to be merged
    if length(p.MaskChannelIndex) > 1 %&& p.UseIntersection
        close(wtBar); error('MaskIntersectionProcess not yet supported in 3D. Please sepecify a single MaskChannelIndex.'); 
%         %Get the indices of any previous mask intersection process
%         iMaskIntProc = movieData.getProcessIndex('MaskIntersectionProcess',1,0);
%         
%         %If the process doesn't exist, create it
%         if isempty(iMaskIntProc)
%             iMaskIntProc = numel(movieData.processes_)+1;
%             movieData.addProcess(MaskIntersectionProcess(movieData,p.OutputDirectory));
%         end
%         maskIntProc = movieData.processes_{iMaskIntProc};
%         
%         %Set up the parameters for mask intersection
%         maskIntParams.ChannelIndex = p.MaskChannelIndex;
%         maskIntParams.SegProcessIndex = p.MaskProcessIndex;
%         
%         parseProcessParams(maskIntProc,maskIntParams);
%         maskIntProc.run;
%         
%         %Then use this mask process
%         maskProc = maskIntProc;    
    end 
    % Get mask directory and names
    maskDir = maskProc.outFilePaths_(p.MaskChannelIndex);
else
    maskProc = [];
end

%Check the input processes if any and get loader handles for each channel
imDirs = cell(1, nChanDet);
imLoader = cell(1, nChanDet);
for j = 1:nChanDet
    if p.InputImageProcessIndex(j) > 0
        %Check the specified input process
        assert(movieData.processes_{p.InputImageProcessIndex(j)}.checkChannelOutput(p.ChannelIndex(j)),['No valid output for input process specified for channel ' num2str(p.ChannelIndex(j))]);
        imDirs{p.ChannelIndex(j)} = movieData.processes_{p.InputImageProcessIndex(j)}.outFilePaths_{p.ChannelIndex(j)};
        imLoader{p.ChannelIndex(j)} = @(f)(movieData.processes_{p.InputImageProcessIndex(j)}.loadChannelOutput(p.ChannelIndex(j),f));                         
    else
        %If raw data specified
        imDirs{p.ChannelIndex(j)} = movieData.channels_(p.ChannelIndex(j)).channelPath_;
        imLoader{p.ChannelIndex(j)} = @(f)(movieData.channels_(p.ChannelIndex(j)).loadStack(f));
        
    end    
end

% Set up the input directories
inFilePaths = cell(1,numel(movieData.channels_));
for j = 1:numel(p.ChannelIndex)
    inFilePaths{1,p.ChannelIndex(j)} = imDirs{p.ChannelIndex(j)};
    if ~isempty(p.MaskProcessIndex) && ~isempty(p.MaskChannelIndex)
        inFilePaths{2,p.ChannelIndex(j)} = maskDir;
    end
end
pointSourceDetProc3D.setInFilePaths(inFilePaths);
    
% Set up the output directories
outFilePaths = cell(1, numel(movieData.channels_));
for i = p.ChannelIndex;    
    %Create string for current directory
    outFilePaths{1,i} = [p.OutputDirectory filesep 'channel_' num2str(i) '.mat'];
    %outFilePaths{2,i} = [p.OutputDirectory filesep 'channel_detectionLabRef' num2str(i) '.mat'];
end
mkClrDir(p.OutputDirectory)
pointSourceDetProc3D.setOutFilePaths(outFilePaths);

%Get ROI mask if any.
%roiMask = movieData.getROIMask;



%% --------------- Add optional auio-estimation of PSF sigma ---------------%%% 


%% --------------- Point source detection ---------------%%% 

disp('Starting detecting diffraction-limited objects');
logMsg = @(chan) ['Please wait, detecting diffraction-limited objects for channel ' num2str(chan)];
timeMsg = @(t) ['\nEstimated time remaining: ' num2str(round(t)) 's'];
tic;
nChan = length(p.ChannelIndex);
nTot = nChan*nFrames;

ZXRatio=movieData.pixelSizeZ_/movieData.pixelSize_;

for i = 1:numel(p.ChannelIndex)

    iChan = p.ChannelIndex(i);
    % Set up for parfor
    detP = splitPerChannelParams(p, iChan);

    processFrames = detP.frameRange(1):detP.frameRange(2);
    labels = cell(1,numel(processFrames));
    movieInfo(numel(processFrames), 1) = struct('xCoord', [], 'yCoord',[],'zCoord', [], 'amp', [], 'int',[]);
                     
    % Log display
    disp(logMsg(iChan))
    disp(imDirs{1,iChan});
    if ~isempty(p.MaskProcessIndex) && ~isempty(p.MaskChannelIndex)
        disp(sprintf('Using mask from: %s', maskDir{1}));
    end
    disp('Results will be saved under:')
    disp(outFilePaths{1,iChan});
    
    
    %Set up parameter structure for detection on this channel
    if (strcmp(detP.algorithmType,'pointSourceAutoSigma') ...
            ||strcmp(detP.algorithmType,'pointSourceAutoSigmaFit') ...
            ||strcmp(detP.algorithmType,'pSAutoSigmaMarkedWatershed') ...
            ||strcmp(detP.algorithmType,'pointSourceAutoSigmaMixture') ...
            ||strcmp(detP.algorithmType,'pointSourceAutoSigmaLM') ...
            ||strcmp(detP.algorithmType,'pointSourceAutoSigmaFitSig') ...
            ||strcmp(detP.algorithmType,'pSAutoSigmaWatershed')) ...
            %&&(isempty(detP.filterSigma))
        volList=[];
        for fIdx=1:length(processFrames)
            volList=[volList {double(movieData.getChannel(iChan).loadStack(fIdx))}];
        end
        sigmasPSF=getGaussianPSFsigmaFrom3DData(volList,'Display',logical(p.showAll));
        arrayfun(@(x) disp(['Estimed scales: ' num2str(x)]), sigmasPSF);
    else
        sigmasPSF = detP.filterSigma;
    end

   
    parfor frameIdx = 1:length(processFrames)    
                
        timePoint = processFrames(frameIdx);
        disp(['Processing time point ' num2str(frameIdx,'%04.f')])
        
        % loading the entire stack
        % Hunter's approach
        detP_pf = detP;
        vol = double(imLoader{iChan}(timePoint));
        
        % Philippe's approach
        % vol = double(movieData.getChannel(iChan).loadStack(timePoint)); #
        volSize = size(vol);
        lab = [];
        
        % find mask offset (WARNING works only for cubic mask)
        if (~isempty(ROI))
            [maskMinX,maskMinY,maskMinZ]=ind2sub(size(ROI), find(ROI,1));
            [maskMaxX,maskMaxY,maskMaxZ]=ind2sub(size(ROI), find(ROI,1,'last'));
        end
        
        if(~isempty(ROI))
            tmp = nan(1+[maskMaxX,maskMaxY,maskMaxZ]-[maskMinX,maskMinY,maskMinZ]);
            tmp(:) = vol(ROI>0);
            vol = tmp;
        end
        
%         %!!!!! fix the mask for 3d data!!!!!%
        if ~isempty(maskProc)
            currMask = maskProc.loadChannelOutput(p.MaskChannelIndex,timePoint); %& roiMask(:,:,frameIdx); %istack is the stack index!!!! test for validity!!!!
            detP_pf.Mask =  currMask;
        else
            detP_pf.Mask = [];
        end

       
        switch detP.algorithmType
            case {'pointSourceLM',...
                  'pointSource',...
                  'pointSourceAutoSigma',...
                  'pointSourceAutoSigmaFit',...
                  'pointSourceFit',...
                  'pSAutoSigmaMarkedWatershed',... 
                  'pointSourceAutoSigmaLM',...     
                  'pSAutoSigmaWatershed'}
                [pstruct, mask, imgLM, imgLoG] = pointSourceDetection3D(vol, sigmasPSF, detP_pf);
                
                switch detP.algorithmType
                      case {'pointSource','pointSourceAutoSigma'}
                       lab = double(mask); 
                        movieInfo(frameIdx) = labelToMovieInfo(double(mask),vol);
                      case {'pointSourceLM','pointSourceAutoSigmaLM'}
                        lab = double(mask); 
                        movieInfo(frameIdx) = pointCloudToMovieInfo(imgLM,vol);  
                      case 'pSAutoSigmaMarkedWatershed'
                        wat = markedWatershed(vol,sigmasPSF,0);
                        wat(mask==0) = 0;       
                        lab = double(wat); 
                        movieInfo(frameIdx) = labelToMovieInfo(double(wat),vol);
                      case 'pSAutoSigmaWatershed'
                        wat = watershed(-vol.*mask);
                        wat(mask==0) = 0;
                        lab = double(wat); 
                        movieInfo(frameIdx) = labelToMovieInfo(double(wat),vol);
                      case {'pointSourceFit','pointSourceAutoSigmaFit'}
                        movieInfo(frameIdx) = pstructToMovieInfo(pstruct);
                        lab = double(mask); %.*imgLoG;
                    otherwise
                end
            
            case {'pointSourceAutoSigmaMixture'}
                detPt = detP_pf;
                detPt.FitMixtures = true;
                [pstruct, mask, imgLM, imgLoG] = pointSourceDetection3D(vol,sigmasPSF,detPt);
                movieInfo(frameIdx) = pstructToMovieInfo(pstruct);
                lab = double(mask); % adjust label

            case {'pointSourceAutoSigmaFitSig'}
                detPt = detP_pf;
                detPt.Mode = 'xyzAcsr';
                [pstruct,mask,imgLM,imgLoG] = pointSourceDetection3D(vol,sigmasPSF,detPt);
                movieInfo(frameIdx) = pstructToMovieInfo(pstruct);
                lab = double(mask); % adjust label

            case {'pointSourceFitSig'}
                detPt = detP_pf;
                detPt.Mode = 'xyzAcsr';
                [pstruct,mask,imgLM,imgLoG] = pointSourceDetection3D(vol,sigmasPSF,detPt);
                movieInfo(frameIdx) = pstructToMovieInfo(pstruct);
                lab = double(mask); % adjust label
                
            case {'pointSourceAutoSigmaFitSig'}
                    detPt = detP_pf;
                    detPt.Mode = 'xyzAcsr';
                    [pstruct,mask,imgLM,imgLoG] = pointSourceDetection3D(vol,sigmasPSF,detPt);
                    movieInfo(frameIdx) = pstructToMovieInfo(pstruct);
                        lab = double(mask); % adjust label                

            % ----------------------------------------------------------------------------------------
            case 'watershedApplegate'
                filterVol = filterGauss3D(vol,detP.highFreq)-filterGauss3D(vol,detP.lowFreq);
                [movieInfo(frameIdx),lab] = detectComets3D(filterVol,detP.waterStep,detP.waterThresh,[1 1 1]);

            case 'watershedApplegateAuto'
                filterVol = filterGauss3D(vol,detP.highFreq)-filterGauss3D(vol,detP.lowFreq);
                thresh = QDApplegateThesh(filterVol,detP.showAll);
                [movieInfo(frameIdx),lab] = detectComets3D(filterVol,detP.waterStep,thresh,[1 1 1]);

            % ----------------------------------------------------------------------------------------
            case 'bandPassWatershed'
                filterVol = filterGauss3D(vol,detP.highFreq) - filterGauss3D(vol,detP.lowFreq);
                label = watershed(-filterVol); 
                label(filterVol < detP.waterThresh) = 0;
                lab = label;
                movieInfo(frameIdx)=labelToMovieInfo(label,filterVol);

            % ----------------------------------------------------------------------------------------
            case 'watershedMatlab'
                label = watershed(-vol); 
                label(vol<p.waterThresh) = 0;
                [dummy,nFeats] = bwlabeln(label);
                lab = label;
                movieInfo(frameIdx) = labelToMovieInfo(label,vol);
                
            case 'markedWatershed'
                label = markedWatershed(vol,sigmasPSF,p.waterThresh);
                movieInfo(frameIdx) = labelToMovieInfo(label, vol);
              % ----------------------------------------------------------------------------------------                
            
            otherwise 
                error('Supported detection algorithm method:');
        end
        
        if isfield(p, 'isoCoord') && p.isoCoord
            movieInfo(frameIdx).zCoord(:,1)=movieInfo(frameIdx).zCoord(:,1)*ZXRatio;
        end

        if(~isempty(ROI))
            tmplab=zeros(volSize);
            tmplab(ROI>0)=lab;
            lab=tmplab; 
            movieInfo(frameIdx).xCoord(:,1)=movieInfo(frameIdx).xCoord(:,1)+maskMinY-1;
            movieInfo(frameIdx).yCoord(:,1)=movieInfo(frameIdx).yCoord(:,1)+maskMinX-1;
            movieInfo(frameIdx).zCoord(:,1)=movieInfo(frameIdx).zCoord(:,1)+maskMinZ-1;
        end    
        labels{frameIdx}=lab;
    end %%%% end parfor (frame loop)
    
    if ~exist('movieInfo','var')
        %in the case that no channels/frames had detected points
        movieInfo = [];
    end
    %     for fIdx=1:length(detectionLabRef)
%         detectionLabRef(fIdx).zCoord(:,1)=detectionLabRef(fIdx).zCoord(:,1)*movieData.pixelSizeZ_/movieData.pixelSize_;
%     end

    save(outFilePaths{1,iChan}, 'movieInfo', 'labels');
    %save(outFilePaths{2,iChan}, 'detectionLabRef');
    
%     clear movieInfo detectionLabRef;
    clear movieInfo;

end %%% channel loop

% Close waitbar
if ishandle(wtBar), close(wtBar); end
disp('Finished detecting diffraction-limited objects!')
movieData.save;


function movieInfo= labelToMovieInfo(label,vol)
    [feats,nFeats] = bwlabeln(label);
    featsProp = regionprops(feats,vol,'Area','WeightedCentroid','MeanIntensity','MaxIntensity','PixelValues');

    movieInfo=struct('xCoord',[],'yCoord',[],'zCoord',[],'amp',[],'int',[]);
    
    % centroid coordinates with 0.5 uncertainties
    tmp = vertcat(featsProp.WeightedCentroid)-1;
    if ~isempty(featsProp)
        xCoord = [tmp(:,1) 0.5*ones(nFeats,1)]; yCoord = [tmp(:,2) 0.5*ones(nFeats,1)]; zCoord = [tmp(:,3) 0.5*ones(nFeats,1)];
        amp=[vertcat(featsProp.MaxIntensity) 0.5*ones(nFeats,1)];
    
        movieInfo.xCoord= xCoord;movieInfo.yCoord=yCoord;movieInfo.zCoord=zCoord;
        movieInfo.amp=amp;
        movieInfo.int=amp;
    end

    % u-track formating


function movieInfo= pointCloudToMovieInfo(imgLM,vol)
    lmIdx = find(imgLM~=0);
    [lmy,lmx,lmz] = ind2sub(size(vol), lmIdx);
    N=length(lmy);
    % centroid coordinates with 0.5 uncertainties
    xCoord = [lmx 0.5*ones(N,1)]; yCoord = [lmy 0.5*ones(N,1)]; zCoord = [lmz 0.5*ones(N,1)];
    amp=[vol(lmIdx) 0.5*ones(N,1)];

    % u-track formating
    movieInfo=struct('xCoord',[],'yCoord',[],'zCoord',[],'amp',[],'int',[]);
    movieInfo.xCoord= xCoord;movieInfo.yCoord=yCoord;movieInfo.zCoord=zCoord;
    movieInfo.amp=amp;
    movieInfo.int=amp;

function mkdir(path)
system(['mkdir -p ' path]);


function movieInfo= pstructToMovieInfo(pstruct)
    if(~isempty(pstruct))
        movieInfo.xCoord = [pstruct.x' pstruct.x_pstd'];
        movieInfo.yCoord = [pstruct.y' pstruct.y_pstd'];
        movieInfo.zCoord = [pstruct.z' pstruct.z_pstd'];
        movieInfo.amp = [pstruct.A' pstruct.A_pstd'];
        movieInfo.int= [pstruct.A' pstruct.A_pstd'];
    else
        movieInfo.xCoord=[];
        movieInfo.yCoord=[];
        movieInfo.zCoord=[];
        movieInfo.amp=[];
        movieInfo.int=[];
    end
%     movieInfo.sigmaX = [pstruct.s' pstruct.s_pstd'];
%     movieInfo.sigmaY = [pstruct.s' pstruct.s_pstd'];
%     movieInfo.sigmaZ = [pstruct.s' pstruct.s_pstd'];
%     movieInfo.bkg = [pstruct.c' pstruct.c_pstd'];

function threshNoise= QDApplegateThesh(filterDiff,show)
    % Perform maximum filter and mask out significant pixels
    thFilterDiff = locmax3d(filterDiff,1);
    threshold = thresholdOtsu(thFilterDiff)/3 + ...
        thresholdRosin(thFilterDiff)*2/3;
    std=nanstd(filterDiff(thFilterDiff<threshold));
    threshNoise= 3*std;

    if(show)
        figure();hist(thFilterDiff,100),vline([threshNoise, threshold],['-b','-r']);
    end 

