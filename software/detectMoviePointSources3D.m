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
% Copyright (C) 2019, Danuser Lab - UTSouthwestern 
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
% Philippe Roudot 2016-2018

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

% precondition / error checking
if isa(pointSourceDetProc3D, 'PointSourceDetectionProcess3DDynROI')
    buildDynROIProcId = movieData.getProcessIndex('BuildDynROIProcess'); % if numel(buildDynROIProcId) > 1, popup window will show and let user to choose which BuildDynROIProcess.
    if isempty(buildDynROIProcId)
        error("BuildDynROIProcess needs to be done before run PointSourceDetectionProcess3DDynROI.")
    elseif ~ismember(1, movieData.processes_{buildDynROIProcId}.funParams_.ChannelIndex)
        error("Channel 1 in BuildDynROIProcess needs to be analyzed before run PointSourceDetectionProcess3DDynROI.")
    end
end

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
        if(~isempty(p.processBuildDynROI)&&p.processBuildDynROI.isSwaped())
            movieData=p.processBuildDynROI.loadSwap(1);
            imDirs{p.ChannelIndex(j)} = movieData.channels_(p.ChannelIndex(j)).channelPath_;
            imLoader{p.ChannelIndex(j)} = @(f)(movieData.channels_(p.ChannelIndex(j)).loadStack(f));
        else
            imDirs{p.ChannelIndex(j)} = movieData.channels_(p.ChannelIndex(j)).channelPath_;
            imLoader{p.ChannelIndex(j)} = @(f)(movieData.channels_(p.ChannelIndex(j)).loadStack(f));
        end
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
    if ~isempty(p.processBuildDynROI) && isa(pointSourceDetProc3D, 'PointSourceDetectionProcess3DDynROI')
        outFilePaths{2,i} = [p.OutputDirectory filesep 'channel_DynROIRef' num2str(i) '.mat'];
    end
end
mkClrDir(p.OutputDirectory);
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
    if(isempty(detP.frameRange))
        detP.frameRange(2)=movieData.nFrames_;
        detP.frameRange(1)=1;
    end
        
    processFrames = detP.frameRange(1):detP.frameRange(2);
    labels = cell(1,numel(processFrames));

    %% saving maps and detection for debug purposes
    debugVols = cell(1,numel(processFrames));

    movieInfo(numel(processFrames), 1) = struct('xCoord', [], 'yCoord',[],'zCoord', [], 'amp', [], 'int',[]);
    debugPos(numel(processFrames))=Detections;
    labelSegPos(numel(processFrames))=Detections;

                     
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
        detP_pf = detP;
        vol = double(imLoader{iChan}(timePoint));
                lab = [];
        debugVol=[];
        energyMap=[];

        %% build cuboid mask from dynROI
        ROI=[];
        maskMinY=[];maskMinX=[];maskMinZ=[];
        maskMaxX=[];maskMaxY=[];maskMaxZ=[];
        if(~isempty(detP.processBuildDynROI))
            tmp=detP.processBuildDynROI.loadFileOrCache(); % try initDynROIs
            dynROICell=tmp{1}.dynROICell;
            dynROI=dynROICell{1};
            [BBmin,BBmax]=dynROI.getBoundingBox();
            maskMinX=BBmin(2); maskMinY=BBmin(1); maskMinZ=ceil(BBmin(3)/ZXRatio); 
            maskMaxX=BBmax(2); maskMaxY=BBmax(1); maskMaxZ=floor(BBmax(3)/ZXRatio); 
        end

        if(~isempty(detP.processBuildDynROI)&&~(detP.processBuildDynROI.isSwaped()))
            % ROI=false(size(vol));
            % ROI(max(maskMinX,1):min(maskMaxX,end),max(1,maskMinY):min(end,maskMaxY),max(maskMinZ,1):min(end,maskMaxZ))=true;
            % disp('getBoundingBox');
            % [maskMinX,maskMinY,maskMinZ]=ind2sub(size(ROI), find(ROI,1));
            % [maskMaxX,maskMaxY,maskMaxZ]=ind2sub(size(ROI), find(ROI,1,'last'));
            % tmp = nan(1+[maskMaxX,maskMaxY,maskMaxZ]-[maskMinX,maskMinY,maskMinZ]);
            % size(tmp)
            % tmp(:) = vol(ROI);
            % origVol=vol;
            % vol = tmp;

            tmp=detP.processBuildDynROI.loadFileOrCache(); % try initDynROIs
            dynROICell=tmp{1}.dynROICell;
            dynROI=dynROICell{1};
            [vol,ROI]=dynROI.cropBoundingVol(vol,ZXRatio);
        end


        samplePos=[];
        if(~isempty(detP.samplePos))  
            samplePos=detP.samplePos(timePoint);
            if (~isempty(ROI))
                samplePos.addOffset(-maskMinY+1,-maskMinX+1,-maskMinZ+1);
            end
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
                  'pSWatershed',...
                  'pSAutoSigmaWatershed'}

                [pstruct, mask,imgLM] = pointSourceDetection3D(vol, sigmasPSF, detP_pf);
                
                switch detP.algorithmType
                      case {'pointSource','pointSourceAutoSigma'}
                       lab = uint8(mask); 
                        movieInfo(frameIdx) = labelToMovieInfo(double(mask),vol);
                      case {'pointSourceLM','pointSourceAutoSigmaLM'}
                        lab = uint8(mask); 
                        movieInfo(frameIdx) = pointCloudToMovieInfo(imgLM,vol);  
                      case 'pSAutoSigmaMarkedWatershed'
                        wat = markedWatershed(vol,sigmasPSF,0);
                        wat(mask==0) = 0;       
                        lab = uint8(wat); 
                        movieInfo(frameIdx) = labelToMovieInfo(double(wat),vol);
                      case {'pSAutoSigmaWatershed','pSWatershed'}
                        wat = watershed(-vol.*mask);
                        wat(mask==0) = 0;
                        lab = uint8(wat); 
                        movieInfo(frameIdx) = labelToMovieInfo(double(wat),vol);
                      case {'pointSourceFit','pointSourceAutoSigmaFit'}
                        movieInfo(frameIdx) = pstructToMovieInfo(pstruct);
                        lab = uint8(mask); %.*imgLoG;
                    otherwise
                end
            
            case {'pointSourceLM','pointSourceAutoSigmaLM'}
                [pstruct, mask, imgLM] = pointSourceDetection3D(vol, sigmasPSF, detP_pf);
                lab = uint8(mask); 
                movieInfo(frameIdx) = pointCloudToMovieInfo(imgLM,vol); 
            case {'pointSourceAutoSigmaMixture'}
                detPt = detP_pf;
                detPt.FitMixtures = true;
                [pstruct, mask, imgLM] = pointSourceDetection3D(vol,sigmasPSF,detPt);
                movieInfo(frameIdx) = pstructToMovieInfo(pstruct);
                lab = uint8(mask); % adjust label

            case {'pointSourceAutoSigmaFitSig'}
                detPt = detP_pf;
                detPt.Mode = 'xyzAcsr';
                [pstruct,mask,imgLM] = pointSourceDetection3D(vol,sigmasPSF,detPt);
                movieInfo(frameIdx) = pstructToMovieInfo(pstruct);
                lab = uint8(mask); % adjust label

            case {'multiscaleDetection'}
                detPt = detP_pf;
                [imgLM,scaleVol]=multiscaleStochasticFiltering(vol,1./ZXRatio,detPt);
                movieInfo(frameIdx) = pointCloudToMovieInfo(imgLM,vol);  
                lab=imgLM;
            case {'multiscaleDetectionDebug'}
                detPt = detP_pf;
                [pos,labelSeg,energyMap]=multiscaleStochasticFiltering(vol,1./ZXRatio,detPt,'deepakImplementation',true,'samplePos',samplePos);
                movieInfo(frameIdx) = pos;  
                lab=labelSeg;
                debugVol={energyMap};

            case {'multiscaleDetectionDebugSeparable'}
                detPt = detP_pf;
                [pos,labelSeg,energyMap]=multiscaleStochasticFiltering(vol,1./ZXRatio,detPt,'deepakImplementation',false,'samplePos',samplePos);
                movieInfo(frameIdx) = pos;  
                lab=labelSeg;
                debugVol={energyMap};
            case {'labelToDetectionStruct'}
                movieInfo(frameIdx) = labelToMovieInfo(vol,vol);
                lab=vol;
                
            case {'fastBlobFinder'}
                detPt = detP_pf;
                [pos,imgLoG]=fastBlobFinder(vol,sigmasPSF, detPt);
                movieInfo(frameIdx) = Detections().initFromPosMatrices(pos,pos).getStruct();
                lab=[];

            case {'pointSourceFitSig'}
                detPt = detP_pf;
                detPt.Mode = 'xyzAcsr';
                [pstruct,mask,imgLM] = pointSourceDetection3D(vol,sigmasPSF,detPt);
                movieInfo(frameIdx) = pstructToMovieInfo(pstruct);
                lab = uint8(mask); % adjust label
                
            case {'pointSourceAutoSigmaFitSig'}
                    detPt = detP_pf;
                    detPt.Mode = 'xyzAcsr';
                    [pstruct,mask,imgLM] = pointSourceDetection3D(vol,sigmasPSF,detPt);
                    movieInfo(frameIdx) = pstructToMovieInfo(pstruct);
                        lab = uint8(mask); % adjust label                

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

        % rawMovieInfo=movieInfo;

        if(~isempty(detP.processBuildDynROI))
            if(~isempty(movieInfo(frameIdx).xCoord))
            movieInfo(frameIdx).xCoord(:,1)=movieInfo(frameIdx).xCoord(:,1)+maskMinY-1;
            movieInfo(frameIdx).yCoord(:,1)=movieInfo(frameIdx).yCoord(:,1)+maskMinX-1;
            movieInfo(frameIdx).zCoord(:,1)=movieInfo(frameIdx).zCoord(:,1)+maskMinZ-1;
            end
            if(~isempty(detP.samplePos))  
                samplePos.addOffset(maskMinY-1,maskMinX-1,maskMinZ-1);
            end
            
        end  


        if(~isempty(detP.saveMaskFilePattern))
            stackWrite(uint16(lab),sprintfPath(detP.saveMaskFilePattern,frameIdx));

        end
        labDet=Detections().initFromPointCloud(lab,lab,1);
        if(~isempty(detP.processBuildDynROI))
            labDet.addOffset(maskMinY-1,+maskMinX-1,maskMinZ-1);
        end

        debugDet=Detections();
        if(~isempty(energyMap))
            debugDet=Detections().initFromPointCloud(energyMap,energyMap,1);
            if(~isempty(detP.processBuildDynROI))
                debugDet.addOffset(maskMinY-1,maskMinX-1,maskMinZ-1);
            end
        end

   
        if isfield(p, 'isoCoord') && p.isoCoord && (~isempty(movieInfo(frameIdx).zCoord))
            movieInfo(frameIdx).zCoord(:,1)=movieInfo(frameIdx).zCoord(:,1)*ZXRatio;
            labDet.zCoord(:,1)=labDet.zCoord(:,1)*ZXRatio;
            if(~isempty(debugDet))
                if(~isempty(debugDet.zCoord))
                    debugDet.zCoord(:,1)=debugDet.zCoord(:,1)*ZXRatio;
                end
            end
        end
        % MDout=dynROI.swapRawBoundingBox(movieData,'testROI')
        % croppedVol=vol;
        % vol = double(imLoader{iChan}(timePoint));
        % testDynROIOverlay(dynROI,croppedVol,vol,rawMovieInfo,movieInfo,frameIdx,ZXRatio);
        
        debugPos(frameIdx)=debugDet;

        labelSegPos(frameIdx)=labDet;

    end %%%% end parfor (frame loop)
    if ~exist('movieInfo','var')
        %in the case that no channels/frames had detected points
        movieInfo = [];
    end


    save(outFilePaths{1,iChan}, 'movieInfo','debugPos','labelSegPos');

    if(~isempty(detP.processBuildDynROI))             
        tmp=detP.processBuildDynROI.loadFileOrCache(); % try initDynROIs
        dynROICell=tmp{1}.dynROICell;
        dynROI=dynROICell{1};
        detDynROIRef=dynROI.getDefaultRef().applyBase(Detections(movieInfo));
        [BBmin,BBmax]=dynROI.getBoundingBox(dynROI.getDefaultRef());
        detDynROIRef.addOffset(-BBmin(1)+1,-BBmin(2)+1,-BBmin(3)+1);
        movieInfoDynROIRef=detDynROIRef.getStruct();
        save(outFilePaths{2,iChan}, 'movieInfoDynROIRef');
    end
    
%     clear movieInfo detectionLabRef;
    clear movieInfo;

end %%% channel loop

% Close waitbar
if ishandle(wtBar), close(wtBar); end
disp('Finished detecting diffraction-limited objects!')
movieData.save;


function [movieInfo,feats]= labelToMovieInfo(label,vol)
    [feats,nFeats] = bwlabeln(label);
    featsProp = regionprops(feats,vol,'WeightedCentroid','MeanIntensity','MaxIntensity');

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



function testDynROIOverlay(dynROI,croppedVol,vol,rawMovieInfo,movieInfo,frameIdx,ZXRatio)
    origVolSize=size(vol);
    vol=zeros(size(vol));
    pos=ceil(size(vol).*rand(500,3));
    vol(sub2ind(size(vol),pos(:,1),pos(:,2),pos(:,3)))=1000;
    pos=pos(:,[2 1 3]);
    scaledPos=pos;
    scaledPos(:,3)=ZXRatio*scaledPos(:,3);
    
    movieInfo=Detections().initFromPosMatrices(arrayfun(@(d) scaledPos,1:numel(movieInfo),'unif',0),arrayfun(@(d) scaledPos,1:numel(movieInfo),'unif',0));

    tracks=movieInfo.buildTracksFromDetection();
    dynROI=TracksROI(tracks(1),20);

    movieInfo=dynROI.mapDetections(movieInfo);




    % Testing if dynROI volume rendering and overlay are correctly aligned
    % Visual feedback
    [Handles,~,F]=setupFigure(4,3,12);

    % Basic display using masked full volume 
    [subVol,ROI]=dynROI.cropBoundingVol(vol,ZXRatio);
    maskedVol=zeros(size(vol));
    maskedVol(:)=min(subVol(:));
    maskedVol(ROI>0)=subVol;
    maskedVol=imresize3(maskedVol,[size(vol,1) size(vol,2) ceil(ZXRatio*size(vol,3))],'nearest');

    det=Detections(movieInfo);
    s=det(frameIdx).getAllStruct();
    X=s.x;
    Y=s.y;
    Z=s.z;

    XY=mat2gray(max(maskedVol,[],3));
    ZY=mat2gray(squeeze(max(maskedVol,[],2)));
    ZX=mat2gray(squeeze(max(maskedVol,[],1)));

    H=Handles(1);
    imshow(XY, 'Parent', H);
    hold(H,'on');
    scatter(H,X,Y,100,'r','linewidth',2);
    hold(H,'off')
    H=Handles(2);
    imshow(ZY, 'Parent', H);
    hold(H,'on');
    scatter(H,Z,Y,100,'r','linewidth',2);
    hold(H,'off')   
    H=Handles(3);
    imshow(ZX, 'Parent', H);
    hold(H,'on');
    scatter(H,Z,X,100,'r','linewidth',2);
    hold(H,'off')  

     % Basic display using normal crop on
    [subVol,ROI,minCoord]=dynROI.cropBoundingVol(vol,ZXRatio);
    subVol=imresize3(subVol,[size(subVol,1) size(subVol,2) ceil(ZXRatio*size(subVol,3))]);

    det=Detections(movieInfo);
    s=det(frameIdx).getAllStruct();
    X=s.x-minCoord(1)+1;
    Y=s.y-minCoord(2)+1;
    Z=s.z-ZXRatio*(minCoord(3)-1);
    % X=[X 0 1];
    % Y=[Y 0 1];
    % Z=[Z 0 1];
    XY=mat2gray(max(subVol,[],3));
    ZY=mat2gray(squeeze(max(subVol,[],2)));
    ZX=mat2gray(squeeze(max(subVol,[],1)));

    H=Handles(4);
    imshow(XY, 'Parent', H);
    hold(H,'on');
    scatter(H,X,Y,100,'r','linewidth',2);
    hold(H,'off')
    H=Handles(5);
    imshow(ZY, 'Parent', H);
    hold(H,'on');
    scatter(H,Z,Y,100,'r','linewidth',2);
    hold(H,'off')   
    H=Handles(6);
    imshow(ZX, 'Parent', H);
    hold(H,'on');
    scatter(H,Z,X,100,'r','linewidth',2);
    hold(H,'off')  


    % Building subvolume, transforming coordinate, scaling and overlay.
    [subVol,minCoord,maxCoord] = dynROI.getSubVol(vol,ZXRatio,frameIdx);
    % MinCoord is the zeros of subVolume
    % MaxCoord is the size(subVol)
    det=Detections(movieInfo);
    ref=dynROI.getDefaultRef();
    detRef=ref.applyBase(det);

    [BBmin,BBmax]=dynROI.getBoundingBox(ref,frameIdx);
    minCoord
    maxCoord
    s=detRef(frameIdx).getAllStruct();
    disp(num2str(maxCoord-minCoord));
    disp(size(subVol));
    sX=size(subVol,2);
    sY=size(subVol,1);
    sZ=size(subVol,3);

    % s.x=[minCoord(1) maxCoord(1) 0];
    % s.y=[minCoord(2) maxCoord(2) 0];
    % s.z=[minCoord(3) maxCoord(3) 0];
    r=(sX-1)/(maxCoord(1)-minCoord(1))
    X=s.x*r+sX-r*maxCoord(1)
    r=(sY-1)/(maxCoord(2)-minCoord(2));
    Y=s.y*r+sY-r*maxCoord(2);
    r=(sZ-1)/(maxCoord(3)-minCoord(3));
    Z=s.z*r+sZ-r*maxCoord(3);
  

    XY=mat2gray(max(subVol,[],3));
    ZY=mat2gray(squeeze(max(subVol,[],2)));
    ZX=mat2gray(squeeze(max(subVol,[],1)));

    isoHandle=setupFigure(1,3,3);
    H=Handles(7);
    H=isoHandle(1);
    imshow(XY, 'Parent', H);
    hold(H,'on');
    scatter(H,X,Y,100,'r','linewidth',2);
    hold(H,'off')
    H=Handles(8);
    H=isoHandle(2);
    imshow(ZY, 'Parent', H);
    hold(H,'on');
    scatter(H,Z,Y,100,'r','linewidth',2);
    hold(H,'off')   
    H=Handles(9);
    H=isoHandle(3);
    imshow(ZX, 'Parent', H);
    hold(H,'on');
    scatter(H,Z,X,100,'r','linewidth',2);
    hold(H,'off')  

    %% cropped Vol
    croppedVol=imresize3(croppedVol,[size(croppedVol,1) size(croppedVol,2) ceil(ZXRatio*size(croppedVol,3))]);

    XY=mat2gray(max(croppedVol,[],3));
    ZY=mat2gray(squeeze(max(croppedVol,[],2)));
    ZX=mat2gray(squeeze(max(croppedVol,[],1)));
    det=Detections(rawMovieInfo);
    s=det(frameIdx).getAllStruct();
    X=s.x;
    Y=s.y;
    Z=ZXRatio*s.z;

    %% Scaling: 
    %% 1) we want (1,1,1) to be (1,1,1) on the volume
    %% 2) we want (size(croppedVolume))

   
    hoffset=9;
    H=Handles(hoffset+1);
    imshow(XY, 'Parent', H);
    hold(H,'on');
    scatter(H,X,Y,100,'r','linewidth',2);
    hold(H,'off')
    H=Handles(hoffset+2);
    imshow(ZY, 'Parent', H);
    hold(H,'on');
    scatter(H,Z,Y,100,'r','linewidth',2);
    hold(H,'off')   
    H=Handles(hoffset+3);
    imshow(ZX, 'Parent', H);
    hold(H,'on');
    scatter(H,Z,X,100,'r','linewidth',2);
    hold(H,'off')  
