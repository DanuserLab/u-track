function projectDynROI(MD,varargin)
% WIPS: 
%   -   respect computeMIPProcess Specs
%   -   using dynROI class instead of ad hoc object
%   -   comment and clarify options.
%   -   shorter function 
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
ip = inputParser;
ip.CaseSensitive = false;
ip.KeepUnmatched = true;
ip.addRequired('MD',@(x) isa(x,'MovieData')); 
ip.addOptional('dynROI',TracksROI());          % This dynROI define the boundary of the current view (default, insetDynROI) only used if 'crop' is set to manifold
ip.addOptional('insetDynROI',[],@(x) isa(x,'DynROI')||isempty(x));     % Every pixel mapped inside the ROI is rendered in a MIP mask
ip.addOptional('renderedChannel',1:length(MD.channels_),@isnumeric);
ip.addOptional('crop','manifold');  % 'manifold': crop around p.dynROI, 'full': show the full volume.
ip.addOptional('FoF',FrameOfRef().genCanonicalRef(MD.nFrames_));           % frame of refenrence for the projection
ip.addOptional('transType','affineOnePass');
ip.addOptional('processFrame',1:MD.nFrames_);
ip.addOptional('channelRender','grayRed');
ip.addOptional('intMinPrctil',[1 1]);
ip.addOptional('intMaxPrctil',[100 100]);
ip.addOptional('contrastIn',[0 1]);
ip.addOptional('contrastOut',[0 1]);
ip.addOptional('gamma',1);
ip.addOptional('insetOnly',false);
ip.addOptional('processSingleProj',[]);
ip.addOptional('suppressROIBorder',true);
ip.addOptional('format','png');
ip.addOptional('rawTIFF',false);
ip.addOptional('processMaskVolume',[])
ip.addOptional('processRenderer',[]);
ip.addOptional('saveVideo',false);
ip.addOptional('maxMIPSize',max([400,MD.imSize_,ceil(MD.zSize_*MD.pixelSizeZ_/MD.pixelSize_)]));
ip.parse(MD,varargin{:});
p=ip.Results;

dynROI=p.dynROI;
insetDynROI=p.insetDynROI;
rawTIFF=p.rawTIFF;
processFrame=p.processFrame;

%% Set normalization value for rendering (should be in a rendering function ...)
minIntensityNorm=zeros(1,numel(MD.channels_));
maxIntensityNorm=zeros(1,numel(MD.channels_));
for chIdx=1:length(MD.channels_)
    vol=MD.getChannel(chIdx).loadStack(1,'ZSize',MD.zSize_); 
    minIntensityNorm(chIdx)=[ prctile(vol(:),p.intMinPrctil(chIdx))];
    maxIntensityNorm(chIdx)=[ prctile(vol(:),p.intMaxPrctil(chIdx))];
end

fringeWidth=0;
if(~isempty(dynROI))
    fringeWidth=dynROI.fringe;
end

insetFringeWidth=0;
if(~isempty(insetDynROI))
    insetFringeWidth=insetDynROI.fringe;
end
%% Define the static Rectangular cuboid that contains the pixel to be projected in the frame of reference.
%% The coordinate of this cube are specified such as the origin of the frame of reference is the zero.

%% In the manifold crop case, the boundaries are given by the transform coordinate along the manifold polygon
%% in the full case, one have to estimate the maximum rectangle cuboid contained that can descibe the extremum coordinate of the original volume
if(strcmp(p.crop,'manifold')&&~isempty(dynROI)&&(~isempty(dynROI.tracks)))
    [minCoord,maxCoord]=dynROI.getBoundingBox(p.FoF);
else
    if(~isempty(p.FoF))
        maxCoord=[MD.getDimensions('X')-p.FoF.origin(1,1),MD.getDimensions('Y')-p.FoF.origin(1,2),MD.getDimensions('Z')*(MD.pixelSizeZ_/MD.pixelSize_)-p.FoF.origin(1,3)];
        minCoord=[1-p.FoF.origin(1,1), 1-p.FoF.origin(1,2),1-p.FoF.origin(1,3)];
    else
        maxCoord=[MD.getDimensions('X'),MD.getDimensions('Y'),MD.getDimensions('Z')*(MD.pixelSizeZ_/MD.pixelSize_)];
        minCoord=[1 1 1];
    end
end
maxXBorder=maxCoord(1);
maxYBorder=maxCoord(2);
maxZBorder=maxCoord(3);
minXBorder=minCoord(1);
minYBorder=minCoord(2);
minZBorder=minCoord(3);

tform=affine3d();

if(~isempty(p.FoF))
    B=p.FoF.getBase(1);
    tform.T(1:3,1:3)=B;
end
format=p.format;
if(rawTIFF)
    format='tif';
end

projectDynROIProcess=p.processSingleProj;
if(~isempty(projectDynROIProcess))
set(projectDynROIProcess,'ref',p.FoF);
set(projectDynROIProcess,'nFrames',length(p.processFrame));
projectDynROIProcess.setBoundingBox( ...
   [minXBorder maxXBorder],...
   [minYBorder maxYBorder],...
   [minZBorder maxZBorder] );
end

if(~isempty(p.processSingleProj))
    %% simulate run to comply with movieViewer requirement
    % Disabling because it triggers an unwanted MD save
    % p.processSingleProj.setFunName((@(x) x));
    % p.processSingleProj.run();
    %% Save the current function run for futur rerun
    p.processSingleProj.setFunName((@(p) projectDynROI(MD,varargin{:},'processSingleProj',p)));
end


format='png';
% Standardized output for processRenderer
outFilePathsRenderer = cell(1, 5);

if(~isempty(p.processRenderer))
    p.processRenderer.emptyCache();
    set(p.processRenderer,'ref',p.FoF);
    set(p.processRenderer,'nFrames',length(p.processFrame));
    p.processRenderer.setBoundingBox( ...
       [minXBorder maxXBorder],...
       [minYBorder maxYBorder],...
       [minZBorder maxZBorder] );
    %% simulate run to comply with movieViewer requirement
    % Disabling because it triggers an unwanted MD save
%     p.processRenderer.setFunName((@(x) x));
%     p.processRenderer.run();
    %% Save the current function run for futur rerun
    p.processRenderer.setFunName((@(p) projectDynROI(MD,varargin{:},'processRenderer',p)));
end


% Optional process for saving sparse volume
if(~isempty(p.processMaskVolume))
    outFilePathsMaskVolume=cell(2,numel(MD.channels_));
    frameNb=length(p.processFrame);
    for i = 1:numel(MD.channels_);    
        outFilePathsMaskVolume{1,i} = [outputDirSingleProj filesep 'volMask'  filesep 'ch' filesep 'frame_nb%04d.mat'];
        outFilePathsMaskVolume{2,1} = [outputDirSingleProj filesep 'volMask' filesep 'limits.mat'];
    end;
    for mIdx=1:numel(p.processFrame)
        mkdirRobust(fileparts(outFilePathsMaskVolume{1,i}));
    end
    save(outFilePathsMaskVolume{2,1},'minXBorder', 'maxXBorder','minYBorder','maxYBorder','minZBorder','maxZBorder','frameNb');
    p.processMaskVolume.setOutFilePaths(outFilePathsMaskVolume);
end

% if(~isempty(insetDynROI))
%     insetDynROI=insetDynROI.tracks;
% end

% warp, crop, fuse and save each time point
% <<<<<<< HEAD
% parfor fIdx = processFrame
% =======
% vol = cell(processFrame,1);

%% Allowing parfor again
renderingCellXY=cell(1,numel(processFrame));
renderingCellZY=cell(1,numel(processFrame));
renderingCellZX=cell(1,numel(processFrame));

parfor fIdx = processFrame
    fprintf('.') 
    % produce a ROI mask using the 1D polygon (segment defined by the extremities of the insetDynROI).
    % todo: N Channel (now 2).
    mask=ones(MD.imSize_(1),MD.imSize_(2),MD.zSize_);
    maskedVol=[];


    if(~isempty(insetDynROI))
        if(isa(insetDynROI,'DynROI'))
            ZRatio=MD.pixelSizeZ_/MD.pixelSize_;
            mask=insetDynROI.getMask([MD.imSize_(1),MD.imSize_(2),MD.zSize_],ZRatio,fIdx);
            mask=mask{fIdx};
        else
            % Collect relative frameIdx
            pIndices=nan(1,length(insetDynROI));
            for polIdx=1:length(insetDynROI)
                F=insetDynROI(polIdx).f;
                pIdx=find(F==fIdx,1);
                if isempty(pIdx)
                    if(fIdx>max(F))   pIdx=length(F);  else   pIdx=1; 
                    end;
                end
                pIndices(polIdx)=pIdx;
            end

            % Building mask in the 1D case
            nextPoint=length(insetDynROI);
            PCurrent=[insetDynROI(1).x(pIndices(1)) insetDynROI(1).y(pIndices(1)) insetDynROI(1).z(pIndices(1))];
            KCurrent=[insetDynROI(nextPoint).x(pIndices(nextPoint)) insetDynROI(nextPoint).y(pIndices(nextPoint)) insetDynROI(nextPoint).z(pIndices(nextPoint))];

            vol = MD.getChannel(1).loadStack(fIdx,'ZSize',MD.zSize_);
            disp(num2str([size(vol) MD.imSize_ MD.zSize_]));

        % Building mask for both channel on the whole volume
        % NOTE: in order to apply fringe isotropically, we need the mask to
        % be isotropized briefly.
        mask=zeros(size(vol,1),size(vol,2),ceil(size(vol,3)*MD.pixelSizeZ_/MD.pixelSize_));
        sampling=100;
        xSeg=round(linspace(PCurrent(1),KCurrent(1),sampling));
        ySeg=round(linspace(PCurrent(2),KCurrent(2),sampling));
        zSeg=round(linspace(PCurrent(3),KCurrent(3),sampling));
        indx=sub2ind(size(mask),ySeg,xSeg,zSeg);
        
        mask(indx)=1;
        
        %mask=imdilate(mask,IMSphere);  %ones(cubeHalfWidth,cubeHalfWidth,round(cubeHalfWidth*MD.pixelSize_/MD.pixelSizeZ_)));
        % If no transform are needed, now to save on bwdist.

        distMap=mask;
        distMap=bwdist(distMap);
        mask(distMap<(insetFringeWidth+1))=1;
        [y x z]=...
        ndgrid( linspace(1,size(mask,1),size(vol,1)),...
            linspace(1,size(mask,2),size(vol,2)),...
            linspace(1,size(mask,3),size(vol,3)));
        mask=interp3(mask,x,y,z);
        end
    end
    mips=cell(3,length(MD.channels_));
    if(~isempty(insetDynROI))
        if(isempty(p.FoF))&&(strcmp(p.crop,'manifold'))
            aminXBorder=max(1,minXBorder);
            aminYBorder=max(1,minYBorder);
            aminZBorder=max(1,minZBorder);
            amaxXBorder=min(size(mask,2),maxXBorder);
            amaxYBorder=min(size(mask,1),maxYBorder);
            amaxZBorder=min(size(mask,3)*MD.pixelSizeZ_/MD.pixelSize_,maxZBorder);
            
            XCropMask=false(1,size(mask,2));
            XCropMask(aminXBorder:amaxXBorder)=true;
            
            YCropMask=false(1,size(mask,1));
            YCropMask(aminYBorder:amaxYBorder)=true;
            
            ZCropMask=false(1,size(mask,3));
            ZCropMask(ceil((aminZBorder:amaxZBorder)*MD.pixelSize_/MD.pixelSizeZ_))=true;
            
            mask(:,:,~ZCropMask)=[];
            mask(~YCropMask,:,:)=[];
            mask(:,~XCropMask,:)=[];
        end
    end
    
    for chIdx=1:length(MD.channels_)
        vol=MD.getChannel(chIdx).loadStack(fIdx,'ZSize',MD.zSize_);
        if(~isempty(dynROI))
            if(isempty(p.FoF))&&(strcmp(p.crop,'manifold'))
                aminXBorder=max(1,minXBorder);
                aminYBorder=max(1,minYBorder);
                aminZBorder=max(1,minZBorder);
                amaxXBorder=min(size(vol,2),maxXBorder);
                amaxYBorder=min(size(vol,1),maxYBorder);
                amaxZBorder=min(size(vol,3)*MD.pixelSizeZ_/MD.pixelSize_,maxZBorder);
                
                XCropVol=false(1,size(vol,2));
                XCropVol(aminXBorder:amaxXBorder)=true;
                
                YCropVol=false(1,size(vol,1));
                YCropVol(aminYBorder:amaxYBorder)=true;
                
                ZCropVol=false(1,size(vol,3));
                ZCropVol(ceil((aminZBorder:amaxZBorder)*MD.pixelSize_/MD.pixelSizeZ_))=true;

                vol(:,:,~ZCropVol)=[];
                vol(~YCropVol,:,:)=[];
                vol(:,~XCropVol,:)=[];
                
            end
            
            maskedVol=vol;
            maskedVol(~mask)=0;
        end
        
        % If needed the map must rotated before cropped (efficiency)
        % Rotation depends on the FrameOfRef associated to the tracks the compose the dynanimc polygon
        % Cropping area according to the polygon OVER TIME plus added vizualiation margin
        % Rotation will use imwarp
        % Can we use imwar for cropping too ?
        
        %% if a FoF is specified, warp and crop data according to the
        tform=affine3d();

        warpedVol=vol;
        warpedMaskedVol=[];
        warpedMask=mask;
        if(~isempty(insetDynROI))
            warpedMaskedVol=maskedVol;
        else
            warpedMaskedVol=zeros(size(warpedVol));
        end
        if(~isempty(p.FoF))
            B=p.FoF.getBase(fIdx);
            tform.T(4,[1 2 3])=(-p.FoF.getOrigAtFrame(fIdx)+p.FoF.origin(1,:))*B;
            tform.T(1:3,1:3)=B;
            %
            tformTransOnly=affine3d();
            tformTransOnly.T(4,[1 2 3])=(-p.FoF.getOrigAtFrame(fIdx));
            
            %
            tformRelTransOnly=affine3d();
            tformRelTransOnly.T(4,[1 2 3])=(-p.FoF.origin(1,:)+p.FoF.getOrigAtFrame(fIdx));
            
            tformRotOnly=affine3d();
            B=p.FoF.getBase(fIdx);
            tformRotOnly.T(1:3,1:3)=B(:,[1 2 3]);
            
            tformRotOnlyInit=affine3d();
            B=p.FoF.getBase(1);
            tformRotOnlyInit.T(1:3,1:3)=B;
            
            orig=p.FoF.getOrigAtFrame(fIdx);
            
            inputRef=imref3d([ MD.getDimensions('Y') MD.getDimensions('X') MD.getDimensions('Z')], ...
                [1 MD.getDimensions('X')],[1 MD.getDimensions('Y')],[1 MD.getDimensions('Z')*MD.pixelSizeZ_/MD.pixelSize_]);
            
            
            switch p.transType
                case 'affineOnePass'
                    maxXBorderFull=MD.getDimensions('X');
                    maxYBorderFull=MD.getDimensions('Y');
                    maxZBorderFull=MD.getDimensions('Z')*(MD.pixelSizeZ_/MD.pixelSize_);
                    minXBorderFull=1;
                    minYBorderFull=1;
                    minZBorderFull=1;
                    
                    minXBorderCurr=minXBorderFull - orig(1); maxXBorderCurr=maxXBorderFull - orig(1);
                    minYBorderCurr=minYBorderFull - orig(2); maxYBorderCurr=maxYBorderFull - orig(2);
                    minZBorderCurr=minZBorderFull - orig(3); maxZBorderCurr=maxZBorderFull - orig(3);
                    
                    inputRef=imref3d([ MD.getDimensions('Y') MD.getDimensions('X') MD.getDimensions('Z')], ...
                        [minXBorderCurr maxXBorderCurr], [minYBorderCurr maxYBorderCurr], [minZBorderCurr maxZBorderCurr]);
                    
                    minXBorderCurr=minXBorder ;
                    maxXBorderCurr=maxXBorder ; 
                    minYBorderCurr=minYBorder ;
                    maxYBorderCurr=maxYBorder ;
                    minZBorderCurr=minZBorder ;
                    maxZBorderCurr=maxZBorder ;
                    
                    rotOutputRef=imref3d([    ceil(maxYBorderCurr-minYBorderCurr) ...
                        ceil(maxXBorderCurr-minXBorderCurr) ...
                        ceil(maxZBorderCurr-minZBorderCurr) ], ...
                        [minXBorderCurr maxXBorderCurr], [minYBorderCurr maxYBorderCurr], [minZBorderCurr maxZBorderCurr]);
                    
                    if(~p.insetOnly)
                            warpedVol=imwarp(vol,inputRef,tformRotOnly,'OutputView',rotOutputRef);
                    end
                    if(~isempty(insetDynROI))
                        if(p.suppressROIBorder&&(~p.insetOnly))
                            warpedMask=imwarp(double(mask),inputRef,tformRotOnly,'OutputView',rotOutputRef);
                            warpedMaskedVol=warpedMask;
                            warpedMaskedVol(warpedMask>0)=warpedVol(warpedMask>0);
                        else
                            warpedMaskedVol=imwarp(maskedVol,inputRef,tformRotOnly,'OutputView',rotOutputRef);
                        end
                    else
                        if(p.insetOnly)
                            error('No inset specified for <p.insetOnly>');
                        else
                            warpedMaskedVol=zeros(size(warpedVol));
                        end
                    end
                    
                case 'translation'
                    disp(num2str(fIdx))
                    minXBorderCurr=minXBorder ;%+ orig(1) - p.FoF.origin(1,1);
                    maxXBorderCurr=maxXBorder ;%+ orig(1) - p.FoF.origin(1,1);
                    minYBorderCurr=minYBorder ;%+ orig(2) - p.FoF.origin(1,2);
                    maxYBorderCurr=maxYBorder ;%+ orig(2) - p.FoF.origin(1,2);
                    minZBorderCurr=minZBorder ;%+ orig(3) - p.FoF.origin(1,3);
                    maxZBorderCurr=maxZBorder ;%+ orig(3) - p.FoF.origin(1,3);
                    
                    %             [xLimitsOut,yLimitsOut,zLimitsOut] = outputLimits(tformTransOnly,[minXBorder maxXBorder], [minYBorder maxYBorder], [minZBorder maxZBorder]);
                    %             minXBorderCurr=xLimitsOut(1); maxXBorderCurr=xLimitsOut(2);
                    %             minYBorderCurr=yLimitsOut(1); maxYBorderCurr=yLimitsOut(2);
                    %             minZBorderCurr=zLimitsOut(1); maxZBorderCurr=zLimitsOut(2);
                    %
                    outputRef=imref3d([ ceil(maxYBorderCurr-minYBorderCurr) ...
                        ceil(maxXBorderCurr-minXBorderCurr) ...
                        ceil(maxZBorderCurr-minZBorderCurr) ], ...
                        [minXBorderCurr maxXBorderCurr], [minYBorderCurr maxYBorderCurr], [minZBorderCurr maxZBorderCurr]);
                    
                    warpedVol=imwarp(vol,inputRef,tformTransOnly,'OutputView',outputRef);
                    warpedMaskedVol=imwarp(maskedVol,inputRef,tformTransOnly,'OutputView',outputRef);
                    
                case 'transCrop'
                    %% to be updated
                    [xLimitsOut,yLimitsOut,zLimitsOut] = outputLimits(tformRelTransOnly,[minXBorder maxXBorder], [minYBorder maxYBorder], [minZBorder maxZBorder]);
                    minXBorderCurr=xLimitsOut(1); maxXBorderCurr=xLimitsOut(2);
                    minYBorderCurr=yLimitsOut(1); maxYBorderCurr=yLimitsOut(2);
                    minZBorderCurr=zLimitsOut(1); maxZBorderCurr=zLimitsOut(2);
                    
                    maskcrop=maskedVol;
                    nullMaskXY=(squeeze(any(maskcrop,3)));
                    YNull=~(squeeze(any(any(mask,3),2)));
                    XNull=~(squeeze(any(any(mask,3),1)));
                    ZNull=~(squeeze(any(any(mask,1),2)));
                    
                    YNull= zeros(1,size(maskcrop,1));
                    YNull(1:minYBorderCurr)=1;
                    YNull(maxYBorderCurr:end)=1;
                    YNull=logical(YNull);
                    
                    XNull= zeros(1,size(maskcrop,2));
                    XNull(1:minXBorderCurr)=1;
                    XNull(maxXBorderCurr:end)=1;
                    XNull=logical(XNull);
                    
                    ZNull= zeros(1,size(maskcrop,3));
                    ZNull(1:ceil(minZBorderCurr*MD.pixelSize_/MD.pixelSizeZ_))=1;
                    ZNull(ceil(maxZBorderCurr*MD.pixelSize_/MD.pixelSizeZ_):end)=1;
                    ZNull=logical(ZNull);
                    
                    maskcrop(:,:,ZNull)=[];
                    maskcrop(YNull,:,:)=[];
                    maskcrop(:,XNull,:)=[];
                    
                    warpedMaskedVol=maskcrop;
                    
                    warpedVol=vol;
                    warpedVol(:,:,ZNull)=[];
                    warpedVol(YNull,:,:)=[];
                    warpedVol(:,XNull,:)=[];
                    
                otherwise
                    error('unknown trans type');
            end
        end

        
        %% Create MIPS for each channel, fuse mask and full volume
        ZRatio=1;
        switch p.transType
            case 'none'
            case 'transCrop'
                ZRatio=MD.pixelSizeZ_/MD.pixelSize_;
            otherwise
                ZRatio=1;
        end;
        if(isempty(p.FoF))
            ZRatio=MD.pixelSizeZ_/MD.pixelSize_;
        end

        [maxXY,maxZY,maxZX,~]=computeMIPs(warpedMaskedVol,ZRatio, ...
                minIntensityNorm(chIdx),maxIntensityNorm(chIdx),'raw',rawTIFF);
        if(p.insetOnly)
            fullmaxXY=zeros(size(maxXY));
            fullmaxZY=zeros(size(maxZY));
            fullmaxZX=zeros(size(maxZX));
        else
            [fullmaxXY,fullmaxZY,fullmaxZX,~]=computeMIPs(warpedVol,ZRatio, ...
                minIntensityNorm(chIdx),maxIntensityNorm(chIdx),'raw',rawTIFF);
        end
            
        [maskXY,maskZY,maskZX,~]=computeMIPs(warpedMask,ZRatio, ...
                minIntensityNorm(chIdx),maxIntensityNorm(chIdx),'raw',true);

%         imdisp({fullmaxZY,maxZY,maskZY});
        % disp(any(maxXY(:)))
        % disp(all(maskXY(:)))
        
        % Create MIP of ROI and context
        if(~isempty(insetDynROI))
            if(~p.insetOnly)
            if(p.suppressROIBorder)
            maxXY(imerode(maskXY,ones(3))<1)=fullmaxXY(imerode(maskXY,ones(3))<1);
            maxZY(imerode(maskZY,ones(3))<1)=fullmaxZY(imerode(maskZY,ones(3))<1);
            maxZX(imerode(maskZX,ones(3))<1)=fullmaxZX(imerode(maskZX,ones(3))<1);   
            else
            % Create MIP of ROI and context
                maxXY(maskXY==0)=fullmaxXY(maskXY==0);
                maxZY(maskXY==0)=fullmaxZY(maskXY==0);
                maxZX(maskXY==0)=fullmaxZX(maskXY==0);                        
            end
            end
        else
            maxXY=fullmaxXY;
            maxZY=fullmaxZY;
            maxZX=fullmaxZX;                        
        end
        
        %% Resize and fuse channel MIPS
        maxMIPSize=p.maxMIPSize;
        [sX,sY,sZ]=size(warpedMaskedVol);
         sZ=sZ*ZRatio;
        resizeScale=maxMIPSize/max([sX,sY,sZ]);
        
        XYMax=imresize(maxXY,resizeScale,'nearest');
        ZYMax=imresize(maxZY,resizeScale,'nearest');
        ZXMax=imresize(maxZX,resizeScale,'nearest');
        
        if(rawTIFF)
            mips{1,chIdx} = uint8((2^8-1)*mat2gray(XYMax,double([minIntensityNorm(chIdx),maxIntensityNorm(chIdx)])));
            mips{2,chIdx} = uint8((2^8-1)*mat2gray(ZYMax,double([minIntensityNorm(chIdx),maxIntensityNorm(chIdx)])));
            mips{3,chIdx} = uint8((2^8-1)*mat2gray(ZXMax,double([minIntensityNorm(chIdx),maxIntensityNorm(chIdx)])));
        else
            mips{1,chIdx} = XYMax;
            mips{2,chIdx} = ZYMax;
            mips{3,chIdx} = ZXMax;
        end
            
       
        if(~isempty(p.processSingleProj))
            projectDynROIProcess.saveFrame(chIdx,fIdx,maxXY,maxZY,maxZX);
        end


        %% save sparse Mask volume
        if(~isempty(p.processMaskVolume))
            sparseMask=ndSparse(double(warpedMaskedVol));
            saveMask(sprintfPath(p.processMaskVolume.outFilePaths_{1,chIdx},fIdx),sparseMask)
        end
    end
    
    if(~isempty(p.processRenderer))
       gamma=repmat(p.gamma,[1 size(mips,2)]);
       for chIdx=1:size(mips,2)
           for pIdx=1:size(mips,1)
                mips{pIdx,chIdx}=imadjust(mips{pIdx,chIdx},p.contrastIn,p.contrastOut,gamma(chIdx));
           end
       end
           
        
        %% fuse volume if 2 channels
        if(length(p.renderedChannel)==2)
            XYProj=renderChannel(mips{1,1},mips{1,2},p.channelRender);
            ZYProj=renderChannel(mips{2,1},mips{2,2},p.channelRender);
            ZXProj=renderChannel(mips{3,1},mips{3,2},p.channelRender);
        else
            XYProj=repmat(mips{1,p.renderedChannel(1)},1,1,3);
            ZYProj=repmat(mips{2,p.renderedChannel(1)},1,1,3);
            ZXProj=repmat(mips{3,p.renderedChannel(1)},1,1,3);
        end
        

        renderingCellXY{fIdx}=XYProj;
        renderingCellZY{fIdx}=ZYProj;
        renderingCellZX{fIdx}=ZXProj;
    end
end


%% write images
for fIdx=processFrame
    if(~isempty(p.processRenderer))
        p.processRenderer.saveFrame(1,fIdx, renderingCellXY{fIdx}, renderingCellZY{fIdx}, renderingCellZX{fIdx});
    end
end

if((~isempty(p.processRenderer))&&(p.saveVideo)) 
    ProjAnimation(p.processRenderer,'ortho').saveVideo([p.processRenderer.getOutputDir()  '.avi']);
end


    
    

function saveMask(filepath,sparseMask)
    save(filepath,'sparseMask');