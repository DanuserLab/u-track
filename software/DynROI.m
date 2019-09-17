classdef DynROI < hgsetget & matlab.mixin.Copyable & handle

    properties (SetAccess = public, GetAccess = public)
        defaultRef; % If the ROI if 1D or 2D, it will generate a default frame of reference. 
    end
    methods
    function obj = setDefaultRef(obj,ref)
    	obj.defaultRef=ref;
    end

    function ref = getDefaultRef(obj)
    	ref=obj.defaultRef;
    end

    function [det,indices,distances] = mapDetections(obj,detections)
        [indices,distances]=obj.mapPosition(detections.getPosMatrix());
        det=detections.getSelectIdx(indices);
    end

    function [maskCell] = getMask(obj,volSize,ZRatio,frames)

        %% In a first version, the mask is built using the mapPosition function.
        %% Every pixel in the boundingbox is transformed as a position 
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
        maskCell=cell(1,numel(frames));
        pos=cell(1,max(frames));
        for fIdx=frames
            [m,M]=obj.getBoundingBox();
            [y x z]=ndgrid(m(2):M(2),m(1):M(1),m(3):M(3));
            posMat=[x(:)  y(:) z(:)];
            pos{fIdx}=posMat;
        end

        %% every position in the dynROI is mapped
        [indices]=obj.mapPosition(pos);

        %% Every mapped position is a pixel in a mask of size <volSize>
        for fIdx=frames   
            mask=false(volSize);
            posMat=pos{fIdx};
            PosZResize=min(max(1,round(posMat(indices{fIdx},3)/ZRatio)),volSize(3));
            mask(sub2ind(volSize,round(posMat(indices{fIdx},2)),round(posMat(indices{fIdx},1)),PosZResize))=true;
            maskCell{fIdx}=mask;
        end
    end

    function [XY,ZY,ZX,minCoord,maxCoord]=getMIP(obj,MD,chIdx,processFrames)
        XY=cell(1,numel(processFrames));
        ZY=cell(1,numel(processFrames));
        ZX=cell(1,numel(processFrames));
        ZXRatio=MD.pixelSizeZ_/MD.pixelSize_;

        stackLoader = {@(f)(MD.channels_(chIdx).loadStack(f))};

        ref=obj.getDefaultRef();
        if(~isempty(ref))
            [minCoord,maxCoord]=obj.getBoundingBox(ref);
            minCoord=floor(minCoord);
            maxCoord=ceil(maxCoord);
        else
            minCoord=[1 1 1];
            maxCoord=[MD.imSize_(2) MD.imSize_(1) MD.zSize_*ZXRatio];
        end

        parfor f=processFrames
            vol=stackLoader{1}(f);
            % For added robustness, works if DynROI is a dummy
            if(~isempty(obj.getDefaultRef()))
                subVol=obj.getSubVol(vol,ZXRatio,f);
            else
                subVol=imresize3(vol,[size(vol,1) size(vol,2) ZXRatio*size(vol,3)]);
            end
            XY{f} = max(subVol,[],3);
            ZY{f} = (squeeze(max(subVol, [], 2)));
            ZX{f} = (squeeze(max(subVol, [], 1)));
        end
    end

    function [MDout,MDFile]=swapRawBoundingBox(obj,MD,subFolder)
        % Swap the box bounding the dynROI in the lab ref. 
        XY=cell(1,MD.nFrames_);
        ZXRatio=MD.pixelSizeZ_/MD.pixelSize_;

        outputDir=[MD.outputDirectory_ filesep subFolder filesep];
        channelList(numel(MD.channels_))=Channel();

        for cIdx=1:numel(MD.channels_)
            stackLoader = arrayfun(@(f) @(f)(MD.channels_(cIdx).loadStack(f)),1:(MD.nFrames_),'unif',0);
            chFolder=[outputDir filesep 'ch' num2str(cIdx) filesep];
            mkClrDir(chFolder);
            filePattern=[chFolder filesep 'subvol-ch-' num2str(cIdx) '-frame-%04d.tif'];

            parfor f=1:MD.nFrames_
                disp(['swap sub volume ' num2str(f)])
                vol=stackLoader{f}(f);

                [BBmin,BBmax]=obj.getBoundingBox();
                maskMinX=BBmin(2); maskMinY=BBmin(1); maskMinZ=ceil(BBmin(3)/ZXRatio); 
                maskMaxX=BBmax(2); maskMaxY=BBmax(1); maskMaxZ=floor(BBmax(3)/ZXRatio); 
                ROI=false(size(vol));
                ROI(max(maskMinX,1):min(maskMaxX,end),max(1,maskMinY):min(end,maskMaxY),max(maskMinZ,1):min(end,maskMaxZ))=true;
                [maskMinX,maskMinY,maskMinZ]=ind2sub(size(ROI), find(ROI,1));
                [maskMaxX,maskMaxY,maskMaxZ]=ind2sub(size(ROI), find(ROI,1,'last'));
                tmp = nan(1+[maskMaxX,maskMaxY,maskMaxZ]-[maskMinX,maskMinY,maskMinZ]);
                tmp(:) = vol(ROI);
                origVol=vol;
                vol = tmp;

                stackWrite(uint16(vol),sprintfPath(filePattern,f));
           end
           channelList(cIdx)=Channel(chFolder);
        end

        tiffReader=TiffSeriesReader({channelList.channelPath_},'force3D',true);
        MDout=MovieData(channelList,[outputDir 'analysis'],'movieDataFileName_','movieData.mat','movieDataPath_',[outputDir filesep 'analysis'], ...
            'pixelSize_',MD.pixelSize_,'pixelSizeZ_',MD.pixelSizeZ_,'timeInterval_',MD.timeInterval_);
        MDout.setReader(tiffReader);
        MDout.sanityCheck(); % the usual fucking performance killer...
        MDout.save();
        MDFile=[outputDir filesep 'analysis' filesep 'movieData.mat'];
    end

    function [MDout,MDFile]=swapDynBoundingBox(obj,MD,subFolder)
        % Swap the box bounding the dynROI in the default ref. 
        ZXRatio=MD.pixelSizeZ_/MD.pixelSize_;

        outputDir=[MD.outputDirectory_ filesep subFolder filesep];
        channelList(numel(MD.channels_))=Channel();

        for cIdx=1:numel(MD.channels_)
            stackLoader = arrayfun(@(f) @(f)(MD.channels_(cIdx).loadStack(f)),1:(MD.nFrames_),'unif',0);
            chFolder=[outputDir filesep 'ch' num2str(cIdx) filesep];
            mkClrDir(chFolder);
            filePattern=[chFolder filesep 'subvol-ch-' num2str(cIdx) '-frame-%04d.tif'];

            parfor f=1:MD.nFrames_
                disp(['swap sub volume ' num2str(f)])
                vol=stackLoader{f}(f);
                vol=obj.getSubVol(vol,ZXRatio,f);
                stackWrite(uint16(vol),sprintfPath(filePattern,f));
           end
           channelList(cIdx)=Channel(chFolder);
        end

        tiffReader=TiffSeriesReader({channelList.channelPath_},'force3D',true);
        MDout=MovieData(channelList,[outputDir 'analysis'],'movieDataFileName_','movieData.mat','movieDataPath_',[outputDir filesep 'analysis'], ...
            'pixelSize_',MD.pixelSize_,'pixelSizeZ_',MD.pixelSize_,'timeInterval_',MD.timeInterval_);
        MDout.setReader(tiffReader);
        MDout.sanityCheck(); % the usual fucking performance killer...
        MDout.save();
        MDFile=[outputDir filesep 'analysis' filesep 'movieData.mat'];
    end
    
    function [subVol,minCoord,maxCoord] = getSubVol(obj,vol,ZXRatio,frameIdx)
        % The intrinsic coordinate values (x,y,z) of the center point of any
        % pixel are identical to the values of the column, row, and plane
        % subscripts for that pixel. For example, the center point of the
        % pixel in row 5, column 3, plane 4 has intrinsic coordinates x = 3.0,
        % y = 5.0, z = 4.0.

        % The order of the coordinate specification (3.0,5.0,4.0) is reversed
        % in intrinsic coordinates relative to pixel subscripts (5,3,4).
        % Intrinsic coordinates are defined on a continuous plane, while the
        % subscript locations are discrete locations with integer values.
        ref=obj.getDefaultRef();

        maxXBorderFull=size(vol,2);
        maxYBorderFull=size(vol,1);
        maxZBorderFull=size(vol,3)*ZXRatio;
        minXBorderFull=1;
        minYBorderFull=1;
        minZBorderFull=ZXRatio;
        
        orig=ref.getOrigAtFrame(frameIdx);

        minXBorderCurr=minXBorderFull - orig(1); maxXBorderCurr=maxXBorderFull - orig(1);
        minYBorderCurr=minYBorderFull - orig(2); maxYBorderCurr=maxYBorderFull - orig(2);
        minZBorderCurr=minZBorderFull - orig(3); maxZBorderCurr=maxZBorderFull - orig(3);
        
        inputRef=imref3d(size(vol), ...
            [minXBorderCurr maxXBorderCurr], [minYBorderCurr maxYBorderCurr], [minZBorderCurr maxZBorderCurr]);
        
        [minCoord,maxCoord]=obj.getBoundingBox(ref);
        maxXBorder=maxCoord(1);
        maxYBorder=maxCoord(2);
        maxZBorder=maxCoord(3);
        minXBorder=minCoord(1);
        minYBorder=minCoord(2);
        minZBorder=minCoord(3);
        
        minXBorderCurr=(minXBorder);
        maxXBorderCurr=(maxXBorder); 
        minYBorderCurr=(minYBorder);
        maxYBorderCurr=(maxYBorder);
        minZBorderCurr=(minZBorder);
        maxZBorderCurr=(maxZBorder);
        
        rotOutputRef=imref3d([    ceil(maxYBorderCurr-minYBorderCurr)+1 ...
            ceil(maxXBorderCurr-minXBorderCurr)+1 ...
            ceil(maxZBorderCurr-minZBorderCurr)+1 ], ...
             [minXBorderCurr maxXBorderCurr], [minYBorderCurr maxYBorderCurr], [minZBorderCurr maxZBorderCurr]);
       
        B=ref.getBase(frameIdx);
        tformRotOnly=affine3d();
        tformRotOnly.T(1:3,1:3)=B(:,[1 2 3]);
        subVol=imwarp(vol,inputRef,tformRotOnly,'OutputView',rotOutputRef);
        minCoord=[minXBorderCurr minYBorderCurr minZBorderCurr];
        maxCoord=[maxXBorderCurr maxYBorderCurr maxZBorderCurr];
    end

    function [subVol,ROI,minCoord] = cropBoundingVol(obj,vol,ZXRatio)
        [BBmin,BBmax]=obj.getBoundingBox();
        BBmin=floor(BBmin);
        BBmax=ceil(BBmax);
        maskMinX=BBmin(1); maskMinY=BBmin(2); maskMinZ=ceil(BBmin(3)/ZXRatio); 
        maskMaxX=BBmax(1); maskMaxY=BBmax(2); maskMaxZ=floor(BBmax(3)/ZXRatio); 
        ROI=false(size(vol));
        ROI(max(maskMinY,1):min(maskMaxY,end),max(1,maskMinX):min(end,maskMaxX),max(maskMinZ,1):min(end,maskMaxZ))=true;
        [maskMinY,maskMaxX,maskMinZ]=ind2sub(size(ROI), find(ROI,1));
        [maskMaxY,maskMaxX,maskMaxZ]=ind2sub(size(ROI), find(ROI,1,'last'));
        tmp = nan(1+[maskMaxY,maskMaxX,maskMaxZ]-[maskMinY,maskMinX,maskMinZ]);
        tmp(:) = vol(ROI);
        subVol = tmp;
        minCoord=[maskMinX,maskMinY,maskMinZ];
    end

    function [vertices,edges]=getBoundingParallelogram(obj,frames)
        [m,M]=obj.getBoundingBoxOptim(obj.getDefaultRef(),frames);
        
        vertices=zeros(8,3);
        vertices(1,:)=m;
        vertices(2,:)=[m(1) m(2) M(3)];
        vertices(3,:)=[m(1) M(2) M(3)];
        vertices(4,:)=[m(1) M(2) m(3)];
        vertices(5,:)=[M(1) m(2) m(3)];
        vertices(6,:)=[M(1) m(2) M(3)];
        vertices(7,:)=[M(1) M(2) M(3)];
        vertices(8,:)=[M(1) M(2) m(3)];

        edges=zeros(12,2);
        edges(1,:)=[1 2];
        edges(2,:)=[2 3];
        edges(3,:)=[3 4];
        edges(4,:)=[1 4];

        edges(5,:)=[5 6];
        edges(6,:)=[6 7];
        edges(7,:)=[7 8];
        edges(8,:)=[8 5];

        edges(9,:)=[1 5];
        edges(10,:)=[2 6];
        edges(11,:)=[3 7];
        edges(12,:)=[4 8];
        
        vertices=obj.getDefaultRef().applyInvBaseToPosPointCloud(vertices,frames);
	end

    function displayBox(obj,varargin)
        ip = inputParser;
        ip.CaseSensitive = false;
        ip.KeepUnmatched=true;
        ip.addOptional('label',{},(@(x) iscell(x)||isnumeric(x)));
        ip.addOptional('handle',[]);
        ip.addOptional('frames',[]);
        ip.addOptional('MarkerSize',50);
        ip.addOptional('timeInterval',0.1);
        ip.addOptional('color',[1 0 0])
        ip.parse(varargin{:});
        p=ip.Results;

        if(isempty(p.handle))
            figure();
            h=gca();
        else
            h=p.handle;
            hold on;
        end

        frames=p.frames;
        if(isempty(frames))
            frames=obj.getDefaultRef().frame;
        end
        % xlim([1 500]);
        % ylim([1 500]);
        % zlim([1 500]);
        [vertices,edges]=obj.getBoundingParallelogram(min(frames));
        sH=scatter3(h,vertices(:,1),vertices(:,2),vertices(:,3),1,'MarkerFaceColor',[.5 0.1 0.1]);
        hold on 
        pH=arrayfun(@(e) plot3(h,vertices(edges(e,:),1),vertices(edges(e,:),2),vertices(edges(e,:),3),'Color', p.color),1:size(edges,1),'unif',0)
        

        % for f=frames 
        %     [vertices,edges]=obj.getBoundingParallelogram(f);

        %     sH.XData = vertices(:,1); 
        %     sH.YData = vertices(:,2); 
        %     sH.ZData=  vertices(:,3);

        %     for eIdx=1:numel(pH)
        %         pH{eIdx}.XData=vertices(edges(eIdx,:),1);
        %         pH{eIdx}.YData=vertices(edges(eIdx,:),2);
        %         pH{eIdx}.ZData=vertices(edges(eIdx,:),3);
        %     end

        %      % pause 2/10 second: 
        %      pause(p.timeInterval);
        %  end
     end
     
    end
    methods (Abstract)
        [minCoord, maxCoord]=getBoundingBox(obj,ref,frameIdx);
        frame=getStartFrame(obj);
        frame=getEndFrame(obj);
        [indices,dists]=mapPosition(obj,positionCell);
        [tracks,indices,distances] = mapTracks(obj,detections);
           

    end
end
