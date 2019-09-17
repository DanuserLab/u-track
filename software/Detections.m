classdef Detections <  handle  & matlab.mixin.Copyable & dynamicprops
    % Data encapsulator for 3D detections
    % General todo: better model.
    %
    % June 2018, Mark Kittisopikul implemented minimal 2D mode
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
    properties %(SetAccess = protected)
       xCoord;    % 
       yCoord;    % 
       zCoord;
       amp;       % Image amplitude (when estimated) 
       labels;    % MxN matric for N labels 
       ref;       % optional FrameOfRef objects. 
    end
    
    methods
        
        function obj = Detections(movieInfo)
            if (nargin>0)
                if(isfield(movieInfo,'zCoord')||any(isprop(movieInfo,'zCoord')))
                    for i=1:length(movieInfo)
                        obj(i).xCoord=movieInfo(i).xCoord;
                        obj(i).yCoord=movieInfo(i).yCoord;
                        obj(i).zCoord=movieInfo(i).zCoord;
                        obj(i).amp=movieInfo(i).amp;
                    end
                else
                    for i=1:length(movieInfo)
                        obj(i).xCoord=movieInfo(i).xCoord;
                        obj(i).yCoord=movieInfo(i).yCoord;
                        obj(i).zCoord=[];
                        obj(i).amp=movieInfo(i).amp;
                    end
                end
            else
                obj.zCoord=[];
                obj.xCoord=[];
                obj.yCoord=[];
            end
        end

        function setLabel(obj,labelName,values)
            if(~iscell(values)) values={values}; end;
            for i=1:numel(values)
                obj(i).labels=setfield(obj(i).labels,labelName,values{i});
,           end
        end

        function values=getLabel(obj,labelName)
            values=cell(1,numel(obj));
            for i=1:numel(values)   
                values{i}=getfield(obj(i).labels,labelName);
            end
         end

         function proc=saveInProcess(obj,MD,processName,channelIdx)
             outputFolder=fullfile(MD.outputDirectory_,processName)
             mkdirRobust(outputFolder);
             outputFilePath=fullfile(outputFolder,['detection.mat']); 
             movieInfo=obj.getStruct();
             save(outputFilePath,'movieInfo','-v7.3');
             proc=PointSourceDetectionProcess3D(MD, outputFolder,UTrackPackage3D.getDefaultDetectionParams(MD,outputFolder));
             outFilCell=cell(1,numel(MD.channels_));
             outFilCell{channelIdx}=outputFilePath;
             proc.setOutFilePaths(outFilCell);
             proc.setProcessTag(processName);
         end

         % function proc=saveAsProcess(obj,MD,processTag,channelIdx)
         %    outputFolder=fullfile(fileparts(MD.outFilePaths_,processTag));
         %    mkdirRobust(outputFolder);
         %    outputFilePath=fullfile(outputFolder,['channel_' num2str(channelIdx) '.mat']); 
         %    movieInfo=obj.getStruct();
         %    save(outputFilePath,'movieInfo');
         %    proc=PointSourceDetectionProcess3D(MD, outputFolder,UTrackPackage3D.getDefaultDetectionParams(MD,outputFolder));
         %    proc.setOutFilePaths({outputFilePath,outputFilePath});
         %    proc.setProcessTag(['detectMapped-' name]);
         % end

         function [plotHandle,h]=scatterPlot(obj,varargin)
             ip = inputParser;
             ip.CaseSensitive = false;
             ip.KeepUnmatched=true;
             ip.addOptional('label',{},(@(x) iscell(x)||isnumeric(x)));
             ip.addOptional('handle',[]);
             ip.addOptional('MarkerSize',50);
             ip.addOptional('Opacity',1);
             ip.addOptional('colormap',parula);
             ip.parse(varargin{:});
             p=ip.Results;

             if(isempty(p.handle))
                 figure();
                 h=gca();
             else
                 h=p.handle;
                 hold on;
             end
             st=obj.getAllStruct();

             if(isempty(p.label))
                 label=st.A;
             else
                 label=p.label;
             end

             if(iscell(label))
                 label=vertcat(label{:})';
             end

             plotHandle=scatter3(h,st.x,st.y,st.z,p.MarkerSize,label,'marker','.','MarkerFaceAlpha',p.Opacity);
             % plotHandle.MarkerFaceAlpha=p.Opacity;
             alpha(plotHandle,p.Opacity);
             daspect([1 1 1]);

             colormap(p.colormap);
             colorbar;
         end
         


        function obj=setFromTracks(obj,tracks)
            mi=tracks.getMovieInfo();
            nanIdx=arrayfun(@(mi) isnan(mi.xCoord(:,1)),mi,'unif',0);
            if(isfield(mi,'zCoord'))
                for i=1:length(mi)
                    obj(i).xCoord=mi(i).xCoord(~nanIdx{i},:);
                    obj(i).yCoord=mi(i).yCoord(~nanIdx{i},:);
                    obj(i).zCoord=mi(i).zCoord(~nanIdx{i},:);
                    obj(i).amp=mi(i).amp(~nanIdx{i},:);
                end
            else
                for i=1:length(mi)
                    obj(i).xCoord=mi(i).xCoord(~nanIdx{i},:);
                    obj(i).yCoord=mi(i).yCoord(~nanIdx{i},:);
                    obj(i).zCoord=[];
                    obj(i).amp=mi(i).amp(~nanIdx{i},:);
                end
            end
        end 



        function obj=concatenateOffset(obj,detections)
            if(~iscell(detections))
                detections={detections};
            end

            rOffset=0;
            cOffset=0;
            for rIdx=1:size(detections,1)
                cOffset=0;
                for cIdx=1:size(detections,2)
                    test=obj.getPosMatrix();
                    maxBound=max(vertcat(test{:}));
                    minBound=min(vertcat(test{:}));
                    if(isempty(maxBound))
                        maxBound=0;
                        minBound=0;
                    end
                    cOffset=cOffset+maxBound(1);
                    test=detections{rIdx,cIdx}.getPosMatrix();
                    minBound=min(vertcat(test{:}));
                    if(minBound(1)<0)
                        cOffset=cOffset-minBound(1);
                    end
                    detections{rIdx,cIdx}.addOffset(cOffset,rOffset,0);
                    obj=obj.concatenate(detections{rIdx,cIdx});
                end
                test=obj.getPosMatrix();
                maxBound=max(vertcat(test{:}));
                minBound=min(vertcat(test{:}));
                rOffset=rOffset+maxBound(2);
            end
        end

        function [playerHandle,h]=dynScatterPlot(obj,varargin)
            ip = inputParser;
            ip.CaseSensitive = false;
            ip.KeepUnmatched=true;
            ip.addOptional('colormap',uint8(255*summer(256)));
            ip.addOptional('colorIndex',{});
            ip.addOptional('handle',[]);
            ip.addOptional('timeInterval',0.2);
            ip.addOptional('detections',[]);
            ip.addOptional('MarkerSize',20);
            ip.addOptional('player',[]);
            ip.addOptional('printFilePattern',[]);
            ip.addOptional('view','XY');
            ip.addOptional('show',true);
            ip.parse(varargin{:});
            p=ip.Results;

            posCell=obj.getPosMatrix();
            if(isempty(p.colorIndex))
                minInt=min(obj.getAllStruct().A);
                maxInt=max(obj.getAllStruct().A);
                NC=size(p.colormap,1);
                colorIndex=arrayfun(@(d) (ceil((NC-1)*mat2gray(d.getAllStruct().A,double([minInt maxInt])))+1),obj,'unif',0);
            else
                colorIndex=p.colorIndex;
            end

            colors=cellfun(@(ci) uint8(p.colormap(ci,:)),colorIndex,'unif',0);


            if(isempty(p.player))            
                maxBound=max(vertcat(posCell{:}));
                minBound=min(vertcat(posCell{:}));
                player=pcplayer([minBound(1) maxBound(1)],[minBound(2) maxBound(2)],[minBound(3) maxBound(3)],'MarkerSize',p.MarkerSize);
                player.Axes.Color='k';
                player.Axes.GridColor=[0.9, 0.9, 0.9];
                player.Axes.XColor=[0.9, 0.9, 0.9];
                player.Axes.YColor=[0.9, 0.9, 0.9];
                player.Axes.ZColor=[0.9, 0.9, 0.9];
                player.Axes.Children.SizeData=10;
                player.Axes.Parent.Color='k';
                player.Axes.Parent.InvertHardcopy = 'off';

            else
                player=p.player;
            end

            switch p.view
                case 'XY'
                    az=0;
                    el=90;
                    view(player.Axes,az,el);
                case 'ZY'
                    az=90;
                    el=180;
                    view(player.Axes,az,el);
                case 'ZX'
                    az=180;
                    el=180;
                    view(player.Axes,az,el);
                otherwise
                    body
            end
            if(~isempty(p.printFilePattern))
                mkdirRobust(fileparts(p.printFilePattern));
            end

            if(p.show)
                firstLoop=true;
            while player.isOpen;
                for i=1:numel(posCell)
                    pcUnif=pointCloud(posCell{i},'Color',colors{i});
                    view(player,pcUnif);
                    pause(p.timeInterval);

                    if(~isempty(p.printFilePattern)&&firstLoop)
                        savePath=sprintfPath(p.printFilePattern,i);
                        print(player.Axes.Parent,savePath,'-dpng');
                    end
                end;
                firstLoop=false;
            end
            end
            h=player.Axes;
            playerHandle=player;
        end

        function tracks=buildTracksFromDetection(obj,trackIndices)
            positions=obj.getPosMatrix();

            if(nargin==1)
                trackIndices=arrayfun(@(n) (1:min(size(positions{n},1)))',1:numel(positions),'unif',0);
            end
            maxIndx=max(vertcat(trackIndices{:}));

            % Uses the indices of the first frame as track indices
            M=nan(maxIndx,8*numel(positions));
            for pIdx=1:numel(positions)
                ti=trackIndices{pIdx};
                pos=positions{pIdx};
                if(~isempty(pos))
                    M(ti,(8*(pIdx-1)+1):(8*pIdx))=[pos(:,1) pos(:,2) pos(:,3) pos(:,3) ...
                                                   pos(:,1) pos(:,2) pos(:,3) pos(:,3)];
                end
            end
            tracks=TracksHandle(M);
            for i=1:maxIndx
                tracks(i).tracksFeatIndxCG=i*ones(1,numel(obj));
            end
        end 

        function notNanIdx=getIndexInTracks(obj,tracks)
            mi=tracks.getMovieInfo();
            notNanIdx=arrayfun(@(mi) ~isnan(mi.xCoord(:,1)),mi, ...
                               'unif',0);
            if((numel(obj)-numel(mi))>0)
                paddingIdx=arrayfun(@(c) false(c,1),obj((numel(mi)+1):end).getCard(),'unif',0);
                notNanIdx=[notNanIdx; paddingIdx'];
            end
        end     

        
        function obj= initFromPointCloud(obj,imgLM,vol,ZXRatio)
            lmIdx = find(imgLM~=0);
            [lmy,lmx,lmz] = ind2sub(size(vol), lmIdx);
            N=length(lmy);
            if(nargin>3)
                lmz=lmz*ZXRatio;
            end
            obj=obj.initFromPosMatrices([lmx,lmy,lmz],[lmx,lmy,lmz]);
            obj.amp=[vol(lmIdx) 0.5*ones(N,1)];
        end

        function obj= initFromSegmentedMovie(obj,vol,ZXRatio)
            %% take the centroid of each connected component as a detection.
            featsProp = regionprops(vol,vol,'WeightedCentroid','MaxIntensity');
            nFeats=numel(featsProp);
            weightedCentroid= vertcat(featsProp.WeightedCentroid);
            weightedCentroid(:,3)=ZXRatio*weightedCentroid(:,3);
            obj=obj.initFromPosMatrices(weightedCentroid,0.5*ones(size(weightedCentroid)));
            try
                obj.setAmp(vertcat(featsProp.MaxIntensity));
            catch
                figure()
                imshow(mat2gray(max(vol,[],3)));
                disp('weighted centroid');
                weightedCentroid
                disp('featProp');
                vertcat(featsProp.MaxIntensity)
            end
        end

        function obj=setPosMatrix(obj,posMatrix,errMatrix) 
            if(nargin<3)
                errMatrix=obj.getErrMatrix();
            end          
            if(~iscell(posMatrix))
                posMatrix={posMatrix};
            end
            if(~iscell(errMatrix))
                errMatrix={errMatrix};
            end
            
            nPosMatrix=numel(posMatrix);
            nDet=numel(obj);
            is3D=max(cellfun(@(p) size(p,2),posMatrix))==3;

            if(is3D)
                for i=1:nPosMatrix
                    if(~isempty(posMatrix{min(i,nDet)}))
                        obj(i).xCoord=[posMatrix{min(i,nDet)}(:,1) errMatrix{min(i,nDet)}(:,1)];
                        obj(i).yCoord=[posMatrix{min(i,nDet)}(:,2) errMatrix{min(i,nDet)}(:,2)];
                        obj(i).zCoord=[posMatrix{min(i,nDet)}(:,3) errMatrix{min(i,nDet)}(:,3)];
                        if(all(size(obj(i).amp)~=size(obj(i).xCoord)))
                            obj(i).amp=zeros(size([posMatrix{min(i,nDet)}(:,1) errMatrix{min(i,nDet)}(:,1)]));
                        end
                    end
                end
            else
                for i=1:nPosMatrix
                    if(~isempty(posMatrix{min(i,nDet)}))
                        obj(i).xCoord=[posMatrix{min(i,nDet)}(:,1) errMatrix{min(i,nDet)}(:,1)];
                        obj(i).yCoord=[posMatrix{min(i,nDet)}(:,2) errMatrix{min(i,nDet)}(:,2)];
                        obj(i).zCoord=[];
                        if(all(size(obj(i).amp)~=size(obj(i).xCoord)))
                            obj(i).amp=zeros(size([posMatrix{min(i,nDet)}(:,1) errMatrix{min(i,nDet)}(:,1)]));
                        end
                    end
                end
            end
        end     

        function obj=setAmp(obj,amp,err) 
            if(nargin<3)
                err=obj.getErrMatrix(4);
            end          
            if(~iscell(amp))
                amp={amp};
            end
            if(~iscell(err))
               err={err};
           end
           nPosMatrix=numel(amp);
           nDet=numel(obj);
           for i=1:nDet
            if(~isempty(amp{min(i,nPosMatrix)}))
                obj(i).amp(:,1)=amp{min(i,nPosMatrix)};
                obj(i).amp(:,2)=err{min(i,nPosMatrix)};
            end
            end
        end     

        function [detIndexCell,objIndexCell]=findCloseDetections(obj,det,radii)
            % Find the point in det that are with in a radius from obj.
            % contract: obj and dets have same size
            warning('Warning: inverted obj and det for compliance');

            detIndexCell=cell(1,numel(det));
            objIndexCell=cell(1,numel(det));
            pos=obj.getPosMatrix();
            posDet=det.getPosMatrix();
            parfor i=1:numel(obj)
                disp('findCloseDetections');
                tic;
                if(~isempty(pos{i})&&~isempty(posDet{i}))
                    closeIndx = KDTreeBallQuery(posDet{i},pos{i},radii);
                    objIndex=arrayfun(@(i) i*ones(size(closeIndx{i})),1:numel(closeIndx),'unif',0);
                    objIndexCell{i}=(vertcat(objIndex{:}));
                    detIndexCell{i}=(vertcat(closeIndx{:}));
                else
                    objIndexCell{i}=[];
                    detIndexCell{i}=[];
                end
                toc;
            end
        end

        function [indexCell,distCell]=findClosestDetections(obj,det,radii)
            % For each point in obj, find the closest position in det
            % contract: obj and dets have same size
            indexCell=cell(1,numel(obj));
            distCell=cell(1,numel(obj));
            pos=obj.getPosMatrix();
            pos2=det.getPosMatrix();
            parfor i=1:numel(obj)
                % disp(['findClosestDetections : ' num2str(i)]);
                if(~isempty(pos{i})&&~isempty(pos2{i}))
                    [I,D]=KDTreeBallQuery(pos2{i},pos{i},radii);
                    noMatch=cellfun(@isempty,D);
                    idx=zeros(1,numel(D)); %% idx(10) is the closest idx in det for the 10th in obj(i)
                    idx(~noMatch)=cellfun(@(idxs) idxs(1),I(~noMatch));
                    dista=-1*ones(1,numel(D));
                    dista(~noMatch)=cellfun(@(d) d(1),D(~noMatch));
                    indexCell{i}=idx;
                    distCell{i}=dista;
                end
            end
        end

        function obj=initFromPosMatrices(obj,posMatrix,errMatrix)
            if(~iscell(posMatrix))
                posMatrix={posMatrix};
                errMatrix={errMatrix};
            end
            nPosMatrix=numel(posMatrix);
            obj(nPosMatrix)=Detections();
            obj.setPosMatrix(posMatrix,errMatrix);
        end     

        function detLabel=getLabelFromTracks(obj,tracks,segmentOrTrackLabel,defaultValue)
            detLabel=cell(1,length(obj));
            for i=1:length(obj)
                detLabel{i}=defaultValue*ones(size(obj(i).xCoord,1),1);
            end
            for i=1:length(tracks)
                fr=tracks(i).f;
                featIdx=tracks(i).tracksFeatIndxCG;
                nzIdx = featIdx(1:end-1) ~= 0;
                nzIdx = find(nzIdx(:));
                for pIdx=nzIdx'
                    if(iscell(segmentOrTrackLabel))
                        detLabel{fr(pIdx)}(featIdx(pIdx))=segmentOrTrackLabel{i}(pIdx);
                    else
                        detLabel{fr(pIdx)}(featIdx(pIdx))=segmentOrTrackLabel(i);
                    end
                end
            end
        end 

        function obj=setFromTracksIndx(obj,tracks)
            N=max([tracks.endFrame]);
            pos=cell(1,N)
            for i=1:length(tracks)
                tr=tracks(i);
                fr=tr.f;
                featIdx=tr.tracksFeatIndxCG;
                nzIdx = featIdx ~= 0;
                nzIdx = find(nzIdx(:));
                P=[tr.x' tr.y' tr.z'];
                for pIdx=nzIdx'
                        pos{fr(pIdx)}(featIdx(pIdx),:)=P(pIdx,:);
                end
            end
 
            obj=obj.initFromPosMatrices(pos,pos);
        end
        
        function ret=is3D(obj)
            ret=~isempty(obj.zCoord);
        end     
        
        function pos=getPosMatrix(obj,XYZIndx)
            if(nargin==1)
                XYZIndx=0;
            end
            if((numel(obj)==1))
                if(~isempty(obj.xCoord))
                switch XYZIndx
                case 0
                    pos={[obj.xCoord(:,1),obj.yCoord(:,1),obj.zCoord(:,1)]};
                case 1
                    pos={obj.xCoord(:,1)};
                case 2
                    pos={obj.yCoord(:,1)};
                case 3
                    pos={obj.zCoord(:,1)};
                case 4
                    pos={obj.amp(:,1)};
                otherwise
                end
                else
                    pos={[]};
                end
            else
                pos=arrayfun(@(d) d.getPosMatrix(XYZIndx),obj,'unif',0);
                pos=[pos{:}];
            end
        end

        function err=getErrMatrix(obj,XYZIndx)
            if(nargin==1)
                XYZIndx=0;
            end
            if(numel(obj)==1)
                if(~isempty(obj.xCoord))
                switch XYZIndx
                case 0
                    err={[obj.xCoord(:,2),obj.yCoord(:,2),obj.zCoord(:,2)]};
                case 1
                    err={obj.xCoord(:,2)};
                case 2
                    err={obj.yCoord(:,2)};
                case 3
                    err={obj.zCoord(:,2)};
                case 4
                    err={obj.amp(:,2)};
                otherwise
                end
                else
                    err={[]};
                end
            else
                err=arrayfun(@(d) d.getErrMatrix(XYZIndx),obj,'unif',0);
                err=[err{:}];
            end
        end


        function N=getCard(obj)
                N=arrayfun(@(d) size(d.xCoord,1),obj);
        end

        function obj=setFromSingleTrack(obj,track)
            for fIdx=1:length(obj)
                det=obj(fIdx);
                for i=1:length(track.x)
                    obj(track.f(i)).xCoord=[track.x(i)' 0.5*ones(length(track.lifetime,1))];
                    obj(track.f(i)).yCoord=[track.y(i)' 0.5*ones(length(track.lifetime,1))];
                    obj(track.f(i)).zCoord=[track.z(i)' 0.5*ones(length(track.lifetime,1))];
                end
            end    
        end

        function [obj,lifetimeLabels,trackIndices]=getTracksCoord(obj,tracks)
            pos=nan(3,numel(tracks),max([tracks.endFrame]));
            lifetime=nan(numel(tracks),max([tracks.endFrame]));
            trackIndices=nan(numel(tracks),max([tracks.endFrame]));
            for tIdx=1:numel(tracks)
                pos(1,tIdx,tracks(tIdx).f)=tracks(tIdx).x;
                pos(2,tIdx,tracks(tIdx).f)=tracks(tIdx).y;
                pos(3,tIdx,tracks(tIdx).f)=tracks(tIdx).z;
                lifetime(tIdx,tracks(tIdx).f)=tracks(tIdx).lifetime;
                trackIndices(tIdx,tracks(tIdx).f)=tIdx;
            end
            posCell=arrayfun(@(f) pos(:,~isnan(pos(1,:,f)),f)',1:max([tracks.endFrame]),'unif',0);
            lifetimeLabels=arrayfun(@(f) lifetime(~isnan(pos(1,:,f)),f),1:max([tracks.endFrame]),'unif',0);
            trackIndices=arrayfun(@(f) trackIndices(~isnan(pos(1,:,f)),f),1:max([tracks.endFrame]),'unif',0);
            obj=obj.initFromPosMatrices(posCell,posCell);
        end

        function movieInfo=getStruct(obj)
            movieInfo(length(obj))=struct('xCoord',[],'yCoord',[],'zCoord',[]);
            for fIdx=1:length(obj)
                %progressText(tIdx/length(EB3tracks),'Loading EB3 spherical coordinates.')
                det=obj(fIdx);
                movieInfo(fIdx).xCoord=det.xCoord;
                movieInfo(fIdx).yCoord=det.yCoord;
                movieInfo(fIdx).zCoord=det.zCoord;
                movieInfo(fIdx).amp=det.amp;
                movieInfo(fIdx).int=det.amp;

            end    
        end

        function allStruct=getAllStruct(obj,timeStart,timeInterval)
            if(nargin==1)
                timeStart=0;
                timeInterval=1;
            end
            pos=obj.getPosMatrix();
            structCell=cell(5,numel(obj));
            for fIdx=1:length(pos)
                if(~isempty(pos{fIdx}))
                structCell{1,fIdx}=pos{fIdx}(:,1)';
                structCell{2,fIdx}=pos{fIdx}(:,2)';
                structCell{3,fIdx}=pos{fIdx}(:,3)';
                structCell{4,fIdx}=fIdx*ones(1,numel(structCell{1,fIdx}));
                structCell{5,fIdx}=obj(fIdx).amp(:,1)';
                end
            end
            allStruct.x=[structCell{1,:}];
            allStruct.y=[structCell{2,:}];
            allStruct.z=[structCell{3,:}];
            allStruct.f=[structCell{4,:}];
            allStruct.A=[structCell{5,:}];

            allStruct.t=timeStart+(allStruct.f-1)*timeInterval;

            if(isprop(obj(1),'azimuth'))
               structCell=cell(3,numel(obj));
               for fIdx=1:length(pos)
                   structCell{1,fIdx}=obj(fIdx).azimuth';
                   structCell{2,fIdx}=obj(fIdx).elevation';
                   structCell{3,fIdx}=obj(fIdx).rho';
               end         
               allStruct.azimuth=[structCell{1,:}];
               allStruct.elevation=[structCell{2,:}];
               allStruct.rho=[structCell{3,:}];
           end
        end
          
        function T=getTable(obj)
          T=Detections.detection2table(obj);
        end
        
        function obj=addSphericalCoord(obj)
            for fIdx=1:length(obj)
                %progressText(tIdx/length(EB3tracks),'Loading EB3 spherical coordinates.');
                
                det=obj(fIdx);
                try
                    det.addprop('azimuth');      
                    det.addprop('elevation');    
                    det.addprop('rho');          
                catch
                end
                [det.azimuth,det.elevation,det.rho]=cart2sph(det.xCoord(:,1), ... 
                                                    det.yCoord(:,1), ...
                                                    det.zCoord(:,1));
    
            end    
        end

        function valueCell=interpMD(obj,MD,channelIdx)
            valueCell=cell(1,MD.nFrames_);  
            for fIdx=1:min(numel(obj),MD.nFrames_)
                vol=MD.getChannel(channelIdx).loadStack(fIdx);
                valueCell{fIdx}=obj(fIdx).interpValue(vol,[1 1 MD.pixelSize_/MD.pixelSizeZ_]);
            end
        end


        function vol=getPointCloudVolume(obj,varargin)
            assert(length(obj)==1); 

            ip = inputParser;
            ip.CaseSensitive = false;
            ip.KeepUnmatched=true;
            ip.addOptional('XRange',[],@isnumeric);
            ip.addOptional('YRange',[],@isnumeric);
            ip.addOptional('ZRange',[],@isnumeric);
            ip.parse(varargin{:});
            p=ip.Results;

            XRange=p.XRange;
            YRange=p.YRange;
            ZRange=p.ZRange;

            if(isempty(XRange))
            pos=obj.getPosMatrix();
            maxBound=max(vertcat(pos{:}));
            minBound=min(vertcat(pos{:}));
            XRange=minBound(1):maxBound(1);
            YRange=minBound(2):maxBound(2);
            ZRange=minBound(3):maxBound(3);
            end

            vol=uint16(zeros(numel(YRange),numel(XRange),numel(ZRange)));
            pos=obj.getPosMatrix();
            A=obj.getAllStruct().A;
            pos=pos{1};
            mX=numel(XRange);
            mY=numel(YRange);
            mZ=numel(ZRange);
            for pIdx=1:size(pos,1)
                [~,X]=min(abs( YRange -pos(pIdx,2)));
                [~,Y]=min(abs( XRange -pos(pIdx,1)));
                [~,Z]=min(abs( ZRange -pos(pIdx,3)));
                if((X~=1)&&(Y~=1)&&(Z~=1)&&(X~=mY)&&(Y~=mX)&&(Z~=mZ))
                    vol(X,Y,Z)=A(pIdx);
                end
            end 
        end

        function [MDout,MDFile]=saveTiffVolume(obj,outputDir,XRange,YRange,ZRange)

            mkClrDir(outputDir);
            filePattern=[outputDir filesep 'mask-Vol-frame-%04d.tif'];

            parfor f=1:numel(obj)
                vol=obj(f).getPointCloudVolume(XRange,YRange,ZRange);
                stackWrite(uint16(vol),sprintfPath(filePattern,f));
            end
            channelList=Channel(outputDir);

            tiffReader=TiffSeriesReader({channelList.channelPath_},'force3D',true);
            MDout=MovieData(channelList,[outputDir],'movieDataFileName_','movieData.mat','movieDataPath_',[outputDir], ...
                'pixelSize_',1,'pixelSizeZ_',1,'timeInterval_',1);
            MDout.setReader(tiffReader);
            MDout.sanityCheck(); % the usual fucking performance killer...
            MDout.save();
            MDFile=[outputDir filesep 'movieData.mat'];
        end

        function obj=buildFromMD(obj,MD,channelIdx,ZXRatio)
            %% Each non zero pixel become a detection
            dCell=cell(1,MD.nFrames_);  
            % ZXRatio=MD.pixelSizeZ_/MD.pixelSize_;
            parfor fIdx=1:MD.nFrames_
                d=Detections();
                vol=MD.getChannel(channelIdx).loadStack(fIdx);
                d=d.initFromPointCloud(vol,vol,ZXRatio);
                dCell{fIdx}=d;
            end
            obj=[dCell{:}];
        end

        function obj=buildFromSegmentedMD(obj,MD,channelIdx)
            %% take the centroid of each connected component
            dCell=cell(1,MD.nFrames_);  
            ZXRatio=MD.pixelSizeZ_/MD.pixelSize_;
            for fIdx=1:MD.nFrames_
                d=Detections();
                vol=MD.getChannel(channelIdx).loadStack(fIdx);
                d=d.initFromSegmentedMovie(vol,ZXRatio);
                dCell{fIdx}=d;
            end
            obj=[dCell{:}];
        end

        function values=interpValue(obj,vol,scale)
            if(numel(obj)==1)
                Xq=scale(1)*double(obj.xCoord(:,1))';
                Yq=scale(2)*double(obj.yCoord(:,1))';
                Zq=scale(3)*double(obj.zCoord(:,1))';
                values = interp3(double(vol),Xq,Yq,Zq);
            else
                values=arrayfun(@(d) d.interpValue(vol),obj,'unif',0);
            end
        end
        
        function obj=scale(obj,MD,varargin)
            ip = inputParser;
            ip.CaseSensitive = false;
            ip.KeepUnmatched=true;
            ip.addRequired('MD');
            ip.addOptional('fromTo','isoPixelToReal',@ischar);
            ip.parse(MD,varargin{:});
            p=ip.Results;
            for fIdx=1:length(obj)
                %progressText(tIdx/length(EB3tracks),'Loading EB3 spherical coordinates.');
                det=obj(fIdx);
                switch(p.fromTo)
                case 'isoPixelToReal'
                    det.xCoord(:,1)=(det.xCoord(:,1)-1)*MD.pixelSize_+1;
                    det.yCoord(:,1)=(det.yCoord(:,1)-1)*MD.pixelSize_+1;
                    det.zCoord(:,1)=(det.zCoord(:,1)-1)*MD.pixelSize_+1;
                case 'anisoPixelToReal'
                    det.xCoord(:,1)=(det.xCoord(:,1)-1)*MD.pixelSize_+1;
                    det.yCoord(:,1)=(det.yCoord(:,1)-1)*MD.pixelSize_+1;
                    det.zCoord(:,1)=(det.zCoord(:,1)-1)*MD.pixelSizeZ_+1;
                case 'anisoPixelToPixel'
                    det.zCoord(:,1)=(det.zCoord(:,1)-1)*MD.pixelSizeZ_/MD.pixelSize_+1;
                end
                
            end
        end

        function obj=addOffset(obj,X,Y,Z)
            ip = inputParser;
            ip.CaseSensitive = false;
            ip.KeepUnmatched=true;  
            p=ip.Results;
            if(~iscell(X))
                Z={Z};
                X={X};
                Y={Y};
            end
            for fIdx=1:length(obj)
                det=obj(fIdx);
                if(~isempty(det.xCoord))
                    det.xCoord(:,1)=det.xCoord(:,1)+X{min(end,fIdx)};
                    det.yCoord(:,1)=det.yCoord(:,1)+Y{min(end,fIdx)};
                    det.zCoord(:,1)=det.zCoord(:,1)+Z{min(end,fIdx)};
                end
            end
        end

        function add(obj,detToAdd)
            for fIdx=1:length(obj)
                %progressText(tIdx/length(EB3tracks),'Loading EB3 spherical coordinates.');
                det=obj(fIdx);
                det2=detToAdd(fIdx);
                det.xCoord(:,1)=det.xCoord(:,1)+det2.xCoord(:,1);
                det.yCoord(:,1)=det.yCoord(:,1)+det2.yCoord(:,1);
                det.zCoord(:,1)=det.zCoord(:,1)+det2.zCoord(:,1);
            end
        end
        
        

        
        function obj=concatenate(obj,cdets)
            if(numel(cdets)>numel(obj))
                obj(numel(cdets))=Detections();
            end
            for fIdx=1:numel(cdets)
                det=obj(fIdx);
                det.xCoord=[det.xCoord;cdets(fIdx).xCoord];
                det.yCoord=[det.yCoord;cdets(fIdx).yCoord];
                det.zCoord=[det.zCoord;cdets(fIdx).zCoord];
                det.amp=[det.amp;cdets(fIdx).amp];
                obj(fIdx)=det;
            end
        end
        
        function det=getSelectIdx(obj,indices)
            det=obj.copy();
            det.selectIdx(indices);
        end
        
        function obj=selectIdx(obj,indices)
            if(~iscell(indices))
                for fIdx=1:numel(obj)
                    det=obj(fIdx);
                    det.xCoord=det.xCoord(indices,:);
                    det.yCoord=det.yCoord(indices,:);
                    det.zCoord=det.zCoord(indices,:);
                    det.amp=det.amp(indices,:);
                    if(~isempty(det.labels))
                        fNames=fieldnames(det.labels);
                        for f=1:numel(fNames)
                            l=det.getLabel(fNames{f});
                            l=l{1};
                            det.setLabel(fNames{f},l(indices));
                        end
                    end
                end
            else
                assert(numel(obj)==numel(indices));
                arrayfun(@(i) selectIdx(obj(i),indices{i}),1:numel(indices),'unif',0);
            end
        end        
%         function setFromMovieInfo(obj,aMovieInfo)
%             if isfield(aMovieInfo,'zCoord')
%                 obj.xyz=[aMovieInfo.xCoord(:,1); aMovieInfo.yCoord(:,1); aMovieInfo.zCoord(:,1)]; 
%                 obj.dxyz=[aMovieInfo.xCoord(:,2); aMovieInfo.yCoord(:,2); aMovieInfo.zCoord(:,2)]; 
%             else 
%                 obj.xyz=[aMovieInfo.xCoord(:,1); aMovieInfo.yCoord(:,1)]; 
%                 obj.dxyz=[aMovieInfo.xCoord(:,2); aMovieInfo.yCoord(:,2)]; 
%             endl
%             obj.amp=aMovieInfo.amp(:,1); obj.dAmp=aMovieInfo.amp(:,2);
%         end
%        
        
        function compare(obj)
        end 
        
        function setFromLabel(obj)
        end


    end
    methods(Static)
    function T=detection2table(detections,timeInterval)
        if(nargin<2)
            timeInterval=1;
        end
        tic;
        varName={'t','X','Y','Z','Amp','dX','dY','dZ','dAmp'};
        tableCell=cell(1,length(detections));
        for fIdx=1:length(detections)
            if(~isempty(detections(fIdx).xCoord))
                if(isfield(detections(fIdx),'zCoord'))
                    Z=detections(fIdx).zCoord(:,1);
                    dZ=detections(fIdx).zCoord(:,2);
                else
                    Z=zeros(size((detections(fIdx).xCoord(:,1))));
                    dZ=Z;
                end
                detArray=[timeInterval*ones(size(Z)), detections(fIdx).xCoord(:,1), detections(fIdx).yCoord(:,1),Z,detections(fIdx).amp(:,1), detections(fIdx).xCoord(:,2), detections(fIdx).yCoord(:,2),dZ(:),detections(fIdx).amp(:,2)];
                tableCell{fIdx}=array2table(detArray,'variableNames',varName);
            else
                tableCell{fIdx}=table();
            end
        end
        toc;
        tic;
        T=vertcat(tableCell{:});
        toc;
    end
end
end


%% SNIPPETS

        
%         function setFromPStruct(obj,aPstruct)
%             if isfield(aPstruct,'z')
%                 obj.xyz=[aPstruct.x; aPstruct.y; aPstruct.z]; 
%                 obj.dxyz=[aPstruct.x_pstd ; aPstruct.y_pstd; aPstruct.z_pstd];
%                 obj.scale=[aPstruct.s; aPstruct.s; aPstruct.s];
%             else 
%                 obj.xyz=[aPstruct.x; aPstruct.y;]; 
%                 obj.dxyz=[aPstruct.x_pstd ; aPstruct.y_pstd;];
%                 obj.scale=[aPstruct.s; aPstruct.s;];
%             end
%             obj.amp=aPstruct.A; obj.dAmp=aPstruct.A_pstd;
%             obj.bg=aPstruct.c; obj.dBg=aPstruct.c_pstd;
%         end


%       function P = addprop(obj,propName)
%             if(~isscalar(obj))
%                 P = arrayfun(@(x) addprop(x,propName),obj,'UniformOutput',false);
%                 P = [P{:}];
%                 P = reshape(P,size(obj));
%             else
%                 P = addprop@dynamicprops(obj,propName);
%             end
%         end
%         
%         function out = subsasgn(A,S,B,varargin)
%             try
%                 if(isempty(A))
%                     A = Detections.empty;
%                 end
%                 out = builtin('subsasgn',A,S,B,varargin{:});
%             catch err
%                 switch(err.identifier)
%                     case 'MATLAB:noPublicFieldForClass'
%                         if(~all(isprop(A,S(1).subs)))
%                             rethrow(err);
%                         end
%                         if(nargin < 4)
%                             % Allow for [tracks.prop] = 5; for dynamic
%                             % properties
%                             out = arrayfun(@(t) subsasgn(t,S,B),A,'UniformOutput',false);
%                             out = [out{:}];
%                             out = reshape(out,size(A));
%                         else
%                             % Allow for
%                             % test = {1,2,3,4,5}
%                             % [tracks.prop] = test{:}
%                             % for dynamic properties
%                             out = arrayfun(@(t,b) subsasgn(t,S,b{1}),A,[{B} varargin],'UniformOutput',false);
%                             out = [out{:}];
%                             out = reshape(out,size(A));
%                         end
%                     otherwise
%                         rethrow(err)
%                 end
%             end
%         end
%         function varargout = subsref(A,S)
%             try
%                 [varargout{1:nargout}] = builtin('subsref',A,S);
%             catch err
%                 switch(err.identifier)
%                     case 'MATLAB:noSuchMethodOrField'
%                         if(all(isprop(A,S(1).subs)))
%                             % Allow for tracks.prop where prop is a dynamic
%                             % property
%                             varargout = arrayfun(@(t) subsref(t,S),A,'Unif',false);
%                         else
%                             rethrow(err);
%                         end
%                     otherwise
%                         rethrow(err);
%                 end
%             end
%         end