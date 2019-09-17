classdef CachedSequence < handle 
%% Simple cache interface for a XYCT sequences.
%% Caching is very basic. Could be improved with advanced libraries.
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

    properties
        outputPath;
        cache; % A 2D+time 3D stack.
        nFrames;
        nCh;
        imgNDims; %% Each channel can be a RGB images, hence the additional possible dimension;
        cached;
    end
    methods 
    function obj = CachedSequence(outputPath,nFrames,nCh)
        ip = inputParser;
        ip.CaseSensitive = false;
        ip.KeepUnmatched = true;
        ip.addRequired('outputPath',@ischar);
        ip.addRequired('nFrames',@isnumeric);
        ip.addRequired('nCh',@isnumeric);
        ip.parse(outputPath,nFrames,nCh);
        p=ip.Results;
        obj.nFrames=nFrames;
        obj.nCh=nCh;
        obj.outputPath=outputPath;
        obj.cached=false;
    end

    function obj=loadCache(obj)
        if(obj.cached)
            if(isempty(obj.cache))
                tmp=load(obj.outputPath);
                obj.cache=tmp.cache;
            end
        end
    end

    function obj=forceLoadCache(obj)
        if(obj.cached)
            tmp=load(obj.outputPath);
            obj.cache=tmp.cache;
        end
    end
    
    function imdisp(obj)
        figure();
        obj.loadCache();
        imdisp(obj.cache,'size',1);
    end

    function obj=saveCache(obj)
       cache=obj.cache;
       mkdirRobust(fileparts(obj.outputPath));
       save(obj.outputPath,'cache');
       obj.cached=true;
    end

    function obj=swapCache(obj)
       if(~obj.cached)
            obj.saveCache();
       end
        obj.cache=[];
    end

    function obj=emptyCache(obj)  
       obj.cache=[]; 
    end

    function img=loadFrame(obj,fIdx,chIdx)
        if(isempty(obj.cache))
            obj.loadCache();
        end
        if(obj.imgNDims==2)
            img=obj.cache(:,:,fIdx,chIdx);
        else
            img=obj.cache(:,:,:,fIdx,chIdx);
        end
    end

    function img=loadAllChannel(obj,fIdx)
        if(isempty(obj.cache))
            obj.loadCache();
        end
        if(obj.imgNDims==2)
            img=obj.cache(:,:,fIdx,:);
        else
            img=obj.cache(:,:,:,fIdx,:);
        end
    end

    function img=saveFrame(obj,img,fIdx,chIdx)
        if(isempty(obj.cache))
            obj.imgNDims=ndims(img);
            obj.cache=repmat(img,[ones(1,ndims(img)) obj.nFrames obj.nCh]);
        end
        if(ndims(img)==2)
            obj.cache(:,:,fIdx,chIdx)=img;
        else
            obj.cache(:,:,:,fIdx,chIdx)=img;
        end
    end

    function n=getFrameNb(obj)
        n=obj.nFrames;
    end
    end
   
    
    methods(Static)
       function p = testPerf()
            imSize=[400 400];
            nFrames=200;
            nCh=2;

            test{1}=CachedSequence('/tmp/testCachedSequence.mat',nFrames,nCh);
            test{2}=CachedSequenceCellCache('/tmp/CachedSequenceCellCache.mat',nFrames,nCh);
            for tIdx=1:numel(test)
            disp(['test write cache ' num2str(tIdx)]); tic;
            for fIdx=1:nFrames
                for cIdx=1:nCh
                    vol=uint16(rand(imSize));
                    test{tIdx}.saveFrame(vol,fIdx,cIdx);
                end
            end
            disp('done writing cache ');toc;
            disp('Swapping cache');tic;
            test{tIdx}.swapCache();
            disp('Done swapping cache');toc;
            end

            for tIdx=1:numel(test)
            disp(['test read from swap ' num2str(tIdx)]); tic;
            for fIdx=1:nFrames
                for cIdx=1:nCh
                    img=test{tIdx}.loadFrame(fIdx,cIdx);
                end
            end
            disp('done reading');toc;
            end

            for tIdx=1:numel(test)
            disp(['test read from cache ' num2str(tIdx)]); tic;
            for fIdx=1:nFrames
                for cIdx=1:nCh
                    img=test{tIdx}.loadFrame(fIdx,cIdx);
                end
            end
            disp('done reading');toc;
            end
       end
     end
end


% if(nargin>1)
%     obj.buildAndSetOutFilePaths([rawProjectDynROIProcess.getOutputDir() filesep 'Rendering' filesep name],1);
%     set(obj,'ref',rawProjectDynROIProcess.ref);
%     set(obj,'nFrames',length(rawProjectDynROIProcess.nFrames));   
%     [BX,BY,BZ]=rawProjectDynROIProcess.getBoundingBox();
%     obj.setBoundingBox(BX,BY,BZ);
% end