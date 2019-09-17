classdef CachedSequenceCellCache < CachedSequence 
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

    methods 
    function obj = CachedSequenceCellCache(outputPath,nFrames,nCh)
        obj=obj@CachedSequence(outputPath,nFrames,nCh);
    end
    
    function img=loadFrame(obj,fIdx,chIdx)
        if(isempty(obj.cache))
            obj.loadCache();
        end
        eImg=find(~cellfun(@isempty,(obj.cache(:,1))));
        if(fIdx<min(eImg))
            fIdx=min(eImg);
        end
        if(fIdx>max(eImg))
            fIdx=max(eImg);
        end
        img=obj.cache{fIdx,chIdx};
    end
    function imdisp(obj)
        figure();
        obj.loadCache();
        eImg=cellfun(@isempty,(obj.cache(:,1)));
        imdisp(obj.cache(~eImg,1),'size',1);
    end

    function img=loadAllChannel(obj,fIdx)
        if(isempty(obj.cache))
            obj.loadCache();
        end
        img=cat(4,obj.cache{fIdx,:});
    end

    function img=saveFrame(obj,img,fIdx,chIdx)
        if(isempty(obj.cache))
            obj.imgNDims=ndims(img);
            obj.cache=cell(obj.nFrames,obj.nCh);
        end
        obj.cache{fIdx,chIdx}=img;
    end

    function n=getFrameNb(obj)
        n=obj.nFrames;
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