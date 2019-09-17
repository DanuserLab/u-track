classdef ImAnimation < Animation
%% This concrete class encapsulate animated images saved consecutive files. 
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
        pathTemplate;
        nFrames;
    end
    methods 
    function obj = ImAnimation(pathTemplate,nFrames,varargin)
        ip = inputParser;
        ip.CaseSensitive = false;
        ip.KeepUnmatched = true;
        ip.addRequired('pathTemplate',@ischar);
        ip.addRequired('nFrames',@isnumeric);
        ip.parse(pathTemplate,nFrames,varargin{:});
        p=ip.Results;
        mkdirRobust(fileparts(pathTemplate));
        obj.pathTemplate=pathTemplate;
        obj.nFrames=nFrames;
    end

    function img=loadView(obj,fIdx)
        img=(imread(sprintfPath(obj.pathTemplate,fIdx)));
    end

    function imdisp(obj)
        figure();
        imdisp(obj.loadView(1));
    end
    function obj=saveView(obj,fIdx,img)
        imwrite(img, sprintfPath(obj.pathTemplate, fIdx), 'Compression', 'none');
    end

    function n=getFrameNb(obj)
        n=obj.nFrames;
    end
end
    methods(Static) 
    function obj=buildFromCache(aCachedAnimation,outputPath)
        if(nargin==1)
            obj=ImAnimation([aCachedAnimation.outputPath filesep 'animFrames' filesep 'animFrame_%04d.png'],aCachedAnimation.getFrameNb());
        else
            obj=ImAnimation([outputPath 'frame_%04d.png'],aCachedAnimation.getFrameNb());
        end
        mkdirRobust(fileparts(obj.pathTemplate));
        for fIdx=1:obj.getFrameNb();
            img=aCachedAnimation.loadView(fIdx);
            obj.saveView(fIdx,img);
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


