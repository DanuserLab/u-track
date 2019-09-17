classdef ProjAnimation < Animation
%% This class provide a view from an orthogonal projection rawProjectDynROIProcess.
%% It reuses ProjectDynROIProcess as a backEnd, since this class already the necessary containers. 
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
        hProjectDynROIProcess;
        viewName;
        chIdx;
    end
    methods 
    function obj = ProjAnimation(hProjectDynROIProcess,viewName,varargin)
        ip = inputParser;
        ip.CaseSensitive = false;
        ip.KeepUnmatched = true;
        ip.addRequired('hProjectDynROIProcess',@(x) isa(x,'ProjectDynROIProcess'));
        ip.addRequired('viewName',@(x) (strcmp(x,'XY')||strcmp(x,'ZY')||strcmp(x,'ZX')||strcmp(x,'ortho')));
        ip.addOptional('chIdx',1,@isnumeric);
        ip.parse(hProjectDynROIProcess,viewName,varargin{:});
        p=ip.Results;

        obj.hProjectDynROIProcess = p.hProjectDynROIProcess;
        obj.viewName=p.viewName;
        obj.chIdx=p.chIdx;
    end

    function img=loadView(obj,fIdx)
        [XY,ZY,ZX,ortho]=obj.hProjectDynROIProcess.loadFrame(obj.chIdx,fIdx);
        switch obj.viewName
        case 'XY'
        img=XY;
        case 'ZY'
        img=ZY;
        case 'ZX'
        img=ZX;
        case 'ortho'
        img=ortho;
        otherwise
            error('Unknown view');
        end
    end

    function n=getFrameNb(obj)
        n=obj.hProjectDynROIProcess.frameNb();
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