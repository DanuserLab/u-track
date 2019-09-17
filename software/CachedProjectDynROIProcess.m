classdef CachedProjectDynROIProcess < ProjectDynROIProcess
%% This class encapsulate a MIP of a view of a 3D movie in a specific frame of reference.
%% This entails the specification of view bounding box and frame of reference.
%% It is designed to be usable for partial rendering "Raw MIP" or "fused channel" in various format,
%% hence the number of rendering channel and frames is taking into account.
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

%% TODO promote method that are compatible with ComputeMIPProcess

    properties (SetAccess = public, GetAccess = public)
       cachedXY;
       cachedZY;
       cachedZX;
       cachedOrtho;
    end

    methods
        function obj = CachedProjectDynROIProcess(owner,varargin)
            ip = inputParser;
            ip.CaseSensitive = false;
            ip.KeepUnmatched = true;
            ip.addRequired('owner',@(x) isa(x,'MovieData'));
            ip.addOptional('name','',@ischar);
            ip.addParameter('nRenderedChannel',length(owner.channels_),@isnumeric);
            ip.addParameter('nFrames',owner.nFrames_,@isnumeric);
            ip.parse(owner,varargin{:});
            p=ip.Results;
            obj = obj@ProjectDynROIProcess(owner,'','nFrames',p.nFrames);
            if(~isempty(p.name))
                buildAndSetOutFilePaths(obj,[owner.outputDirectory_ filesep 'cachedProjections' filesep p.name],p.nRenderedChannel,p.nFrames);
            end
        end

        function outputDir=getOutputDir(obj);
            outputDir=obj.outputDir;
        end

        function obj=imdisp(obj)
            obj.cachedOrtho.imdisp();
        end
        
        function obj=swapCache(obj)
            obj.cachedXY.swapCache();
            obj.cachedZY.swapCache();
            obj.cachedZX.swapCache();
            obj.cachedOrtho.swapCache();
        end

        function obj=emptyCache(obj)
            obj.cachedXY.emptyCache();
            obj.cachedZY.emptyCache();
            obj.cachedZX.emptyCache();
            obj.cachedOrtho.emptyCache();
        end

        function  obj=loadCache(obj)
            obj.cachedXY.loadCache();
            obj.cachedZY.loadCache();
            obj.cachedZX.loadCache();
            obj.cachedOrtho.loadCache();
        end

        function saveFrame(obj,chIdx,fIdx,maxXY,maxZY,maxZX,ortho)

            obj.cachedXY.saveFrame(maxXY,fIdx,chIdx);
            obj.cachedZY.saveFrame(maxZY,fIdx,chIdx);
            obj.cachedZX.saveFrame(maxZX,fIdx,chIdx);

            %% Use Z to index image line (going up)
            if(nargin<7)
                maxZY=permute(maxZY,[2 1 3]);
                maxZX=permute(maxZX,[2 1 3]);
                ortho=projMontage(maxXY,maxZX,maxZY,obj.Zup,obj.ZRight);
            end
            obj.cachedOrtho.saveFrame(ortho,fIdx,chIdx);
        end

        function [maxXY,maxZY,maxZX,three]=loadFrame(obj,chIdx,fIdx)
            maxXY=obj.cachedXY.loadFrame(fIdx,chIdx);
            maxZY=obj.cachedZY.loadFrame(fIdx,chIdx);
            maxZX=obj.cachedZX.loadFrame(fIdx,chIdx);
            three=obj.cachedOrtho.loadFrame(fIdx,chIdx);
        end

        %% path building
        function buildAndSetOutFilePaths(obj,outputDirProj,nChannels,nFrames)
            % Save files in a standardized output for computeMIPProcess 
            % with addtional ROI info and proper file spec
            owner=obj.getOwner();
            obj.cachedXY=CachedSequenceCellCache([outputDirProj filesep 'XY.mat'],nFrames,nChannels);
            obj.cachedZY=CachedSequenceCellCache([outputDirProj filesep 'ZY.mat'],nFrames,nChannels);
            obj.cachedZX=CachedSequenceCellCache([outputDirProj filesep 'ZX.mat'],nFrames,nChannels);
            obj.cachedOrtho=CachedSequenceCellCache([outputDirProj filesep 'Ortho.mat'],nFrames,nChannels);
            obj.outputDir=outputDirProj;
            obj.setOutFilePaths({[outputDirProj filesep 'XY.mat',[outputDirProj filesep 'ZY.mat'],[outputDirProj filesep 'ZX.mat'],[outputDirProj filesep 'Ortho.mat']]});
        end

        function res=frameNb(obj)
            res=obj.nFrames;
        end


    end

    methods (Static)
        function name = getName()
            name = 'DynROI Maximum Intensity Projection';
        end
    end
end
