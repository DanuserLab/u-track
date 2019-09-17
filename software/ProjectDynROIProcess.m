classdef ProjectDynROIProcess <  NonSingularProcess & handle
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
%% WARNING: dependency to <ComputeMIPProcess> was suppressed, ComputeMIPProcess systematic sanity check hit perfs (.

    properties (SetAccess = public, GetAccess = public)
       minXBorder
       maxXBorder 
       minYBorder 
       maxYBorder
       minZBorder
       maxZBorder
       outputDir
       ref
       nFrames
       Zup
       ZRight
    end

    methods
        function obj = ProjectDynROIProcess(owner,varargin)
            ip = inputParser;
            ip.CaseSensitive = false;
            ip.KeepUnmatched = true;
            ip.addRequired('owner',@(x) isa(x,'MovieData'));
            ip.addOptional('name','',@ischar);
            ip.addParameter('nRenderedChannel',length(owner.channels_),@isnumeric);
            ip.addParameter('outputDir',[owner.outputDirectory_ filesep 'projections']);
            ip.addParameter('nFrames',owner.nFrames_,@isnumeric);
            ip.parse(owner,varargin{:});
            p=ip.Results;
            % obj = obj@ComputeMIPProcess(owner,[p.outputDir filesep p.name]);
            obj = obj@NonSingularProcess(owner);
            obj.funName_=(@(p) projectDynROI(owner,'processSingleProj',p));
            obj.maxXBorder=owner.getDimensions('X');
            obj.maxYBorder=owner.getDimensions('Y');
            obj.maxZBorder=owner.getDimensions('Z')*(owner.pixelSizeZ_/owner.pixelSize_);
            obj.minXBorder=1;
            obj.minYBorder=1;
            obj.minZBorder=1;
            obj.nFrames=p.nFrames;
            obj.outputDir=[p.outputDir filesep p.name];
            obj.setProcessTag(p.name);
            obj.Zup=false;
            obj.ZRight=false;

            if(~isempty(p.name))
                buildAndSetOutFilePaths(obj,obj.outputDir,p.nRenderedChannel);
            end
        end

        function outputDir=getOutputDir(obj);
            outputDir=obj.outputDir;
        end


        %% Setter and getter that are not handle by hgsetget inheritance
        function setBoundingBox(obj,minmaxXBorder, minmaxYBorder,minmaxZBorder)
            obj.minXBorder= minmaxXBorder(1);
            obj.maxXBorder= minmaxXBorder(2);
            obj.minYBorder= minmaxYBorder(1);
            obj.maxYBorder= minmaxYBorder(2);
            obj.minZBorder= minmaxZBorder(1);
            obj.maxZBorder= minmaxZBorder(2);
        end
        function [minmaxXBorder, minmaxYBorder,minmaxZBorder]=getBoundingBox(obj)
            minmaxXBorder=[obj.minXBorder obj.maxXBorder];
            minmaxYBorder=[obj.minYBorder obj.maxYBorder];
            minmaxZBorder=[obj.minZBorder obj.maxZBorder];
        end

        function frame=getStartFrame(obj)
            if(~isempty(obj.ref))
                frame=min(obj.ref.frame);
            else
                frame=1;
            end
        end

        function frame=getEndFrame(obj)
            if(~isempty(obj.ref))
                frame=max(obj.ref.frame);
            else
                frame=obj.nFrames;
            end
        end 

        function saveFrame(obj,chIdx,fIdx,maxXY,maxZY,maxZX)
            outFilePaths=obj.outFilePaths_;
            XYFilesPattern = outFilePaths{6, chIdx};
            YZFilesPattern = outFilePaths{7, chIdx};
            XZFilesPattern = outFilePaths{8, chIdx};
            imwrite(maxXY, sprintfPath(XYFilesPattern, fIdx), 'Compression', 'none');
            imwrite(maxZY, sprintfPath(YZFilesPattern, fIdx), 'Compression', 'none');
            imwrite(maxZX, sprintfPath(XZFilesPattern, fIdx), 'Compression', 'none');

            %% Use Z to index image line (going up)
            ThreeFilesPattern = outFilePaths{9, chIdx};
            maxZY=permute(maxZY,[2 1 3]);
            maxZX=permute(maxZX,[2 1 3]);
            three=projMontage(maxXY,maxZX,maxZY,obj.Zup);
            imwrite(three, sprintfPath(ThreeFilesPattern, fIdx), 'Compression', 'none');
        end
        
        function obj=imdisp(obj)
            outFilePaths=obj.outFilePaths_;
            ThreeFilesPattern = outFilePaths{9, 1};
            figure();
            imdisp(imstream(sprintfPath(ThreeFilesPattern, 1)),'Size',1);
        end
        
        function [maxXY,maxZY,maxZX,three]=loadFrame(obj,chIdx,fIdx)
            maxXY=(imread(sprintfPath(obj.outFilePaths_{6,chIdx},fIdx)));
            maxZY=(imread(sprintfPath(obj.outFilePaths_{7,chIdx},fIdx)));
            maxZX=(imread(sprintfPath(obj.outFilePaths_{8,chIdx},fIdx)));
            three=(imread(sprintfPath(obj.outFilePaths_{9,chIdx},fIdx)));
        end
        function obj = importFromDeprecatedExternalProcess(obj,externalProcess)
            % Save files in a standardized output for computeMIPProcess 
            % with addtional ROI info and proper file spec
            outFilePathsOld=externalProcess.outFilePaths_;

            outFilePaths = cell(10, size(outFilePathsOld,2));
            for i =1:size(outFilePaths,2)    
                outFilePaths{1,i} = fileparts(outFilePathsOld{1});
                outFilePaths{2,i} = fileparts(outFilePathsOld{2});
                outFilePaths{3,i} = fileparts(outFilePathsOld{3});
                outFilePaths{4,i} = fileparts(outFilePathsOld{4});
                outFilePaths{5,i} = fileparts(fileparts(outFilePathsOld{4}));
                outFilePaths{6,i} = [outFilePathsOld{1}];
                outFilePaths{7,i} = [outFilePathsOld{2}];
                outFilePaths{8,i} = [outFilePathsOld{3}];
                outFilePaths{9,i} = [outFilePathsOld{4}];
                outFilePaths{10,i}= [fileparts(fileparts(outFilePathsOld{4}));];

                for mIdx=[6:9]
                    mkClrDir(fileparts(outFilePaths{mIdx,i}),false);
                end
                obj.setOutFilePaths(outFilePaths);
            end
        end
        %% path building
        function buildAndSetOutFilePaths(obj,outputDirProj,nChannels)
            % Save files in a standardized output for computeMIPProcess 
            % with addtional ROI info and proper file spec
            owner=obj.getOwner();
            outFilePaths = cell(6, nChannels);
            for i =1:size(outFilePaths,2)    
                outFilePaths{1,i} = [outputDirProj filesep 'ch' num2str(i) filesep 'XY'];
                outFilePaths{2,i} = [outputDirProj filesep 'ch' num2str(i) filesep 'ZY'];
                outFilePaths{3,i} = [outputDirProj filesep 'ch' num2str(i) filesep 'ZX'];
                outFilePaths{4,i} = [outputDirProj filesep 'ch' num2str(i) filesep 'three'];
                outFilePaths{5,i} = [outputDirProj filesep 'ch' num2str(i)];
                outFilePaths{6,i} = [outputDirProj filesep 'ch' num2str(i) filesep 'XY' filesep 'XY_frame_nb%04d.tif'];
                outFilePaths{7,i} = [outputDirProj filesep 'ch' num2str(i) filesep 'ZY' filesep 'ZY_frame_nb%04d.tif'];
                outFilePaths{8,i} = [outputDirProj filesep 'ch' num2str(i) filesep 'ZX' filesep 'ZX_frame_nb%04d.tif'];
                outFilePaths{9,i} = [outputDirProj filesep 'ch' num2str(i) filesep 'three' filesep 'Three_frame_nb%04d.tif'];
                outFilePaths{10,i}= [outputDirProj];
                for mIdx=[6:9]
                    mkClrDir(fileparts(outFilePaths{mIdx,i}),false);
                end
                obj.setOutFilePaths(outFilePaths);

            end
        end

        function res=frameNb(obj)
            res=obj.nFrames;
        end
        function res=getFrameNb(obj)
            res=obj.nFrames;
        end

        function res=getChannelNb(obj)
            res=size(obj.obj.outFilePaths_,2);
        end

        function renderingProcess=renderFused(obj)
            renderingProcess=ProjectDynROIRendering(obj,'fused');
            for fIdx=1:obj.frameNb()
                chCell=cell(4,obj.getChannelNb());
                for cIdx=1:obj.getChannelNb()
                    % [maxXY,maxZY,maxZX,three]=loadFrame(obj,chIdx,fIdx);
                    [chCell{1,cIdx},chCell{2,cIdx},chCell{3,cIdx},chCell{4,cIdx}]=loadFrame(obj,cIdx,fIdx);
                end
                fusedCell=cell(1,obj.getChannelNb());
                for rIdx=1:3
                    fusedCell{rIdx}=sc(chCell(rIdx,:),'stereo');
                end
                renderingProcess.saveFrame(1,fIdx,fusedCell{1},fusedCell{2},fusedCell{3})
            end            
        end


    end

    methods (Static)
        function name = getName()
            name = 'DynROI Maximum Intensity Projection';
        end
        function  getDefaultParams()
            disp('That s silly');
        end
    end
end
