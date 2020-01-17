classdef MaskProcessingProcess < MaskProcess
    %Generic class definition for processes which use / post-process / edit
    %masks which have already been created, but do not themselves directly
    %create masks.
    %
%
% Copyright (C) 2020, Danuser Lab - UTSouthwestern 
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
    
    
    methods(Access = public)
        
        function obj = MaskProcessingProcess(owner,name,funName, funParams,...
                                             inFilePaths,outFilePaths)
           % Constructor of class MaskProcessingProcess
           if nargin == 0
              super_args = {};
           else
               super_args{1} = owner;
               super_args{2} = name;                
           end
           if nargin > 2
               super_args{3} = funName;
           end
           if nargin > 3
               super_args{4} = funParams;
           end
           if nargin > 5
               super_args{5} = outFilePaths;
           end
           
           % Call the superclass constructor - these values are private
           obj = obj@MaskProcess(super_args{:});
           
           if nargin > 4               
              if ~isempty(inFilePaths) && numel(inFilePaths) ...
                      ~= numel(owner.channels_) || ~iscell(inFilePaths)
                 error('lccb:set:fatal','Mask paths must be a cell-array of the same size as the number of image channels!\n\n'); 
              end
              obj.inFilePaths_ = inFilePaths;              
           else
               obj.inFilePaths_ = cell(1,numel(owner.channels_));               
           end
        end
%         function sanityCheck(obj) % throws exception
%             % Sanity Check
%             % 1. Check that non-empty mask paths exist
%             % 2. Check number of mask = number of raw images
%             
%             for i = 1:length(obj.inFilePaths_)
%                 if ~isempty(obj.inFilePaths_{i})
%                     if ~exist(obj.inFilePaths_{i}, 'dir')
%                         error('lccb:set:fatal','Cannot find mask paths %s.\n\n', ...
%                                 obj.inFilePaths_{i});
%                     end                                        
%                     
%                     % Number of mask image files
%                     maskFileNames = imDir(obj.inFilePaths_{i}, true);
%                     
%                     if length(maskFileNames) ~= obj.owner_.nFrames_;
%                         error('lccb:set:fatal', 'Invalid mask input directory: The number of masks in %s is inconsistent with the number of input images in %s.\n\n',...
%                             obj.inFilePaths_{i}, obj.owner_.getChannelPaths(i);
%                     end
%                 end
%             end
%         end        
        function setInMaskPath(obj,iChan,maskPaths)           
            if all(obj.checkChanNum(iChan))
                nChan = numel(iChan);
                if ~iscell(maskPaths)
                    maskPaths = {maskPaths};
                end
                for j = 1:nChan
                    if exist(maskPaths{j},'dir') && ....
                            length(imDir(maskPaths{j})) == obj.owner_.nFrames_;
                        obj.inFilePaths_{iChan(j)} = maskPaths{j};
                    else
                        error(['lccb:set:fatal','Invalid input mask path for channel ' num2str(iChan(j))]);
                    end
                end
            else
                error('lccb:set:fatal','Invalid mask channel number for mask path!\n\n'); 
            end
        end
        function fileNames = getInMaskFileNames(obj,iChan)
            if all(obj.checkChanNum(iChan))
                fileNames = cellfun(@(x)(imDir(x)),obj.inFilePaths_(iChan),'UniformOutput',false);
                fileNames = cellfun(@(x)(arrayfun(@(x)(x.name),x,'UniformOutput',false)),fileNames,'UniformOutput',false);
                nIm = cellfun(@(x)(length(x)),fileNames);
                if ~all(nIm == obj.owner_.nFrames_)                    
                    error('Incorrect number of masks found in one or more channels!')
                end                
            else
                error('Invalid channel numbers! Must be positive integers less than the number of image channels!')
            end    
            
            
        end
       
   
    end

end