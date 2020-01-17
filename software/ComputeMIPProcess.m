classdef  ComputeMIPProcess < ImageProcessingProcess & NonSingularProcess
    % Concrete class for a computing Maximum Intensity Projections (MIP)
    % Andrew R. Jamieson Aug. 2017
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
    
    methods
        function obj = ComputeMIPProcess(owner, varargin)
            
            if nargin == 0
                super_args = {};
            else
                % Input check
                ip = inputParser;
                ip.CaseSensitive = false;
                ip.KeepUnmatched = true;
                ip.addRequired('owner',@(x) isa(x,'MovieData'));
                ip.addOptional('outputDir',owner.outputDirectory_,@ischar);
                ip.addOptional('funParams',[],@isstruct);
                ip.parse(owner,varargin{:});
                outputDir = ip.Results.outputDir;
                funParams = ip.Results.funParams;

                super_args{1} = owner;
                super_args{2} = ComputeMIPProcess.getName;
                super_args{3} = @computeMovieMIP;
                if isempty(funParams)
                    funParams = ComputeMIPProcess.getDefaultParams(owner,outputDir);
                end
                super_args{4} = funParams;
            end
            obj = obj@ImageProcessingProcess(super_args{:});
            obj.is3Dcompatible_ = false; % outputs are 2D
        end

        function h = draw(obj, varargin)
            % Function to draw process output
            outputList = obj.getDrawableOutput();  
                               
            ip = inputParser;
            ip.addRequired('obj',@(x) isa(x,'Process'));
            ip.addRequired('iChan',@isnumeric);
            ip.addOptional('iFrame',[],@isnumeric);
            ip.addOptional('iZ',[], @(x) ismember(x,1:obj.owner_.zSize_));
            ip.addParameter('output', [], @(x) all(ismember(x,{outputList.var})));
            ip.KeepUnmatched = true;
            ip.parse(obj, varargin{:});

            iChan = ip.Results.iChan;
            iOutput = find(cellfun(@(y) isequal(ip.Results.output,y),{outputList.var}));

            % Initialize display method
            if isempty(obj.getDisplayMethod(iOutput,iChan))
                obj.setDisplayMethod(iOutput, iChan,...
                    outputList(iOutput).defaultDisplayMethod(iChan));
            end

            % Display all channels
            switch ip.Results.output
                case {'merged_XY', 'merged_ZY', 'merged_ZX', 'merged_all_three'}

                    switch ip.Results.output
                        case 'merged_XY'
                            iOutput = 1; % see obj.outputfilePaths_
                        case 'merged_ZY'
                            iOutput = 2; % see obj.outputfilePaths_
                        case 'merged_ZX'
                            iOutput = 3; % see obj.outputfilePaths_
                        case 'merged_all_three'
                            iOutput = 4; % see obj.outputfilePaths_
                    end
                        
                    numChan = numel(obj.owner_.channels_);
                    if numChan > 1, cdim=3; else cdim=1; end
                    imData = obj.loadChannelOutput(1, ip.Results.iFrame, 'iOutput', iOutput, 'outputIs3D', false);
                    data = zeros(size(imData,1),size(imData,2), cdim);
                    data(:,:,1) = outputList(iOutput).formatData(imData);
                    
                    if numChan > 1
                        for iChan = 2:numChan
                            imData = obj.loadChannelOutput(iChan, ip.Results.iFrame, 'iOutput', iOutput, 'outputIs3D', false); 
                            data(:,:,iChan) = outputList(iOutput).formatData(imData);
                        end
                    end                  

                    try
                        assert(~isempty(obj.displayMethod_{iOutput,1}));
                    catch ME
                        obj.displayMethod_{iOutput, 1} = ...
                            outputList(iOutput).defaultDisplayMethod();
                    end

                    % Create graphic tag and delegate drawing to the display class
                    tag = ['process' num2str(obj.getIndex()) ip.Results.output 'Output'];
                    h = obj.displayMethod_{iOutput, 1}.draw(data, tag, ip.Unmatched);
                
                case {'XY','ZY','ZX','three'}

                    imData = obj.loadChannelOutput(iChan, ip.Results.iFrame, 'iOutput', iOutput, 'outputIs3D', false);
                    data = outputList(iOutput).formatData(imData);

                    try
                        assert(~isempty(obj.displayMethod_{iOutput,iChan}));
                    catch ME
                        obj.displayMethod_{iOutput, iChan} = ...
                            outputList(iOutput).defaultDisplayMethod();
                    end

                    tag = ['process' num2str(obj.getIndex()) '_channel' num2str(iChan) '_output' num2str(iOutput)];
                    h = obj.displayMethod_{iOutput, iChan}.draw(data, tag, ip.Unmatched);
                
                otherwise
                    error('Incorrect Output Var type');
            end
        end
        
        function output = getDrawableOutput(obj, varargin)
        
            n = 1;
            output(n).name = 'XY';
            output(n).var = 'XY';
            output(n).formatData = @mat2gray;
            output(n).defaultDisplayMethod = @ImageDisplay;
            output(n).type = 'image';
            
            n = length(output)+1;
            output(n).name = 'ZY';
            output(n).var = 'ZY';
            output(n).formatData = @mat2gray;
            output(n).defaultDisplayMethod = @ImageDisplay;
            output(n).type = 'image';
            
            n = length(output)+1;
            output(n).name = 'ZX';
            output(n).var = 'ZX';
            output(n).formatData = @mat2gray;
            output(n).defaultDisplayMethod = @ImageDisplay;
            output(n).type = 'image';

            n = length(output)+1;
            output(n).name = 'three';
            output(n).var = 'three';
            output(n).formatData = @mat2gray;
            output(n).defaultDisplayMethod = @ImageDisplay;
            output(n).type = 'image';

            if numel(obj.owner_.channels_) > 1 
                n = length(output)+1;
                output(n).name = 'merged_all_three';
                output(n).var = 'merged_all_three';
                output(n).formatData = @mat2gray;
                output(n).defaultDisplayMethod = @ImageDisplay;
                output(n).type = 'image';

                n = length(output)+1;
                output(n).name = 'Merged_XY';
                output(n).var = 'merged_XY';
                output(n).formatData = @mat2gray;
                output(n).defaultDisplayMethod = @ImageDisplay;
                output(n).type = 'image';

                n = length(output)+1;
                output(n).name = 'Merged_ZY';
                output(n).var = 'merged_ZY';
                output(n).formatData = @mat2gray;
                output(n).defaultDisplayMethod = @ImageDisplay;
                output(n).type = 'image';

                n = length(output)+1;
                output(n).name = 'Merged_ZX';
                output(n).var = 'merged_ZX';
                output(n).formatData = @mat2gray;
                output(n).defaultDisplayMethod = @ImageDisplay;
                output(n).type = 'image';
            end
        end
    end
    
    methods (Static)
        function name = getName()
            name = 'Maximum Intensity Projection';
        end

        function h = GUI(varargin)
            h = @noSettingsProcessGUI;
        end
        
        function funParams = getDefaultParams(owner, varargin)
            % Input check
            ip=inputParser;
            ip.addRequired('owner',@(x) isa(x,'MovieData'));
            ip.addOptional('outputDir', owner.outputDirectory_, @ischar);
            ip.parse(owner, varargin{:})
            outputDir = ip.Results.outputDir;
            
            % Set default parameters
            funParams.ChannelIndex = 1:numel(owner.channels_);
            funParams.ProcessIndex = []; % empty implies raw images, index used for image input to MIP
            funParams.OutputDirectory = [outputDir  filesep 'MIP'];
        end
    end
end