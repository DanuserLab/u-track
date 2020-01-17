classdef ExternalProcess < NonSingularProcess & NameableProcess
    % ExternalProcess A minimalist implementation of a Process that can run
    % an arbitrary function with the ExternalProcess handle as an argument
    %
    % Constructor
    % -----------
    % ExternalProcess(owner,name,fun)
    % owner: A MovieObject such as MovieData
    % name:  (Optional) A char string describing the process
    % fun:   (Optional) A function handle taking the process as an
    %                   argument
    %
    % Static Methods
    % .getName()          : 'Dummy process'
    % .getDefaultParams() : struct()
    %
    % Example
    % -------
    % MD = MovieData;
    % process = ExternalProcess(MD,'Say something',@(p) disp(p.getParameters().text));
    % process.setParameters(struct('text','Hello world!'))
    % process2 = ExternalProcess(MD,'Say something',@(p) disp(p.getParameters().text));
    % process2.setParameters(struct('text','Good bye!'))
    % MD.addProcess(process)
    % MD.addProcess(process2);
    % cellfun(@run,MD.processes_);
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
    properties
        % drawableOutput - struct array
        % .name - name of output used by movieViewer
        % .var  - output requested as a parameter from
        %         loadChannelOutputFcn_
        % .type - one of image, overlay, or graph (see also graphViewer)
        % .formatData - takes the output of loadChannelOutputFcn_ and feeds
        %               it to .defaultDisplayMethod
        % .defaultDisplayMethod - function_handle for constructor of
        %                         MovieDataDisplay. For options, type
        %                         'help movieManagement/Display'
        drawableOutput_ = struct( ...
            'name','something', ...
            'var','output', ...
            'type','image', ...
            'formatData',@(x) x, ...
            'defaultDisplayMethod',@ImageDisplay);
        % loadChannelOutputFcn_ function_handle
        % Takes the process as the first argument and usually the parameter
        % output (see drawableOutput.var)
        loadChannelOutputFcn_ = @(proc,varargin) zeros(proc.getOwner().imSize_);
        % checkChannelOutputFcn_ function_handle
        % Takes the process as the first argument and the channel number
        checkChannelOutputFcn_ = @(proc,iChan) true;
    end
    
    methods(Access = public)
        
        function obj = ExternalProcess(owner, varargin)
            %ExternalProcess(owner,name,fun)
            %owner: A MovieObject such as MovieData
            %name:  (Optional) A char string describing the process
            %fun:   (Optional) A function handle taking the process as an
            %                  argument
            
            % Input check
            ip = inputParser;
            ip.addRequired('owner',@(x) isa(x,'MovieObject'));
            ip.addOptional('name',ExternalProcess.getName(),@ischar);
            ip.addOptional('fun',@(x) x,@(f) validateattributes(f,{'function_handle','char'},{}));
            ip.addParameter('parameters',ExternalProcess.getDefaultParams(owner,varargin{:}), @isstruct);
            ip.addParameter('drawableOutput',[]);
            ip.addParameter('loadChannelOutputFcn',[]);
            ip.addParameter('checkChannelOutputFcn',[]);
            ip.parse(owner,varargin{:});
            
            % Constructor of the ExternalProcess
            super_args{1} = owner;
            super_args{2} = ip.Results.name;
            obj = obj@NonSingularProcess(super_args{:});
            obj.funName_ = ip.Results.fun;

            
            obj.funParams_ = ip.Results.parameters;
            if(ischar(obj.funName_))
                obj.funName_ = str2func(obj.funName_);
            end
            if(~isempty(ip.Results.drawableOutput))
                obj.drawableOutput_ = ip.Results.drawableOutput;
            end
            if(~isempty(ip.Results.loadChannelOutputFcn))
                obj.loadChannelOutputFcn_ = ip.Results.loadChannelOutputFcn;
            end
            if(~isempty(ip.Results.loadChannelOutputFcn))
                obj.checkChannelOutputFcn_ = ip.Results.checkChannelOutputFcn;
            end
                
            
        end
        function output = getDrawableOutput(obj)
            output = obj.drawableOutput_;
        end
        function out = loadChannelOutput(obj,varargin)
            out = obj.loadChannelOutputFcn_(obj,varargin{:});
        end
        function status = checkChannelOutput(obj,iChan)
           nChanTot = numel(obj.owner_.channels_);
           if nargin < 2 || isempty(iChan), iChan = 1:nChanTot; end
           assert(all(obj.checkChanNum(iChan)));
           status = true(size(iChan));
        end


    end


    methods (Static)
        function name = getName()
            name = 'Dummy process';
        end
        function funParams = getDefaultParams(varargin)
            funParams = struct();
        end
        function func = GUI(varargin)
            func = @cliGUI;
        end
        
    end
end
