classdef inputParserRetrofit < handle
    %inputParserRetrofit Facilitates the conversion of pre-R2007a functions
    %to use inputParser style arguments by allowing arguments to be specified
    %either as options or as parameters
    %
    % The parameter parsing method is tried first. If that fails due to bad
    % parameter names, then option parsing method is tried.
    %
    % Character inputs should have strict validation functions to avoid
    %
    % Properties
    % EmptyMeansDefault - Interpret empty inputs as desiring default values
    %                     (default: true)
    %
    % Read-Only Properties
    % optionParser - delegate inputParser that will parse args as options
    % paramParser - delegate inputParser that will parse args as parameters
    % useOptions - logical flag to use the optionParser if true
    % 
    %
    % Methods
    % addArgument - Sets up an option in optionParser and a parameter in
    %               paramParser
    %
    % Example
    % 
    % ip = inputParserRetrofit;
    % ip.addArgument('salutation','Hello',@ischar);
    % ip.addArgument('object','World',@ischar);
    % ip.addArgument('punctuation','!',@(x) ismember(x,{'.','!','?'}));
    %
    % saySomething = @(x) disp([ x.salutation ' ' x.object x.punctuation]);
    % 
    % % The following all work
    % ip.parse();
    % saySomething(ip.Results);
    % % Output: Hello World!
    %
    % ip.parse('Goodbye');
    % saySomething(ip.Results);
    % % Output: Goodbye World!
    %
    % ip.parse([],[],'.')
    % saySomething(ip.Results);
    % % Output: Hello World.
    %
    % ip.parse('punctuation','?');
    % saySomething(ip.Results);
    % % Output: Hello World?
    %
    % ip.parse('Bonjour','tout le Monde');
    % saySomething(ip.Results);
    % % Output: Bonjour tout le Monde!
    %
    % ip.parse('object','Labmates','salutation','Goodbye');
    % saySomething(ip.Results);
    % % Output: Goodbye Labmates!
    %
    % Parameter and option styles should not be mixed:
    % ip.parse('Goodbye','object','Labmates');
    % % Error using inputParserRetrofit/parse (line 163)
    % % The value of 'punctuation' is invalid. It must satisfy the function: @(x)isempty(x)||varargin{3}(x).
    %
    % This works because of input validation:
    % ip.parse([],[],'punctuation','?');
    % saySomething(ip.Results);
    % % Output: Hello World?
    %
    % See also inputParser
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
    
    % Mark Kittisopikul, December 2014
    
    properties
        EmptyMeansDefault;
    end
    
    properties(SetAccess = protected)
        optionParser;
        paramParser;
        useOptions;
        numRequired;
    end
    
    methods
        function obj = inputParserRetrofit(varargin)
        % See inputParser constructor
            obj.optionParser = inputParser(varargin{:});
            obj.paramParser = inputParser(varargin{:});
            obj.useOptions = true;
            obj.numRequired = 0;
            obj.EmptyMeansDefault = true;
        end
        function [varargout] = subsref(obj,S)
        % overrides subsref to delegate functions to child inputParsers
            try
                % try to use the native subsref first for this class
                [varargout{1:nargout}] = builtin('subsref',obj,S);
            catch error
                % otherwise proxy to the delegate parsers
                if(strcmp(error.identifier,'MATLAB:noSuchMethodOrField') || ...
                   strcmp(error.identifier,'MATLAB:noPublicFieldForClass'))
               
                    % return the value from the active parser
                    parser = obj.getParser;
                    [varargout{1:nargout}] = subsref(parser,S);
                    
                    % try the other parser, might fail if the other parser
                    % did not parse the input and a Result field is
                    % requested
                    otherParser = obj.getInactiveParser;
                    try
                        subsref(otherParser,S);
                    catch
                        % it's normal for otherParser to fail sometimes
                    end
                else
                    rethrow(error);
                end
            end
        end
        function obj = subsasgn(obj,S,B)
        % overrides subsasgn to delegate assignments to child inputParsers
            try
                % try to assign to this class first
                obj = builtin('subsasgn',obj,S,B);
            catch error
                % make the assignment to both delegate parsers
                if(strcmp(error.identifier,'MATLAB:noSuchMethodOrField') || ...
                   strcmp(error.identifier,'MATLAB:noPublicFieldForClass'))
                    % assignments must be made to both parsers equally
                    obj.optionParser = subsasgn(obj.optionParser,S,B);
                    obj.paramParser = subsasgn(obj.paramParser,S,B);
                else
                    rethrow(error);
                end
            end
        end
        function obj = addArgument(obj,varargin)
        % addArgument adds an option to optionParser and
        %                a parameter to paramParser
        % NB: One can still use addOptional and addParamValue to
        % explicitly add one or the other to both delegate parsers
        %
        % See also inputParser.addParameter and inputParser.addOptional
            obj.paramParser.addParamValue(varargin{:});
            % add empty as a valid input for options
            if(obj.EmptyMeansDefault && nargin > 3)
                varargin{3} = @(x) isempty(x) || varargin{3}(x);
            end
            obj.optionParser.addOptional( varargin{:});
        end
        function obj = addOptional(obj,varargin)
        % see inputParser.addOptional
            % add empty as a valid input for options
            if(obj.EmptyMeansDefault && nargin > 3)
                varargin{3} = @(x) isempty(x) || varargin{3}(x);
            end
            obj.paramParser.addOptional(varargin{:});
            obj.optionParser.addOptional( varargin{:});
        end
        function obj = addRequired(obj,varargin)
        % see inputParser.addRequired
            obj.numRequired = obj.numRequired + 1;
            obj.paramParser.addRequired(varargin{:});
            obj.optionParser.addRequired(varargin{:});
        end
        function parse(obj,varargin)
        % see inputParser.parse
            try
                % try parameter parser first
                obj.paramParser.parse(varargin{:});
                obj.useOptions = false;
            catch error
                % if parameter parse fails because of illegal parameter
                % names then the caller wants to use options
                if(strcmp(error.identifier,'MATLAB:InputParser:ParamMustBeChar') || ...
                   strcmp(error.identifier,'MATLAB:InputParser:UnmatchedParameter') || ...
                   strcmp(error.identifier,'MATLAB:InputParser:ParamMissingValue'))
                    obj.optionParser.parse(varargin{:});
                    obj.useOptions = true;
                else
                    rethrow(error);
                end
            end
            parser = obj.getParser;
            if(obj.EmptyMeansDefault)
                % determine defaults
                p = obj.paramParser.copy;
                p.parse(varargin{1:obj.numRequired});
                defaults = struct2cell(p.Results);

                % set empty values to defaults
                sc = struct2cell(parser.Results);
                useDefaults = cellfun('isempty',sc);
                sc(useDefaults) = defaults(useDefaults);
                sc = cell2struct(sc,fieldnames(parser.Results));
                % Depends on StructExpand, consider alternative
                % Also pass Unmatched through in case KeepUnmatched active
                parser.parse(varargin{1:obj.numRequired},sc,parser.Unmatched);
            end
        end
        function parser = getParser(obj)
        % getParser returns the currently active inputParser
        % only relevant after parse has been called
            if(obj.useOptions)
                parser = obj.optionParser;
            else
                parser = obj.paramParser;
            end
        end
        function inactiveParser = getInactiveParser(obj)
        % getInactiveParser returns the currently inactive parser
        % only relevant after parse has been called
            if(obj.useOptions)
                inactiveParser = obj.paramParser;
            else
                inactiveParser = obj.optionParser;
            end
        end
    end
    
end

