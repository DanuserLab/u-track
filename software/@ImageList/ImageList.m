classdef ImageList < MovieObject
    % Concrete implementation of MovieObject for a list of ImageData objects
%
% Copyright (C) 2026, Danuser Lab - UTSouthwestern 
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
    
    % Qiongjing (Jenny) Zou, May 2026

    properties
        imageListFileName_
        imageListPath_
    end

    properties (SetAccess = protected)
        imageDataFile_
    end

    properties (Transient = true)
        images_
    end

    methods
        function obj = ImageList(images, outputDirectory, varargin)
            if nargin > 0
                if nargin < 2 || isempty(outputDirectory)
                    outputDirectory = pwd;
                end

                if iscellstr(images)
                    obj.imageDataFile_ = images(:)';
                elseif iscell(images) && all(cellfun(@(x) isa(x,'ImageData'), images))
                    obj.imageDataFile_ = cellfun(@getFullPath, images, 'UniformOutput', false);
                    obj.images_ = images(:)';
                elseif isa(images, 'ImageData')
                    obj.imageDataFile_ = arrayfun(@getFullPath, images, 'UniformOutput', false);
                    obj.images_ = num2cell(images);
                else
                    error('lccb:iml:constructor', ...
                        'Images should be a cell array of paths, a cell array of ImageData, or an array of ImageData');
                end

                obj.outputDirectory_ = outputDirectory;
                if ~isempty(obj.outputDirectory_)
                    obj.setPath(obj.outputDirectory_);
                end

                nVarargin = numel(varargin);
                if nVarargin > 1 && mod(nVarargin,2)==0
                    for i=1 : 2 : nVarargin-1
                        obj.(varargin{i}) = varargin{i+1};
                    end
                end
                obj.createTime_ = clock;
            end
        end

        function set.imageListPath_(obj, path)
            endingFilesepToken = [regexptranslate('escape',filesep) '$'];
            path = regexprep(path,endingFilesepToken,'');
            obj.checkPropertyValue('imageListPath_',path);
            obj.imageListPath_ = path;
        end

        function set.imageListFileName_(obj, filename)
            obj.checkPropertyValue('imageListFileName_',filename);
            obj.imageListFileName_ = filename;
        end

        function images = getImages(obj,varargin)
            ip = inputParser;
            allIndex = 1 : obj.getSize();
            ip.addOptional('index',allIndex,@(x) all(ismember(x,allIndex)));
            ip.parse(varargin{:});

            if isempty(obj.images_), obj.sanityCheck; end
            images = obj.images_(ip.Results.index);
        end

        function image = getImage(obj, i)
            assert(isscalar(i) && ismember(i, 1 : obj.getSize()));
            image = obj.images_{i};
        end

        function images = getMovies(obj,varargin)
            % Compatibility alias for MovieList-oriented GUI helpers.
            images = obj.getImages(varargin{:});
        end

        function image = getMovie(obj, i)
            % Compatibility alias for MovieList-oriented GUI helpers.
            image = obj.getImage(i);
        end

        function size = getSize(obj)
            size = numel(obj.imageDataFile_);
        end

        function imageException = sanityCheck(obj, varargin)
            askUser = sanityCheck@MovieObject(obj, varargin{:});

            imageException = cell(1, obj.getSize());
            for i = 1 : obj.getSize()
                try
                    if i <= numel(obj.images_) && ~isempty(obj.getImage(i))
                        [imagePath,imageName,imageExt] = fileparts(obj.imageDataFile_{i});
                        obj.getImage(i).sanityCheck(imagePath,[imageName imageExt], askUser);
                    else
                        fprintf(1,'Loading image %g/%g\n',i,obj.getSize());
                        obj.images_{i} = ImageData.load(obj.imageDataFile_{i}, askUser);
                    end
                catch ME
                    imageException{i} = ME;
                    continue
                end
            end

            if ~all(cellfun(@isempty,imageException))
                ME = MException('lccb:iml:sanitycheck','Failed to load image(s)');
                for i=find(~cellfun(@isempty,imageException))
                    ME = ME.addCause(imageException{i});
                end
                throw(ME);
            end

            disp('Saving image list');
            obj.save();
        end

        function attachImages(obj, images, varargin)
            fullrange = 1:obj.getSize();
            ip = inputParser;
            ip.addRequired('images', @(x) isa(x, 'ImageData'));
            ip.addOptional('range', fullrange, @(x) all(ismember(x, fullrange)));
            ip.parse(images, varargin{:});

            for image = images
                canAttach = strcmp(image.getFullPath(), ...
                    obj.imageDataFile_(ip.Results.range));
                if any(canAttach)
                    j = ip.Results.range(find(canAttach, 1));
                    fprintf(1,'Attaching image %g/%g\n',j,obj.getSize());
                    obj.images_{j} = image;
                end
            end
        end

        function relocate(obj,oldRootDir,newRootDir,full)
            relocate@MovieObject(obj,oldRootDir,newRootDir);
            if nargin<4 || ~full, return; end

            fprintf(1,'Relocating images from %s to %s\n',oldRootDir,newRootDir);
            for i = 1 : obj.getSize()
                obj.imageDataFile_{i} = relocatePath(obj.imageDataFile_{i},oldRootDir,newRootDir);
            end
        end

        function save(ImL,varargin)
            fullPath = ImL.getFullPath();
            assert(~isempty(fullPath), 'Invalid path');

            try
                if exist(fullPath,'file')
                    movefile(fullPath,[fullPath(1:end-3) 'old' '.' datestr(now,'yyyymmddTHHMMSS')],'f');
                end
                save(fullPath, 'ImL');
            catch err
                disp(getReport(err))
            end
        end
    end

    methods(Static)
        function obj = load(varargin)
            if nargin == 0
                [filename,pathname] = uigetfile('*.mat','Select .mat file containing an ImageList object');
                if filename
                    varargin{1} = [pathname filesep filename];
                else
                    obj = ImageList.empty;
                    return;
                end
            end

            assert(strcmpi(varargin{1}(end-3:end), '.mat'), ...
                'Input must be a MAT file');
            [obj, filepath] = ImageList.loadMatFile(varargin{1});
            [imagePath,imageName,imageExt]=fileparts(filepath);
            obj.sanityCheck(imagePath,[imageName imageExt], varargin{2:end});
        end

        function [obj, filepath] = loadMatFile(filepath)
            [obj, filepath] = MovieObject.loadMatFile('ImageList', filepath);
            obj.images_=cell(1,length(obj.imageDataFile_));
            for imIdx=1:length(obj.imageDataFile_)
                if exist(obj.imageDataFile_{imIdx},'file')
                    obj.images_{imIdx}=ImageData.loadMatFile(obj.imageDataFile_{imIdx});
                end
            end
        end

        function status=checkValue(property,value)
            ip = inputParser;
            ip.addRequired('property',@(x) ischar(x) || iscell(x));
            ip.parse(property);
            if iscell(property)
                ip.addRequired('value',@(x) iscell(x)&&isequal(size(x),size(property)));
                ip.parse(property,value);
                status=cellfun(@(x,y) ImageList.checkValue(x,y),property,value);
                return
            end

            validator=ImageList.getPropertyValidator(property);
            propName = regexprep(regexprep(property,'(_\>)',''),'([A-Z])',' ${lower($1)}');
            assert(~isempty(validator),['No validator defined for property ' propName]);
            status = isempty(value) || validator(value);
        end

        function validator = getPropertyValidator(property)
            validator = getPropertyValidator@MovieObject(property);
            if ~isempty(validator), return; end
            if ismember(property, {'imageListPath_','imageListFileName_'})
                validator=@ischar;
            end
        end

        function propName = getPathProperty()
            propName = 'imageListPath_';
        end

        function propName = getFilenameProperty()
            propName = 'imageListFileName_';
        end
    end
end
