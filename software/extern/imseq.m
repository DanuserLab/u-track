%IMSEQ Open an image sequence as if it were a video stream.
%
%    h = imseq(filename)
% 
% This class creates an interface to a sequence of images similar to that
% of MATLAB's VideoReader class. However, the read method can only read an
% image at a time.
%
% There are two ways of defining the files to be used:
%  Wildcard - e.g. filename = 'images/*.jpg' will use all files fitting
%             that pattern (as found by DIR), and sort them alphabetically.
%             Note that this doesn't work for image sequences containing
%             (for example) 1.jpg and 10.jpg, as such files will not be
%             sorted numerically.
%  Numbered - e.g. filename = 'images/im.1.jpg' will use all files fitting
%             that numbered pattern, in numerical order, starting at the
%             number given. The format of the filename is assumed to
%             contain an integer, and the  sequence is assumed to increment
%             that integer by one for each consecutive frame, e.g.:
%                 input.98.jpg, input.99.jpg, input.100.jpg, etc.
%             The integer can also be zero padded, e.g.:
%                 0000.png, 0001.png, 0002.png, etc.
%
% IN:
%    filename - string containing the full or partial path to the first
%               frame in an image sequence, or a filename pattern with
%               wildcards.
%
% OUT:
%    h - handle to the stream.
%
%    See also VIDEOREADER, IMSTREAM.
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

% Copyright (C) Oliver Woodford 2011

classdef imseq < hgsetget
    properties (Hidden = true, GetAccess = private, SetAccess = private)
        % Internal properties
        nameFunc;
    end
    properties (SetAccess = private)
        % Visible properties
        NumberOfFrames;
        Duration;
        FrameRate;
        Height;
        Width;
        BitsPerPixel;
        Name;
        VideoFormat;
        Type;
    end
    properties
        % Writable properties
        CurrentTime;
    end
    properties
        % Editable properties
        Tag;
        UserData;
    end
    
    % Functions - copy of VideoReader functions
    methods
        % Contructor
        function this = imseq(sname)
            if any(sname == '*')
                % Directory listing
                [fpath, fname, fext] = fileparts(sname);
                if isempty(fpath)
                    fpath = cd();
                else
                    fpath = cd(cd(fpath));
                end
                fpath = [fpath '/'];
                list = dir([fpath fname fext]);
                list = sort({list(:).name});
                this.NumberOfFrames = numel(list);
                this.nameFunc = @(n) [fpath list{n}];
            else
                % Format string
                this.Name = sname;
                [fpath, fname, fext] = fileparts(which(sname));
                if isempty(fname)
                    [fpath, fname, fext] = fileparts(sname);
                    if isempty(fpath)
                        fpath = cd();
                    else
                        fpath = cd(cd(fpath));
                    end
                end
                % Find the last set of consecutive digits
                [start, finish] = regexp(fname, '[0-9]+');
                if isempty(start)
                    error('No image index found.');
                end
                start = start(end);
                finish = finish(end);
                format_string = sprintf('%s/%s%%.%dd%s%s', fpath, fname(1:start-1), finish-start+1, fname(finish+1:end), fext);
                zero_index = str2double(fname(start:finish)) - 1;
                this.nameFunc = @(n) sprintf(format_string, n+zero_index);
                % Compute sequence length
                fnum = 0;
                while 1
                    fnum = fnum + 1;
                    % Check we can open the file for reading
                    fh = fopen(this.nameFunc(fnum), 'r');
                    if fh == -1
                        break;
                    end
                    fclose(fh);
                end
                this.NumberOfFrames = fnum - 1;
            end
            % Set frame rate and duration
            this.FrameRate = 30;
            this.Duration = this.NumberOfFrames / this.FrameRate;
            % Compute frame info
            A = read(this, 1);
            this.Width = size(A, 2);
            this.Height = size(A, 1);
            switch class(A)
                case {'uint8', 'int8'}
                    this.BitsPerPixel = 8;
                case {'uint16', 'int16'}
                    this.BitsPerPixel = 16;
                case {'uint32', 'int32', 'single'}
                    this.BitsPerPixel = 32;
                case {'uint64', 'int64', 'double'}
                    this.BitsPerPixel = 64;
            end
            this.BitsPerPixel = size(A, 3) * this.BitsPerPixel;
            str = {'Gray%d', '%d', 'RGB%d', 'CMYK%d'};
            this.VideoFormat = sprintf(str{size(A, 3)}, this.BitsPerPixel);
            this.Type = 'imseq';
            this.Tag = '';
            this.CurrentTime = 0;
        end
        % Destructor
        function this = delete(this)
        end
        % Read - same as VideoReader
        function A = read(this, fnum)
            if fnum == Inf
                % Seek to the end
                fnum = this.NumberOfFrames;
            end
            this.CurrentTime = fnum / this.FrameRate;
            if fnum < 1 || fnum > this.NumberOfFrames
                error('Frame %d is not in the range of allowed frames: [1 %d].', fnum, this.NumberOfFrames);
            end
            [A, map] = imread(this.nameFunc(fnum));
            if ~isempty(map)
                A = reshape(map(uint32(A)+1,:), [size(A) size(map, 2)]); % Assume indexed from 0
            end
        end
        % hasFrame - same as VideoReader
        function tf = hasFrame(this)
            tf = this.CurrentTime < this.Duration;
        end
        % readFrame - same as VideoReader
        function A = readFrame(this)
            A = read(this, round(this.CurrentTime * this.FrameRate) + 1);
        end
    end
    % Other VideoReader functions
    methods(Static)
        function b = isPlatformSupported()
            b = true; % Always supported
        end
        function formats = getFileFormats()
            extensions = {'bmp', 'tif', 'tiff', 'jpeg', 'jpg', 'png', 'ppm', 'pgm', 'pbm', 'gif'};
            formats = audiovideo.FileFormatInfo.empty();
            for a = 1:numel(extensions)
                formats(a) = audiovideo.FileFormatInfo(extensions{a}, [upper(extensions{a}) ' file sequence'], true, false);
            end
        end
    end
end