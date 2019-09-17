%CACHE Wrapper class for caching slow-to-load data
%
%    buffer = cache(load_func, [cache_len])
%    obj = buffer{key}
%
% This class implements a cache, which can improve efficiency when loading
% slow-to-load objects several times.
%
% IN:
%   load_func - Handle to a function which takes a scalar key as input and
%               returns an object.
%   cache_len - scalar indicating how many objects can be stored in the
%               cache. Default: 1 (cache only the most recent object).
%   key - scalar key which is passed to read_fun to load an object.
%
% OUT:
%   buffer - handle to the cache.
%   obj - cached object.
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

% Copyright (C) Oliver Woodford 2015

classdef cache < handle
    properties (Hidden = true, SetAccess = private)
        load_func; % Function to read in an object
        % Image cache stuff
        buffer;
        cache_indices;
        cache_count;
        load_count;
    end
    
    methods
        % Constructor
        function this = cache(load_fun, buf_size)
            this.load_func = load_fun;
            if nargin < 2
                buf_size = 1; % Default number of images to keep cached
            end
            buf_size = max(buf_size, 1);
            this.buffer = cell(buf_size, 1);
            this.cache_indices = NaN(buf_size, 1);
            this.cache_count = zeros(buf_size, 1);
            this.load_count = 0;
        end
        % The main function - get
        function A = get(this, n)
            if nargin < 2 || ~isscalar(n)
                error('Only one object can be got at a time');
            end
            % Check if buffered
            ind = find(this.cache_indices == n, 1);
            if isempty(ind)
                % Cache the frame
                % Find the least recently used slot
                [ind, ind] = min(this.cache_count);
                % Read in the frame
                this.cache_indices(ind) = n;
                this.buffer{ind} = this.load_func(n);
            end
            % Retrieve the cached frame
            A = this.buffer{ind};
            % Update the count and frame number
            this.load_count = this.load_count + 1;
            this.cache_count(ind) = this.load_count;
        end
        % Forward calls like cache(a) to get
        function A = subsref(this, frame)
            switch frame(1).type
                case {'()', '{}'}
                    if numel(frame(1).subs) ~= 1
                        error('Only one dimensional indexing supported');
                    end
                    A = get(this, frame(1).subs{1});
                case '.'
                    if strcmp(frame(1).subs, 'get')
                        % Forward these references to the relevant method
                        A = builtin('subsref', this, frame);
                    else
                        error('%s is not a public property or method of the cache class.', frame(1).subs);
                    end
            end
        end
    end
end
