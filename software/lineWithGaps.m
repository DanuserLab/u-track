function [ h ] = lineWithGaps( varargin )
%lineWithGaps Plots a single low-level line with gaps between columns in
%the XYZ data
%
% lineWithGaps(X,Y,[Z], parameters ...)
% 
% Parameters:
% -----------
% Handle - Line handle to reuse instead of creating a new line
%
% Delimeter - Delimeter to use between columns. Default: NaN.
%   Possible gap Delimeters are NaN, +Inf, or -Inf.
%
% RemoveOtherGaps - logical value whether to remove other gaps. 
%                   Default: 'false'.
% 
% LineSpec parameters are forwarded to the built-in line function.
%
% Performance: When Delimeter is empty and RemoveOtherGaps is true, this
% function will try to find a delimeter which does not exist in the data.
% This may thus hurt performance.
%
% See also line, matlab.graphics.primitive.Line
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


ip = inputParser;
ip.addRequired('X',@isnumeric);
ip.addRequired('Y',@isnumeric);
ip.addOptional('Z',[],@isnumeric);
ip.addParameter('Handle',-1,@isscalar);
ip.addParameter('Delimeter',NaN,@(x) ~isfinite(x));
ip.addParameter('RemoveOtherGaps',false,@islogical);
ip.KeepUnmatched = true;
ip.parse(varargin{:});

in = ip.Results;
forward = ip.Unmatched;

[forward.X, forward.Y, forward.Z] = joinColumns(in.Delimeter,in.X,in.Y,in.Z);

if(in.RemoveOtherGaps)
    forward.X = removeOtherGaps(forward.X,size(in.X,1));
    forward.Y = removeOtherGaps(forward.Y,size(in.Y,1));
    forward.Z = removeOtherGaps(forward.Z,size(in.Z,1));
end

if(isempty(forward.Z))
    forward = rmfield(forward,'Z');
end

% Reuse handle if valid
if(ishandle(in.Handle))
    h = in.Handle;
    set(h,forward);
else
    h = line(forward);
end

end
function data = removeOtherGaps(data,nRows)
% removeOtherGaps remove non-finite values except for the values every
% nRows+1
    isDelimeter = false(size(data));
    nRows = nRows + 1;
    isDelimeter(nRows:nRows:end) = true;
    data = data(isfinite(data) | isDelimeter);
end