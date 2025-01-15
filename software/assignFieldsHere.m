function assignFieldsHere(S,varargin)
%assignFieldsHere assign fields from S to the current workspace
%
% INPUT
% S - struct containing fields to assign in the current workspace
% fieldnames - list of strings naming fields to assign
%
% OUTPUT
% None (Fields are assigned in the current workspace)
%
% EXAMPLE
% S.a = 1;
% S.b = 2;
% assignFieldsHere(S);
% S.a = -1;
% assignFieldsHere(S,'a');
%
% Copyright (C) 2025, Danuser Lab - UTSouthwestern 
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
    if(nargin < 2)
        fields = fieldnames(S);
    else
        fields = varargin;
        assert( all(isfield(S,fields)), ...
            'Field names must be fields of the structure');
    end

    for i=1:length(fields)
        assignin('caller',fields{i},S.(fields{i}));
    end
end

