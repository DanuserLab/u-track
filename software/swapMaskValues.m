function [M2]=swapMaskValues(M1,oldValues,newValues)
% SWAPMASKVALUES swaps or replaces values of a mask/image
%
% DESCRIPTION: In the case of a matrix with only two values, the function
%              simply swaps the values if oldValues and newValues are not
%              given.  If they are given, the function replaces
%              oldValues(i) with newValues(i).
%
%              e.g.) say you have a binary mask and want to replace the 0's
%              with 2's and the 1's by NaN's.  Then you can run
%              swapMaskValues(M1,[0 1],[2 NaN]).  This function
%              conveniently deals with multiple replacements.
%
% SYNOPSIS: [M2]=swapMaskValues(M1,oldValues,newValues)
%
% INPUT: M1             : starting mask (or any array)
%        oldValues (opt): n-vector containing values to change
%        newValues (opt): n-vector containing corresponding new values
%
% OUTPUT: M2: matrix updated with new values
%
% MATLAB VERSION (originally written on): 7.2.0.232 (R2006a) Windows_NT
%
% USERNAME: kathomps
% DATE: 13-Jun-2007
%
%
%
% Copyright (C) 2018, Danuser Lab - UTSouthwestern 
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




if nargin==2 || nargin>3 % need 1 or 3 arguments
    error('Wrong number of input arguments');
end
if isempty(M1)
    M2=M1;
    return
end

if nargin==1 % assume user wants a simple swap of 2 values
    U=unique(M1);
    UNaNs=sum(isnan(U)); % NaN's get counted separately with unique.m
    if length(U)==1
        if U==1 % assume binary swap is desired
            M2=zeros(size(M1));
            return
        elseif U==0
            M2=ones(size(M1));
            return
        else
            error('Not simple swap: M1 contains one or more than two values')
        end
        
    elseif length(U)==2 || length(U)-UNaNs==1 % either 2 numbers, or 1 number and a nan
        oldValues=[U(1) U(2)];
        newValues=[U(2) U(1)];
    else
        error('Not simple swap: M1 contains one or more than two values')
    end
end
if nargin==3
    if ~isvector(oldValues) || ~isvector(newValues)
        error('oldValues and newValues must both be in vector format and be nonempty');
    end
    if length(oldValues)~=length(newValues)
        error('oldValues and newValues must have same length')
    end
end

M2=double(M1); % initialize output matrix

nSwaps=length(oldValues);

for i=1:nSwaps
    if ~isnan(oldValues(i))
        indexList=find(M1==oldValues(i));
    else
        indexList=find(isnan(M1));
    end
    M2(indexList)=newValues(i);
end