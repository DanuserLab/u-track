function fImg = locmax2d(img, mask, keepFlat)
%LOCALMAX searches for local maxima in an image
%
%    SYNOPSIS fImg = locmax2d(img, mask, keepFlat)
%
%    INPUT    img    image matrix
%             mask   EITHER a scalar that defines the window dimensions
%                    OR a vector [m n] that defines the window dimensions
%                    OR a binary (0/1) structural element (matrix).
%                    Structural elements such as discs can be defined
%                    using the built-in matlab function "strel".
%                    The input matrix must have an odd number of columns and rows. 
%             keepFlat Optional input variable to choose whether to remove
%                      "flat" maxima or to keep them. Default is 0, to remove them.
%
%    OUTPUT   fImg   image with local maxima (original values) and zeros elsewhere.
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% PARAMETER CHECK
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin<2
   error('Please define all parameters');
end

if numel(mask)==1
    if mod(mask,2)==0
        mask = mask + 1;
    end
    rows = mask;
    cols = mask;
    mask = ones(mask);
    numEl = rows*cols;
    
% If mask is a vector with two entries, then, these two entries define the
% number of rows and columns of a rectangular mask filled with ones:
elseif length(mask(:))==2
    % make sure the mask elements are odd numbers (only then, the 
    % local max operator is properly defined)
    indx = find(~mod(mask,2));
    mask(indx) = mask(indx) + 1;
    rows=mask(1);
    cols=mask(2);
    % number of non-zero elements:
    numEl=prod(mask);
    
    % generate the real mask:
    mask=ones(mask);    
else
    % mask is a flat structural element, with an odd number of cols and
    % rows. Note that this excludes the only possible overlap case [1 1] 
    % which would be treated as identical operation above.
    
    % first check that mask has the right size:
    [rows cols]=size(mask);
    if ~mod(rows,2) || ~mod(cols,2)
        % There is no simple way of extending a general mask, thus:
        error('Mask must have an odd number of rows and columns!');
    end
    
    % get matrix dimensions:
    [rows cols]=size(mask);
    
    % Check that the matrix contains only 0/1. This check could be omitted:
    checkMat = (mask==1 | mask==0);
    if sum(checkMat(:))<rows*cols
        % There is no simple way of extending a general mask, thus:
        error('The mask may contain binary values only!');
    end    
    
    % number of non-zero elements:
    numEl=sum(mask(:));
end

if nargin < 3 || isempty(keepFlat)
    keepFlat = 0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% DEFINITIONS
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% apply a max filter
fImg = ordfilt2(img,numEl,mask);
if keepFlat == 0 %change made by KJ
    fImg2 = ordfilt2(img,numEl-1,mask);
    fImg(fImg2==fImg)=0;
end

% take only those positions where the max filter and the original image value
% are equal -> this is a local maximum
fImg(fImg ~= img) = 0;

% set image border to zero
b = (cols-1)/2;
fImg(:,[1:b end-b+1:end]) = 0;
b = (rows-1)/2;
fImg([1:b end-b+1:end],:) = 0;