function [outputFileList]=getFileStackNames(firstfilename)
% getFileStackNames returns a cell array containing all file names (with path) belonging to a stack.
%
% SYNOPSIS [outputFileList]=getFileStackNames(firstFileName)
%
% INPUT    firstFileName: name of the first greyvalue image to be read
%                    including the full path
%                    the actual filename must consist of 
%                    - alphanumeric body
%                    - numeric number
%                    - extension 
%                    - tolarates missing files e.g. 1,2,5,..
%
% OUTPUT   outputFileList: names of all files belonging to the stack
%                          defined by firstFileName. The output is sorted
%                          with respect to the file number.
%
% DEPENDENCES
%
% Aaron Ponti, October 4th, 2002
% modified by Achim Besser, February 2nd, 2010
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

oldDir = [];

% Output
outputFileList = {};

[fpath,fname,fno,fext]=getFilenameBody(firstfilename);

if(isempty(fname) || isempty(fno) || isempty(fext) )
   error('invalid first filename specified');
end;


if(~isempty(fpath))
	% change to stack directory
   oldDir = cd(fpath);
else
   %check if it is in the matlab search path
   tempName=which(firstfilename);
   if ~isempty(tempName)
      [fpath,fname,fno,fext]=getFilenameBody(tempName);
      oldDir = cd(fpath);
   end
end

dirListing = dir;

% get all relevant filenames
iEntry = 1;
fileList = {};
for i = 1:length(dirListing)
   if(~dirListing(i).isdir)
      fileList(iEntry) = {dirListing(i).name};
      iEntry = iEntry + 1;
   end
end

nEntries = 0;

%Identify all relevant files and store them in unsortedOutputFileList:
for fileIndex=1:length(fileList)
    [dummy,currentfname,currentfno,currentfext]=getFilenameBody(fileList{fileIndex});
    %Here strcmpi is case insesitive:
    if strcmpi(currentfname,fname) && str2double(currentfno)>=str2double(fno) && strcmpi(currentfext,fext)
        nEntries=nEntries+1;
        unsortedOutputFileList(nEntries)={strcat(fpath,filesep,fileList{fileIndex})};
        frameNoList(nEntries)=str2double(currentfno);
    end
end

%The outputFileList might be unsorted, this is fixed in the following:
[dummy,newIndx] = sort(frameNoList);
outputFileList=unsortedOutputFileList(newIndx);

% change back to original directory
if(~isempty(oldDir))
   cd(oldDir);
end;