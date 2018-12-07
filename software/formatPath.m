function [newPath]=formatPath(oldPath)
% FORMATPATH converts between linux and window paths and vice versa
%
% this function attempts to create a directory name from the input path as
% well as the currently working directory. if the path generated is not an
% existing directory, the function prompts the user to select another
% directory. the user should then select any directory above the
% input directory, or the input directory itself (making sure it does
% in fact exist).  the reason this may occur is the current working
% directory may not be pointing to the same server where the input
% directory exists.
%
% Kathryn Applegate 2008
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


% switch direction of fileseps
if ispc
    temp=strrep(oldPath, '/', '\');       
else
    temp=strrep(oldPath, '\', '/');
end
% if isequal(temp,oldPath)
%     % OS didn't change, nothing to do.
%     newPath=oldPath;
%     return
% end

% check to make sure the input path doesn't contain white space
whiteSpaceIdx=regexp(temp,'\s','start')';
if ~isempty(whiteSpaceIdx)
    error('formatPath: input directory name must not include spaces')
end

doneFlag=0;
% look at current directory
currentDir=[pwd filesep];
% find oldPath's filesep locations
tempFilesepIdx=strfind(temp,filesep);
tryNum=1;
while doneFlag==0
    % find current directory's filesep locations
    currFsepIdx=strfind(currentDir,filesep);

    % concat the first part of the current directory and the second part of
    % the target directory
    k=[]; 
    for i=1:length(currFsepIdx)-1
        finalStr=currentDir(currFsepIdx(end-i)+1:currFsepIdx(end-i+1)-1);
        j=strfind(currentDir,finalStr)-1; % last letter of current dir before match
        k = strfind(temp,finalStr); % first letter of match in temp
        if ~isempty(k)
            break
        end
    end
    
    if isempty(k)
        newPath=[];
    else
        newPath=[currentDir(1:j) temp(k:end)];
    end

    % check if the created path is actually a directory. if not,the root
    % was wrong. ask the user to select a new one.
    if isdir(newPath)
        doneFlag=1;
    else
        if tryNum<=3
            currentDir=uigetdir(pwd,['Select a directory above project directory ' oldPath]);
            cd(currentDir)
            currentDir=[currentDir filesep];
        else
            error('formatPath: data not found. either wrong server or permission denied.')
        end
    end
    tryNum=tryNum+1;
end

