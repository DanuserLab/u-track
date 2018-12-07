function [groupList, dirName]=combineGroupListFiles(saveResult)
% combineGroupListFiles allows selection of multiple groupList files to combine
%
% SYNOPSIS: [groupList]=combineGroupListFiles(saveResult)
%
% INPUT   : saveResult: 1 to save new groupList file, 0 to return it to
%                       workspace only
% OUTPUT  : groupList  : see plusTipPickGroups.m for details
%
% 2009/08 - Kathryn Applegate, Matlab R2008a
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

if nargin<1
    saveResult=0;
end

groupList=[];
temp=[];
userEntry='Yes';
while strcmp(userEntry,'Yes')
    [fileName,pathName] = uigetfile('*.mat','Select groupList.mat file');
    if fileName==0
        return
    end
    load([pathName filesep fileName]);
    temp=[temp; groupList];
    disp(['Selected: ' pathName fileName])
    userEntry = questdlg('Select another groupList.mat file?');
end
clear groupList;
groupList=temp;

% format path for current OS
nProj=length(groupList);
curDir=pwd;
groupList(:,2)=cellfun(@(x) formatPath(groupList{x,2}),mat2cell([1:nProj]',ones(nProj,1),1),'uniformOutput',0);
cd(curDir)

if saveResult==1
    temp=inputdlg({'Enter file name:'},'',1);
    dirName=uigetdir(pwd,'Select output directory.');
    save([dirName filesep temp{1}],'groupList')
end