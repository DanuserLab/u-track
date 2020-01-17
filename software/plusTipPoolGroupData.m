function plusTipPoolGroupData(groupData,varargin)
% plusTipPoolGroupData pools plus tip data from multiple projects in groups
%
% SYNOPSIS:  [groupData]=plusTipPoolGroupData(groupList,saveDir,doBtw,doWtn,doPlot,remBegEnd)
%
% INPUT:
% groupList : output of plusTipPickGroups, nProj x 2 cell array where the
%             first column contains the group identifier for each group,
%             and the second column contains the project path
% saveDir   : path to output directory
% doBtw     : 1 to pool data from projects in the same group, 0 to skip
% doWtn     : 1 to save text file of distributions from all projects in
%             each group in a "withinGroupComparisons" folder
% doPlot    : 1 to make histograms and boxplots for within and/or between
%             group data
% remBegEnd : 1 to remove tracks existing at the beginning
%             or end of the movie
%
% OUTPUT:
% groupData : structure containing group information and fields for 9
%             distributions: growth speed (gs), fgap speed (fs), bgap speed
%             (bs), growth lifetime (gl), fgap lifetime (fl), bgap lifetime
%             (bl), growth displacement (gd), fgap displacement (fd), and
%             bgap displacement (bd).
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

% Input check
ip=inputParser;
ip.addRequired('groupData',@(x) isstruct(x) || iscell(x) || isempty(x));
ip.addOptional('saveDir',[],@ischar);
ip.addOptional('doWtn',1,@isscalar);
ip.addOptional('doPlot',1,@isscalar);
ip.parse(groupData,varargin{:});
saveDir=ip.Results.saveDir;
doWtn=ip.Results.doWtn;
doPlot=ip.Results.doPlot;

% Call extract groupData function if empty or string input
if isempty(groupData) || iscell(groupData)
   groupData=plusTipExtractGroupData(groupData); 
end

% Launch interface if empty directory
if isempty(saveDir)
    saveDir=uigetdir(pwd,'Select output directory for pooled group data.');
end

% Set up output directory
if isdir(saveDir), rmdir(saveDir,'s'); end
mkdir(saveDir);

% Within-group comparison
if doWtn
    nGroups = numel(groupData.M);
    for iGroup = 1:nGroups
        % Set up output directory
        wtnDir = [saveDir filesep 'withinGroupComparisons' filesep...
            groupData.names{iGroup}];
        if isdir(wtnDir), rmdir(wtnDir, 's'); end
        mkdir(wtnDir);
        plusTipWithinGroupComparison(groupData, iGroup, wtnDir, doPlot)
    end
end

% Write stats results into a text file
statsFile = [saveDir filesep 'Stats.txt'];
statsNames = fieldnames(groupData.pooledStats{1});
pooledDataStats = cellfun(@struct2cell,groupData.pooledStats,'Unif',false);
pooledDataStats =horzcat(pooledDataStats{:});
fid=fopen(statsFile,'w+');
fprintf(fid,'\t%s',groupData.names{:});
for i=1:numel(statsNames)
    fprintf(fid,'\n%s\t',statsNames{i});
    fprintf(fid,'%g\t',pooledDataStats{i,:});
end
fclose(fid);

if doPlot
    % save histograms of pooled distributions from iGroup
    plusTipMakeHistograms(groupData.M, saveDir, 'labels', groupData.names);
    
    % make between-group boxplots (show pooled data)
    pooledDataMat = cellfun(@(x) vertcat(x{:}),groupData.dataMat,'UniformOutput',false);
    plusTipMakeBoxplots(pooledDataMat',groupData.names',saveDir);
end
