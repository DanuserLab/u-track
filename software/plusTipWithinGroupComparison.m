function plusTipWithinGroupComparison(groupData, iGroup, varargin)
% plusTipPoolGroupData pools plus tip data from multiple projects in groups
%
% SYNOPSIS:  plusTipWithinGroupComparison(groupData, saveDir, doPlot)
%
% INPUT:
% groupList : output of plusTipExtractGroupData, nProj x 2 cell array
% saveDir   : path to output directory
% doPlot    : 1 to make histograms and boxplots for within and/or between
%             group data
%
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
ip.addRequired('groupData', @(x) isstruct(x) || iscell(x) || isempty(x));
nGroups = numel(groupData.M);
ip.addRequired('iGroup', @(x) isscalar(x) || ismember(x, 1 : nGroups));
ip.addOptional('saveDir', '', @ischar);
ip.addOptional('doPlot', 1, @isscalar);
ip.parse(groupData, iGroup, varargin{:});
saveDir = ip.Results.saveDir;
doPlot = ip.Results.doPlot;

% Call extract groupData function if empty or string input
if isempty(groupData) || iscell(groupData)
   groupData = plusTipExtractGroupData(groupData); 
end

% Launch interface if empty directory
if isempty(saveDir)
    saveDir = uigetdir(pwd,'Select output directory.');
end

% Set up output directory
if isdir(saveDir), rmdir(saveDir, 's'); end
mkdir(saveDir);

% write out speed/lifetime/displacement distributions into a text file
stackedM =  vertcat(groupData.M{iGroup}{:});
dlmwrite([saveDir filesep 'gs_fs_bs_gl_fl_bl_gd_fd_bd_' ...
    groupData.names{iGroup} '.txt'], stackedM, 'precision', 3,...
    'delimiter', '\t','newline', 'pc');

% Names for each movie in iGroup
if ~isempty(groupData.movieNames{iGroup})
    wtnGrpNames = groupData.movieNames{iGroup};
else
    wtnGrpNames = arrayfun(@(x) [groupData.names{iGroup} '_' num2str(x)],...
        1:numel(groupData.dataMat{iGroup}),'Unif',0);
end

% Write stats results into a text file
statsFile = [saveDir filesep 'Stats.txt'];
statsNames = fieldnames(groupData.stats{iGroup}{1});
statsData= cellfun(@struct2cell,groupData.stats{iGroup},'Unif',false);
statsData =horzcat(statsData{:});

pooledStatsNames = fieldnames(groupData.pooledStats{iGroup});
pooledStatsData = struct2cell(groupData.pooledStats{iGroup});
assert(all(ismember(pooledStatsNames, statsNames)));

% Save pooled stats
fid = fopen(statsFile,'w+');
fprintf(fid, '\t%s', wtnGrpNames{:});
fprintf(fid, '\tPooled Data');
for i=1:numel(pooledStatsNames)
    iStat = find(strcmp(pooledStatsNames{i},statsNames),1);
    fprintf(fid,'\n%s\t',statsNames{iStat});
    fprintf(fid,'%g\t',statsData{iStat,:});
    fprintf(fid,'%g',pooledStatsData{i});
end
fclose(fid);

if doPlot==1
    % save histograms of pooled distributions from iGroup
    plusTipMakeHistograms(groupData.M{iGroup}, saveDir);
    
    % make within-group boxplots (show each movie in iGroup)
    plusTipMakeBoxplots(groupData.dataMat{iGroup}', wtnGrpNames', saveDir);
    
    % Plot comets as a function of time (Krek lab request)
    plotDetectedCometsNumber(groupData.detection{iGroup}, saveDir)
end

function plotDetectedCometsNumber(data,saveDir)

maxFrame=max(cellfun(@numel,data));
nProjects = numel(data);
colors=hsv(nProjects);

% define small and large fonts
sfont = {'FontName', 'Helvetica', 'FontSize', 18};
lfont = {'FontName', 'Helvetica', 'FontSize', 22};

% plot
saveFig=figure('PaperPositionMode', 'auto'); % enable resizing
hold on;
arrayfun(@(i) plot(data{i},'-','Color',colors(i,:),'LineWidth',2),1:nProjects);

% Set thickness of axes, ticks and assign tick labels
set(gca, 'LineWidth', 1.5, sfont{:}, 'Layer', 'top','Xlim',[0 maxFrame]);
xlabel('Frame number', lfont{:});
ylabel('Comet number', lfont{:});

% remove box around figure
box off;

saveas(saveFig,[saveDir filesep 'DetectedCometsCount.tif']);
close(saveFig)