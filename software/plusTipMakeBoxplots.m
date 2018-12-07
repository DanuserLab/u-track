function plusTipMakeBoxplots(data,setName,saveDir)
% plusTipMakeBoxplots saves speed, lifetime, and displacement boxplots
% for growth, fgap, and bgap populations from different movie groups
%
% SYNOPSIS: plusTipMakeBoxplots(data,setName,saveDir)
%
% Called by plusTipPoolGroupData.  See function for input format.
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

visible = 'off';

if ~isdir(saveDir)
    mkdir(saveDir)
end

% make within-group boxplots
gIdxInd=cellfun(@(x) find(x(:,5)==1),data,'UniformOutput',0);
fIdxInd=cellfun(@(x) find(x(:,5)==2),data,'UniformOutput',0);
bIdxInd=cellfun(@(x) find(x(:,5)==3),data,'UniformOutput',0);
for iParam=1:3
    switch iParam
        case 1
            if isempty(cell2mat(gIdxInd))
                continue
            end
            ms=cellfun(@(x,y) x(y,4),data,gIdxInd,'UniformOutput',0)';
            ml=cellfun(@(x,y) x(y,6),data,gIdxInd,'UniformOutput',0)';
            md=cellfun(@(x,y) x(y,7),data,gIdxInd,'UniformOutput',0)';
            titleStr='growth';
        case 2
            if isempty(cell2mat(fIdxInd))
                continue
            end
            ms=cellfun(@(x,y) x(y,4),data,fIdxInd,'UniformOutput',0)';
            ml=cellfun(@(x,y) x(y,6),data,fIdxInd,'UniformOutput',0)';
            md=cellfun(@(x,y) x(y,7),data,fIdxInd,'UniformOutput',0)';
            titleStr='fgap';
        case 3
            if isempty(cell2mat(bIdxInd))
                continue
            end
            ms=cellfun(@(x,y) x(y,4),data,bIdxInd,'UniformOutput',0)';
            ml=cellfun(@(x,y) x(y,6),data,bIdxInd,'UniformOutput',0)';
            md=cellfun(@(x,y) x(y,7),data,bIdxInd,'UniformOutput',0)';
            titleStr='bgap';
    end
    maxSize=cellfun(@(x) length(x),ms);
    allMS=nan(max(maxSize,[],2),length(ms));
    for i=1:length(ms)
        allMS(1:maxSize(i),i) = ms{1,i};
    end
    allML=nan(max(maxSize,[],2),length(ml));
    for i=1:length(ml)
        allML(1:maxSize(i),i) = ml{1,i};
    end
    allMD=nan(max(maxSize,[],2),length(md));
    for i=1:length(md)
        allMD(1:maxSize(i),i) = md{1,i};
    end
    % for each value, put group name into matrix
    grpVar=repmat(setName',[max(maxSize'),1]);

    saveFig = figure ('Visible',visible); % boxplot of speeds
    boxplot(allMS(:),grpVar(:),'notch','on','orientation','horizontal');
    title([titleStr ' speeds'])
    set(gca,'YDir','reverse')
    xlabel('speed (microns/min)')
    print(saveFig,'-depsc2',[saveDir filesep 'boxplot_speed_' titleStr '.eps'])
    set(saveFig,'Visible','on'); % So that the .fig is in visible state
    saveas(saveFig,[saveDir filesep 'boxplot_speed_' titleStr '.fig'])
    close(saveFig)

    saveFig = figure ('Visible', visible); % boxplot of lifetimes
    boxplot(allML(:),grpVar(:),'notch','on','orientation','horizontal');
    title([titleStr ' lifetimes'])
    set(gca,'YDir','reverse')
    xlabel('lifetimes (sec)')
    print(saveFig,'-depsc2',[saveDir filesep 'boxplot_lifetime_' titleStr '.eps'])
    set(saveFig,'Visible','on'); % So that the .fig is in visible state
    saveas(saveFig,[saveDir filesep 'boxplot_lifetime_' titleStr '.fig'])
    close(saveFig)

    saveFig = figure ('Visible', visible); % boxplot of displacements
    boxplot(allMD(:),grpVar(:),'notch','on','orientation','horizontal');
    title([titleStr ' displacements'])
    set(gca,'YDir','reverse')
    xlabel('displacement (microns)')
    print(saveFig,'-depsc2',[saveDir filesep 'boxplot_displacement_' titleStr '.eps'])
    set(saveFig,'Visible','on'); % So that the .fig is in visible state
    saveas(saveFig,[saveDir filesep 'boxplot_displacement_' titleStr '.fig'])
    close(saveFig)
end