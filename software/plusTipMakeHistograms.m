function plusTipMakeHistograms(rawData,saveDir,varargin)
% plusTipMakeHistograms saves speed, lifetime, and displacement histograms
% for growth, fgap, and bgap populations
%
% Input:
% 
%   rawData - A matrix containing the .
%
%   saveDir - a string containing the 
%
%   plotStd - 0 or 1 if the user wants to plot the standard deviation as an
%   errorbar
%
%   plotSde - 0 or 1 if the user wants to plot the standard error as an
%   errorbar
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

% Called by plusTipPoolGroupData.  See function for input format.
% Kathryn Applegate, Jan 2010
% Sebastien Besson, June 2011

ip = inputParser;
ip.addRequired('rawData',@(x) iscell(x) || isnumeric(x));
ip.addRequired('saveDir',@ischar);
ip.addParamValue('plotStd',0,@isscalar);
ip.addParamValue('plotSte',1,@isscalar);
ip.addParamValue('labels',{},@iscell);
ip.parse(rawData,saveDir,varargin{:});
plotStd=ip.Results.plotStd;
plotSte=ip.Results.plotSte;
labels=ip.Results.labels;

% Convert input to generic cell array of size Nx1 where N is the number of
% groups. Each cell contains M(i)x1 cells where M(i) is the number of
% samples for the i-th group. Each sample cell should be a Lx9 matrix.
if isnumeric(rawData), 
    rawData={{rawData}}; 
elseif isnumeric(rawData{1})
    rawData={rawData};
end
nGroups = numel(rawData);
nCells = cellfun(@numel,rawData);

if ~isdir(saveDir), mkdir(saveDir); end

% Define structures for each quantity/events
eventVals = [{'growth','fgap','bgap'};{[1 0 0],[0 0 1],[0 1 0]}];
eventType = cell2struct(eventVals,{'name','color'});

dataVals = [{'speed','lifetime','displacement'};...
    {'microns/min','sec','microns'};{1:3,4:6,7:9}];
dataType=cell2struct(dataVals,{'name','units','cols'});

plotData(3,3) = struct('rawData',[],'avgBins',[],'stdBins',[],'steBins',[]);

% Create anonymous function to parse raw data for each group/cell
parseSampleRawData = @(groupData,i,j) cellfun(@(x) x(:,dataType(i).cols(j)),...
    groupData,'UniformOutput',false);
parseRawData = @(i,j) cellfun(@(x) parseSampleRawData(x,i,j),rawData,...
    'UniformOutput',false);

% Read raw data for each couple quantity event.
for k=1:9
    [i,j]=ind2sub([3,3],k);
    plotData(i,j).rawData=parseRawData(i,j);
end

binSize = zeros(1,3);
binCenters=cell(1,3);
nBins=25;
% Create 
for i=1:3
    % Read the minimum and maximum value and deduce the bin size
    allCellsData=vertcat(plotData(i,:).rawData);
    allData = vertcat([allCellsData{:}]);
    minData = min(0,nanmin(vertcat(allData{:})));
    maxData = nanmax(vertcat(allData{:}));
    binSize(i) = (maxData-minData)/nBins;
    
    % Create arrays of bin centers
    minEdge = floor(minData/binSize(i))*binSize(i);
    maxEdge = round(maxData/binSize(i))*binSize(i);
    edges = minEdge:binSize(i):maxEdge;
    binCenters{i} = edges(1:end-1)+binSize(i)/2;
end


for k=1:9
    [i,j]=ind2sub([3,3],k);
    % Count data per bin for each cell in each group
    binData = cellfun(@(x) cellfun(@(y) hist(y,binCenters{i})/numel(y),...
        x,'UniformOutput',false),plotData(i,j).rawData,'UniformOutput',false);
    plotData(i,j).binData = cellfun(@(x) vertcat(x{:}),binData,'UniformOutput',false)';
    % Average the results per groups
    plotData(i,j).avgBins = cell2mat(cellfun(@(x) mean(x,1),plotData(i,j).binData,'Unif',false));
    plotData(i,j).stdBins = cell2mat(cellfun(@(x) std(x,0,1),plotData(i,j).binData,'Unif',false));
    plotData(i,j).steBins = cell2mat(cellfun(@(x) std(x,0,1)/sqrt(size(x,1)),plotData(i,j).binData,'Unif',false));
end

for k=1:9
    [i,j]=ind2sub([3,3],k);
    label=[regexprep(dataType(i).name,'(\<[a-z])','${upper($1)}') ' ('  ...
        dataType(i).units ')'];
    
    % Make individual plots for non-empty event data
    %     validBinData = find(~cellfun(@isempty,data(i).avgBins));
    %     for j=validBinData
    saveFig =figure;hold on;
    if nGroups>1, colors=hsv(nGroups); else colors= eventType(j).color; end
    
    h=zeros(1,nGroups);
    for iGroup=1:nGroups
        x = binCenters{i}-binSize(i)/2+(2*iGroup-1)*binSize(i)/(2*nGroups);
        hold on;
        h(iGroup)=bar(x,plotData(i,j).avgBins(iGroup,:),...
            'BarWidth',1/nGroups,'FaceColor',colors(iGroup,:));
        
        % Overlay standard error
        validSte = plotData(i,j).steBins(iGroup,:)~=0;
        if ~isempty(find(validSte,1)) && plotSte
            errorbar(x(validSte),plotData(i,j).avgBins(iGroup,validSte),...
                plotData(i,j).steBins(iGroup,validSte),'.k');
        end
        
        % Overlay standard error
        validStd = plotData(i,j).stdBins(:,iGroup)~=0;
        if ~isempty(find(validStd,1)) && plotStd
            errorbar(x(validStd),plotData(i,j).avgBins(iGroup.validSte),...
                plotData(i,j).stdBins(iGroup,validStd),'.k');
        end
    end

    % Add title, label, limits and save figure
    title([eventType(j).name ' ' dataType(i).name ' distribution'])
    xlabel(label);
    ylabel('Frequency of tracks');
    yLim = get(gca,'YLim');
    set(gca,'YLim',[0 yLim(2)]);
    if nGroups>1, legend(h,strrep(labels,'_',' ')); end
    
    % Save histograms under EPS format
    filename = ['histogram_' dataType(i).name '_' eventType(j).name];
    epsfilename = [filename '.eps'];
    if numel(fullfile(saveDir, filename)) > 128
        % User home directory should be writeable
        homeDir = char(java.lang.System.getProperty('user.home'));
        print(saveFig, '-depsc2', fullfile(homeDir, epsfilename));
        movefile(fullfile(homeDir, epsfilename), saveDir ,'f');
    else
        print(saveFig, '-depsc2', fullfile(saveDir, epsfilename));
    end
     
    close(saveFig)

end

%     figure
%     bar(n,M,'stack')
%     colormap(vertcat(event.color))
%     legend({event.name},'Location','best')
%     title(['Stacked ' data(i).name ' distributions'])
%     xlabel(label);
%     ylabel('Frequency of tracks');
%     saveas(gcf,[saveDir filesep 'histogram_' data(i).name '_stacked.tif'])
%     close(gcf)