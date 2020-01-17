function plusTipPerCellPlot(dataStruct,mean_stdOverAllCellsInGroup,fieldName,groupNames)
%PLUSTIPPERCELLPLOT plots data from the dataStructfor all groups
%  
% Maria Bagonis, 2011
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


nGroups = length(dataStruct);

% properties for plotting
nRep=ceil(nGroups/3);
mrkTpe=repmat(['o';'^';'s'],nRep,1); % we'll make this into a cell later
prop_name(1) = {'Marker'};
prop_name(2) = {'MarkerFaceColor'};
prop_name(3) = {'MarkerEdgeColor'};

cM=jet(nGroups);
cMap=mat2cell(cM,ones(nGroups,1),3);

figure; hold on
for iGroup = 1:nGroups
    
    
    nC=length(dataStruct(iGroup).(fieldName));
    prop_values(1:nC,1) = {mrkTpe(iGroup)};
    prop_values(1:nC,2) = cMap(iGroup,:);
    prop_values(1:nC,3) = cMap(iGroup,:);
    
    % use plot instead of scatter so more flexibility with properties. to
    % do this, make 2 x nPoints matrix where the second row is all NaNs and
    % then use the plot function
    if iGroup == 1
        xStart = 1;
        xEnd = 0;
    else
        xStart = xCoord(1,end) + 1;
        xEnd = xCoord(1,end);
    end
    
    x1 = (xStart:length(dataStruct(iGroup).(fieldName))+xEnd);
    x2 = nan(size(x1));
    
    xCoord = [x1 ; x2];
    xValues(iGroup) = mean(xCoord(1,:));
    y1 = abs(dataStruct(iGroup).(fieldName));
    y2 = nan(size(y1));
    yCoord= [y1 ; y2];
    
    %subplot(2,1,1)
    %hold on
    h =  plot(xCoord,yCoord,'.');
    
    set(h,prop_name,prop_values,'MarkerSize',16)
    % store the handle for the first member of each group for the legend
    h1(iGroup)=h(1);
    
    % these have to be defined each time
    clear prop_values
    
end

hold on

for iGroup = 1:length(groupNames)
    
    yValue = cell2mat(mean_stdOverAllCellsInGroup.(fieldName)(1+iGroup,2));
    
    
    plot(xValues(iGroup), yValue,'Marker','+','MarkerSize',16,'MarkerEdgeColor',cMap{iGroup},'MarkerFaceColor',cMap{iGroup});
    
end

legend(h1,strrep(groupNames,'_',' '),'location','BestOutside');

fieldName = strrep(fieldName, '_',' '); 
ylabel(fieldName); 
xlabel('Cell Count'); 
title(['Scatter Plot of Cellular Distribution for Dynamic Parameter: ' fieldName]); 
end

