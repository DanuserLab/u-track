function [discrimMat]=plusTipTestDistrib(groupData,saveDir,stringency,testID1,testID2,varargin)
% plusTipTestDistrib returns discrimination matrices for plus tip distributions
%
% SYNOPSIS : [discrimMat]=plusTipTestDistrib(groupData,saveDir,stringency,testID1,testID2)
%
% INPUT : 
% groupData : output of plusTipExtractGroupData
% saveDir   : path to output directory
% stringency: stringency to be applied on the result of the statistical
%             tests
% testID1   : id of the first statistical test to be performed on the data
%             (see discriminationMatrix)
% testID2   : id of the second statistical test to be performed on the data
%             (see discriminationMatrix)
%
% OUTPUT :
%   disrimMat : structure containing matrices with p-values above the
%               diagonal from the permutation test for the means, and with
%               percent confidence below the diagonal from the bootstrapped
%               distribution test, which gives the percent confidence that
%               the two distributions are different.
%
%               the fields with the suffix "cell" are identical to the
%               matrix form except they also include group labels in case
%               you want to write them out into an Excel file.
%  
% Maria Bagonis, April 2011
% Sebastien Besson, Sep 2011
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
ip.addRequired('saveDir',@ischar);
ip.addRequired('stringency',@isscalar);
ip.addRequired('testID1',@isscalar);
ip.addRequired('testID2',@isscalar);
ip.parse(groupData,saveDir,stringency,testID1,testID2);


% Call extract groupData function if empty or string input
if isempty(groupData) || iscell(groupData)
   groupData=plusTipExtractGroupData(groupData); 
end

% Launch interface if empty directory
if isempty(saveDir)
    saveDir=uigetdir(pwd,'Select output directory for pooled data statistics.');
end

% List distribution to be tested
distribNames = {'growth_speed';'fgap_speed';'bgap_speed';...
    'growth_lifetime';'fgap_lifetime';'bgap_lifetime';...
    'growth_displacement';'fgap_displacement';'bgap_displacement'};




% Fill the dataStruct structure array
% remove any fields that cannot be compared because they 
% do not subtracks in one of the groups for that param
toremove ={};
nGroups=length(groupData.names);
dataStruct(nGroups,1)=struct();
count = 1;
for iGroup=1:nGroups
    M =  vertcat(groupData.M{iGroup}{:});
    for iDistrib=1:numel(distribNames)
        field=distribNames{iDistrib};
       
        data = M(~isnan(M(:,iDistrib)),iDistrib); 
        dataStruct(iGroup,1).(field)=data;
   
        
        if (isempty(data) || length(data) == 1)   
           toremove{count,1} = field; 
           count = count +1; 
      
        end 
   
    end
    
    
end


    


% Fill the testStruct structure with corresponding test
for iDistrib = 1:numel(distribNames)
    testStruct.(distribNames{iDistrib}) = [testID1 testID2];
end


if ~isempty(toremove) 
   toremove = unique(toremove); 
   dataStruct=  rmfield(dataStruct,toremove); 
   testStruct = rmfield(testStruct,toremove); 
   
end


% run the test
discrimMat=discriminationMatrix(dataStruct,testStruct);

distribNames = fieldnames(testStruct); 

if ~isempty(toremove)
   discrimMat.NotEnoughTracksForAnalysis = toremove;
end 
% convert the matrices to cell arrays with labels
hits={};
for i=1:numel(distribNames)
    field=distribNames{i};
    
    % add stat name in upper left and labels along rows and columns
    discrimMat.([field '_cell'])=[[field groupData.names];...
        [groupData.names' num2cell(discrimMat.(field))]];

    hitsIdx = find(discrimMat.(field)(:,1) <stringency);

    for iHit = hitsIdx'
        hits{end+1,1} = groupData.names{iHit};
        hits{end,2} = field;
    end
    matLength = length(discrimMat.([field '_cell'])(1,:));
    
    space = cell(1,matLength);
    if i ==1
        discrimMat.All = [discrimMat.([field '_cell']); space];
    else
        discrimMat.All =[discrimMat.All ;discrimMat.([field '_cell'])];
        discrimMat.All = [discrimMat.All; space];
    end
end

% Save results
if isempty(hits), hits = 'No Significant Differences'; end 
discrimMat.hits = hits;
save([saveDir filesep 'discrimMat.mat'],'discrimMat');


end  