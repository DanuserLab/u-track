function [ dataMatCrpSecMic,projData] = plusTipRemBegEnd(dataBeforeCrop, projData,comp)
%Decided it would be useful to separate this bit of code from the
%reclassification scheme so can call for subregional analysis where this 
% has already been performed- might decide to change later(MB)
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

%% Remove Growths Initiated in First Frame or Ending in Last Frame From Stats
% do this so one does not bias growth lifetime/displacement data (might not
% be what we want for 

if nargin < 3 || isempty(comp) 
comp = 0; 
end 
dataMat = dataBeforeCrop; 

subIdx2rem=[];
% get index of growth and following fgap or bgap (if it exists) that
% begin in the first frame

sF = projData.detectionFrameRange(1,1);
eF = projData.detectionFrameRange(1,2);

%sF=min(dataMat(:,2));

% compound track IDs of all subtracks with nminimum starting frame number

fullIdx2rem=unique(dataMat(dataMat(:,2)==sF,1)); 
for iTr=1:length(fullIdx2rem)
    subIdx=find(dataMat(:,1)==fullIdx2rem(iTr));
    if (length(subIdx)>1) && (dataMat(subIdx(2),5) > 1) % if there is a forward backward gap linked to growth in first frame
        subIdx2rem=[subIdx2rem; subIdx(1:2)]; % Don't remove entire compound track only the fgap or bgap it's linked to
    else
        subIdx2rem=[subIdx2rem; subIdx(1)];
    end
end

 

% get index of growth and preceeding fgap or bgap
% (if it exists) that end in the last frame

%eF=max(dataMat(:,3));
fullIdx2rem=unique(dataMat(dataMat(:,3)==eF,1));
for iTr=1:length(fullIdx2rem)
    subIdx=find(dataMat(:,1)==fullIdx2rem(iTr));
    if (length(subIdx)>1) && (dataMat(subIdx(end-1),5) > 1)
        subIdx2rem=[subIdx2rem; subIdx(end-1:end)]; % take out the last two of list  
    else
        subIdx2rem=[subIdx2rem; subIdx(end)];
    end
end

if comp == 0 
    
% calculate percentages removed  for output (if not a comp/single mat)
lostDataMat = dataMat(subIdx2rem,:); 
numLostGrowth = length(find(lostDataMat(:,5)==1)); 
numLostFgap = length(find(lostDataMat(:,5) == 2)); 
numLostBgap = length(find(lostDataMat(:,5) == 3)); 

totGrowth = length(find(dataMat(:,5)==1)); 
totBgap = length(find(dataMat(:,5)==3)); 
totFgap = length(find(dataMat(:,5) ==2)); 

projData.percentGrowthAtStartOrEnd = 100*(numLostGrowth/totGrowth); 
projData.percentBgapFlankingBegAndEnd = 100*(numLostBgap/totBgap);
projData.percentFgapFlankingBegAndEnd = 100*(numLostFgap/totFgap);

end 

% remove both classes for statistics
dataMat(subIdx2rem,:)=[];

 

dataMatCrpSecMic=dataMat; % NOTE: Kathyrn makes these all absolute values 
% I think for stats it is better to keep sign (MB) 

end

