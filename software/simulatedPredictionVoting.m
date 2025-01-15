function [trackabilityCost,sampledPredictions,sampleLabel,votingLabel,trackabilityCostFull,predExpectation]=simulatedPredictionVoting( ...
    predStat,costFunc,costMatParam,dynCostMatParam,varargin)
    %% Using normally distributed sampling to test the answer of the LAP for any given cost matrix.
%
% Copyright (C) 2025, Danuser Lab - UTSouthwestern 
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
    ip = inputParser;
    ip.CaseSensitive = false;
    ip.KeepUnmatched = true;
    ip.addParameter('sampleNb',20);
    ip.addParameter('debugMode',[false]);
    ip.parse(varargin{:});
    p=ip.Results;

    %calculate cost matrix
    % -- USER DEFINED FUNCTION -- %
    % loading all the dynamic parameter that adapts to the density, lifetime, etc...    
    movieInfo=dynCostMatParam.movieInfo;
    iFrame=dynCostMatParam.iFrame;
    prevCostStruct=dynCostMatParam.prevCostStruct;
    probDim=dynCostMatParam.probDim;
    nnDistFeatures=dynCostMatParam.nnDistFeatures;
    movieInfo=dynCostMatParam.movieInfo;
    featLifetime=dynCostMatParam.featLifetime;
    trackedFeatureIndx=dynCostMatParam.trackedFeatureIndx;
    % eval(['[costMat,propagationScheme,kalmanFilterInfoTmp,nonlinkMarker]'...
    %     ' = ' costFunc '(movieInfo,predStat,'...
    %     'costMatParam,nnDistFeatures,'...
    %     'probDim,prevCostStruct,featLifetime,trackedFeatureIndx,iFrame);'])
    [costMat,propagationScheme,kalmanFilterInfoTmp,nonlinkMarker] = costFunc(...
                movieInfo,predStat,costMatParam,nnDistFeatures,...
                probDim,prevCostStruct,featLifetime,trackedFeatureIndx,iFrame); % Updated by Carmen Klein Herenbrink and Brian Devree

    if any(costMat(:)~=nonlinkMarker) %if there are potential links
         %link features based on cost matrix, allowing for birth and death
        [origLink12,origLink21] = lap(costMat,nonlinkMarker,0);
    end


    nPred=size(predStat.noiseVar,3);
    samples=cell(1,p.sampleNb);
    voteCell=cell(1,p.sampleNb);

    predPos=predStat.stateVec(:,1:probDim);
    noiseVar=predStat.noiseVar(1:probDim,1:probDim,:);
    noiseStdV=sqrt([squeeze(noiseVar(1,1,:)) squeeze(noiseVar(2,2,:)) squeeze(noiseVar(3,3,:))]);

    % PostPrePos is only used for display purposes, hence we will only keep the first propagation scheme.
    postPredPos=kalmanFilterInfoTmp.stateVec(:,1:3,1);
    for sIdx=1:p.sampleNb
        %% We are re sampling data pre-propagation to take advantage of u-track modularity
        %% It is less rigorous but enable trackability estimation on all type of cost matrices. 
        %R = mvnrnd(predPos,sqrt(noiseVar));
        R = predPos + (noiseStdV).*randn(size(predPos));

        samPredStat=predStat;
        samPredStat.stateVec(:,1:probDim)=R;
        % %%
        % eval(['[costMat,propagationScheme,kalmanFilterInfoTmp,nonlinkMarker]'...
        %     ' = ' costFunc '(movieInfo,samPredStat,'...
        %     'costMatParam,nnDistFeatures,'...
        %     'probDim,prevCostStruct,featLifetime,trackedFeatureIndx,iFrame);'])

        [costMat,propagationScheme,kalmanFilterInfoTmp,nonlinkMarker] = costFunc(...
                movieInfo,samPredStat,costMatParam,nnDistFeatures,...
                probDim,prevCostStruct,featLifetime,trackedFeatureIndx,iFrame); % Updated by Carmen Klein Herenbrink and Brian Devree

        % Sample post propagation reflected in cost mat
        samplePostProg=kalmanFilterInfoTmp.stateVec(:,:,1);
        samplePostProg=samplePostProg(:,1:probDim);
        link12=[];
        if any(costMat(:)~=nonlinkMarker) %if there are potential links
             %link features based on cost matrix, allowing for birth and death
            [link12,~] = lap(costMat,nonlinkMarker,0);
        end
        
        %% Voting: 1 : same vote, 2: other particle 3: new Death
        vote=double(link12==origLink12);
        vote(vote==0)=2;
        vote((vote==2)&(link12>nPred))=3;
        vote((nPred+1):end)=[]; % suppress the vote by qdummy variable

        %% storage
        voteCell{sIdx}=vote;
        samples{sIdx}=samplePostProg;
    end
    allSample=vertcat(samples{:});
    sampledPredictions=Detections();
    sampledPredictions.setPosMatrix(allSample,ones(size(allSample)));
    sampleLabel=repmat((1:size(samples{1},1))',length(samples),1);
    votingLabel=vertcat(voteCell{:});
    allVotes=horzcat(voteCell{:});
    trackabilityCost=sum(allVotes==1,2)/size(allVotes,2);
    trackabilityCostFull=[sum(allVotes==1,2) sum(allVotes==2,2) sum(allVotes==3,2)];
    predExpectation=Detections();
    predExpectation.setPosMatrix(postPredPos,ones(size(postPredPos)));

end




