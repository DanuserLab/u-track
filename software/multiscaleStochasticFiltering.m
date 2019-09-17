function [pos,labelSeg,energyMap]=multiscaleStochasticFiltering(vol,XZRatio,varargin)
    ip = inputParser; 
    ip.CaseSensitive = false; 
    ip.KeepUnmatched=true;
    ip.addParameter('scales',[2:0.5:4]);
    ip.addParameter('debug',false);
    ip.addParameter('samplePos',[])
    ip.addParameter('version','');
    ip.addParameter('alpha',0.05);
    ip.addParameter('deepakImplementation',false);
    ip.parse(varargin{:});
    p=ip.Results;

% Created by Philippe Roudot, 2017.  
% Variables cleared out to streamline memory usage by Kevin Dean, 2017.
%
%% Set Scale Range for Detection
% Adhesions setup - [1.1:0.1:4]
%
% Copyright (C) 2019, Danuser Lab - UTSouthwestern 
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
scales=p.scales;

% Myosin setup
% scales=[1.2:0.1:8];

%% Preallocate Cell Arrays and Measure Response in Parallel
disp('computing scales');tic;
masks=cell(1,length(scales));
scaledLoGs=cell(1,length(scales));
imScaledMaps=cell(1,length(scales));
sampleTestMap=~isempty(p.samplePos);
testMaps=cell(1,length(scales));
parfor sIdx=1:length(scales)
    [mask, imgLM, imgLoG,imgLoGScales,testMap]=pointSourceStochasticFiltering(vol,[scales(sIdx) XZRatio*scales(sIdx)],'Alpha',p.alpha);
    masks{sIdx}=mask;

    if(p.deepakImplementation)
        flagNormalizeGaussian=true;
        flagUseGPU=false;
        [ imgLoGScales ] = filterLoGND( vol, scales(sIdx), ... 
                                        'spacing', [1 1 1/XZRatio], ...
                                        'borderCondition', 'symmetric', ...
                                        'UseNormalizedDerivatives', true, ...
                                        'UseNormalizedGaussian', flagNormalizeGaussian, ...
                                        'UseGPU', flagUseGPU );
        scaledLoGs{sIdx}=-imgLoGScales;
    else
        scaledLoGs{sIdx}=imgLoGScales;
    end
    if(p.debug)
        imScaledMaps{sIdx}=[{imgLoGScales},testMap];
    end
    if(sampleTestMap)
    testMaps{sIdx}=testMap;
    end
end
toc;


scaleSpace=cat(4,scaledLoGs{:});
[maxResponseMap,maxResponseScale]=nanmax(scaleSpace,[],4);
[sumResponseMap]=nansum(scaleSpace,4);

voteMap=single(cat(4,(masks{:})));
voteMap=nansum(voteMap,4);
% voteMap(voteMap==1)=0;
% maxResponseScale(maxResponseScale==length(scales))=0;
scaleVol=maxResponseScale;

switch p.version
    case 'useVoteMap'
        energyMap=voteMap;
        % energyMap(energyMap==1)=0;

    case 'useMaxScale'
        energyMap=-maxResponseScale;
        energyMap(voteMap==0)=0;
        % energyMap(maxResponseScale==1)=0;
        % energyMap(energyMap==length(scales))=0;
    case 'useMaxResponse'
        energyMap=maxResponseMap;
        energyMap(voteMap==0)=0;
        % energyMap(maxResponseScale==numel(scales))=0;

    case 'useMaxResponseAndVote'
        energyMap=maxResponseMap.*voteMap;
        energyMap(voteMap==1)=0;
    otherwise
        disp('Using default energy map: maximum Laplacian of Gaussian response');
        energyMap=maxResponseMap;
end

bassinVol=-double(energyMap);
bassinVol(energyMap==0)=inf;
labelSeg = watershed(bassinVol);
labelSeg(energyMap==0) = 0;
energyMapPlain=energyMap;
energyMap(labelSeg==0)=0;



%% In situ debugging through location sampling
if(~isempty(p.samplePos))
volSample=p.samplePos.interpValue(vol,[1 1 XZRatio]);
p.samplePos.setLabel('vol',volSample);

samples=p.samplePos.interpValue(maxResponseScale,[1 1 XZRatio]);
p.samplePos.setLabel('maxResponseScale',samples);

samples=p.samplePos.interpValue(voteMap,[1 1 XZRatio]);
p.samplePos.setLabel('voteMap',samples);

samples=p.samplePos.interpValue(energyMap,[1 1 XZRatio]);
p.samplePos.setLabel('energyMap',samples);

val=p.samplePos.interpValue(maxResponseMap,[1 1 XZRatio]);
p.samplePos.setLabel('maxResponseMap',val);

% val=cell(1,numel(scaledLoGs));
% for sIdx=1:numel(scaledLoGs)
%     val{sIdx}=p.samplePos.interpValue(scaledLoGs{sIdx},[1 1 XZRatio]);
% end
% val=[val{:}];
% p.samplePos.setLabel('scaledLoGs',val);

% val=cell(1,numel(testMaps));
% for sIdx=1:numel(testMaps)
%     val{sIdx}=p.samplePos.interpValue(testMaps{sIdx},[1 1 XZRatio]);
% end
% val=[val{:}];
% p.samplePos.setLabel('testMaps',val);
end


%% Volume rendering.
pos=labelToMovieInfo(labelSeg,vol);
if(p.debug)
    disp('debug');
%     figure();
%     imseriesmaskshow(vol, masks,'maskColors',hsv(numel(masks)));
    unique(maxResponseScale(:))
    [rgbScale] = label2rgbND(maxResponseScale+1, [0 0 0; parula(max(maxResponseScale(:)));[1 0 0]]);

    bins=linspace(min(maxResponseMap(:)),max(maxResponseMap(:)),1000);
    maxResponseMapD=discretize(maxResponseMap,bins);   
    [rgbResponseMap] = label2rgbND(maxResponseMapD+1, [0 0 0; parula(numel(bins)+1)] );
    
    bins=linspace(min(energyMapPlain(:)),max(energyMapPlain(:)),1000);
    energyMapD=discretize(energyMapPlain,bins);   
    [rgbEnergy] = label2rgbND(energyMapD+1, [0 0 0; parula(numel(bins)+1)] );

    [rgbLabel] = label2rgbND(labelSeg+1, [0 0 0; (prism(numel(unique(labelSeg))+1))] );
    imseriesmaskshowrgb(vol,{rgbScale,rgbResponseMap,rgbEnergy,rgbLabel},'displayRange',[0 max(vol(:))]);
    if(iscell(imScaledMaps{1}))
        for i=1:numel(imScaledMaps{1})
            allScaleMap=cellfun(@(s) s{i},imScaledMaps,'unif',0);
            cutMin=min(allScaleMap{1}(allScaleMap{1}>0));
            cutMax=max(allScaleMap{1}(:));
            contVol=cutMin+mat2gray(vol)*(cutMax-cutMin);
            contScale=cutMin+mat2gray(maxResponseScale)*(cutMax-cutMin);
            allScaleMapReal=[contVol contScale allScaleMap{:}];
            bins=0:10:max(allScaleMapReal(:));
            allScaleMap=discretize(allScaleMapReal,bins);
            allScaleMap(isnan(allScaleMap))=0;
            [rgbMask2] = label2rgbND(allScaleMap+1, [0 0 0; parula(numel(bins)+1)] );
            imseriesmaskshowrgb([zeros(size(allScaleMapReal))],[rgbMask2]);
            % imseriesshow([contVol allScaleMap],'displayRange',[cutMin cutMax]);
        end
    else
        imseriesshow([imScaledMaps{:}],'displayRange',[min(imScaledMaps{1}(imScaledMaps{1}>0)) max(imScaledMaps{1}(:))]);
    end
end





toc;













%%%%%% OLD SNIPPETS


%% The scale volume is a map of scale. 
%% The scale is the largest scale that is a true local maximum in the scale space (function of scale)
% LoGPrevMax=LoGs{1};
% LoGLocaLMax=zeros(size(LoGPrevMax));
% scaleVol=ones(size(LoGPrevMax));
% scaleVol(LoGPrevMax==0)=0;
% for sIdx=2:(length(scales)-1)
%     % keep the mask of maximum response 
%     [maxInd]=LoGPrevMax<LoGs{sIdx};
%     LoGPrevMax(maxInd)=LoGs{sIdx}(maxInd);
%     
%     % Keep a mask of local maximum scale (function of scale evolution)
%     lmaxInd=(LoGs{sIdx}>LoGs{sIdx+1})&maxInd;   % is this local maximum ?
%     LoGLocaLMax(lmaxInd)=LoGs{sIdx}(lmaxInd);
% 
%     scaleVol(lmaxInd)=sIdx;
% end
% LoGmask=LoGLocaLMax;

% if(p.debug)
% %     imseriesmaskshowrgb(vol,{LMBU,rgbMask});
%     [rgbMask2] = label2rgbND( scaleVol+1, [0 0 0; parula(numel(p.scales)+2)] );
%     imseriesmaskshowrgb(vol,rgbMask2,'displayRange',prctile(vol(:),[20,99.999]));
% end

%%
% LMLoc=vertcat(LMLoc{:});
% LMLoG=vertcat(LMLoG{:});
% LMScaleIdx=vertcat(LMScaleIdx{:});

%%
% disp('Local scale maximum');tic;
% R = 10;
% D = createSparseDistanceMatrix(LMLoc, LMLoc, R);
% bestLMScaleIdx=zeros(size(LMScaleIdx));
% for lmIdx=1:size(D,1)
%     localLMIdx=find(D(lmIdx,:)~=0);
%     [~,bestLocalLMIdx]=max(LMLoG(localLMIdx));
%     bestLMScaleIdx(localLMIdx(bestLocalLMIdx))=LMScaleIdx(localLMIdx(bestLocalLMIdx));
% end
% clear LMLoG LMScaleIdx LMScaleIdx lmIdx D localLMIdx bestLocalLMIdx R
% toc

%%
% bestLMLoc=zeros(sum(bestLMScaleIdx~=0),3);
% bestLMLoc(:)=LMLoc(bestLMScaleIdx~=0,:);
% bestLMScaleIdx=bestLMScaleIdx(bestLMScaleIdx~=0);
% clear LMLoc



%%
% disp('buildLMvolume');tic;
% LM=ones(size(vol));
% LM(sub2ind(size(vol),bestLMLoc(:,1),bestLMLoc(:,2),bestLMLoc(:,3)))=bestLMScaleIdx+1;
% % clear bestLMScaleIdx bestLMLoc 

% 
% % Dilate LM
% if(p.version==1)
%     [x,y,z] = ndgrid(-4:4,-4:4,-3:3);
%     se = strel(sqrt(x.^2 + y.^2 + z.^2) <=4);
%     LM=imdilate(LM,se);
%     clear se x y z
% end
% 
% LMBU=[];
% if(p.debug)
%     [rgbMask] = label2rgbND( LM, [0 0 0; parula(numel(p.scales)+3)] );
%     LMBU=rgbMask;
%     imseriesmaskshowrgb(vol,rgbMask);
% end
%LM(LoGmask==0)=0;
