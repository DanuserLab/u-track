function keptIdx=overlayProjDetectionMovie(processProj,varargin)
ip = inputParser;
ip.CaseSensitive = false;
ip.KeepUnmatched = true;
ip.addRequired('processProj');
ip.addOptional('detections',[]);
ip.addOptional('colormap',[]);
ip.addOptional('process',[]);
ip.addOptional('processFrames',[]);
ip.addOptional('cumulative',false);
ip.addOptional('saveVideo',false);
ip.addOptional('printVectorFilePattern','');
ip.addOptional('detectionFrameIdx',[]);  % Useful for decimation: specify the frame associated to each detection 
ip.addOptional('colorIndx',[],@iscell);
ip.addOptional('colorLabel',[],@iscell);
ip.addOptional('minMaxLabel',[],@isnumeric);
ip.addOptional('radius',2);
ip.addOptional('name','detections');
ip.addOptional('showNumber',false);
ip.parse(processProj,varargin{:});
p=ip.Results;

detections=p.detections;
cumulative=p.cumulative;
%% testing imwarp to crop the image
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

projDataIdx=5;
ref=[];


if(isa(processProj,'ExternalProcess'))
  processProjDynROI=ProjectDynROIRendering();
  processProjDynROI.importFromDeprecatedExternalProcess(processProj);
  try % Handle Project1D/ProjDyn different outFilePaths_spec (need to be defined through a class...)
    projData=load(processProj.outFilePaths_{projDataIdx},'minXBorder', 'maxXBorder','minYBorder','maxYBorder','minZBorder','maxZBorder','frameNb');
  catch
    projDataIdx=4;
  end
  projData=load(processProj.outFilePaths_{projDataIdx},'minXBorder', 'maxXBorder','minYBorder','maxYBorder','minZBorder','maxZBorder','frameNb');
  processProjDynROI.setBoundingBox( [projData.minXBorder projData.maxXBorder], [projData.minYBorder projData.maxYBorder], [projData.minZBorder projData.maxZBorder]);
  processProj=processProjDynROI;
end


ref=get(processProj,'ref');
if(~isempty(ref))
    % Ugly hack to fix lack of detection frameID support
    if(~isempty(p.detectionFrameIdx))
        D(max(p.detectionFrameIdx))=Detections();
        if(numel(detections)<numel(p.detectionFrameIdx))
          detections(numel(p.detectionFrameIdx))=Detections();
        end
        D(p.detectionFrameIdx)=detections;
        detections=D;
    end    
    detections=ref.applyBase(detections,'')
    if(~isempty(p.detectionFrameIdx))
        detections=detections(p.detectionFrameIdx);
    end    
end
projData=processProj;


frameNb=min([projData.frameNb,length(detections)]);
processFrames=p.processFrames;
if(isempty(processFrames))
    processFrames=processProj.getStartFrame():processProj.getEndFrame();
end

frameNb=min([projData.frameNb,length(detections),length(processFrames)]);


%% create projection process saving independant projection location
% if(~isempty(p.process))
%   processRenderer = ProjRendering(processProj,p.name);
% end

colorIndx=p.colorIndx;

colormapSize=size(p.colormap,1);
if(isempty(p.colormap))
    colormapSize=256;
end

if(~isempty(p.colorLabel))

  if(isempty(p.minMaxLabel))
    try
      allLabel=vertcat(p.colorLabel{:});
    catch
      allLabel=[p.colorLabel{:}];
    end;
    minLabel=min(allLabel);
    maxLabel=max(allLabel);
  else
    minLabel=p.minMaxLabel(1);
    maxLabel=p.minMaxLabel(2);    
  end   
  try
   colorIndx=cellfun(@(d) ceil((colormapSize-1)*mat2gray(reshape(d,numel(d),1),[minLabel,maxLabel]))+1,p.colorLabel,'unif',0);
 catch
  colorIndx=cellfun(@(d) ones(size(p.colorLabel)),p.colorLabel,'unif',0);
  end
end

if(isempty(colorIndx))
    colorIndx=arrayfun(@(d) ones(1,size(d.zCoord,1)),detections,'unif',0);
end

acolormap=p.colormap;
if(isempty(acolormap))
    acolormap=255*hsv(colormapSize);
end
  
keptIdx=cell(1,numel(processFrames));

%% Going around parfor limitation 
%% TODO: slicing project Class and use spmd for more control
pCell1=cell(1,numel(processFrames));
pCell2=cell(1,numel(processFrames));
pCell3=cell(1,numel(processFrames));
oCell1=cell(1,numel(processFrames));
oCell2=cell(1,numel(processFrames));
oCell3=cell(1,numel(processFrames));
for fIdx=processFrames
  [pCell1{fIdx},pCell2{fIdx},pCell3{fIdx}]=processProj.loadFrame(1,fIdx);
end

[xBound,yBound,zBound]=projData.getBoundingBox();

radii=p.radius;
detFrame=p.detectionFrameIdx;
if(isempty(detFrame))
  detFrame=1:numel(detections);
end

for fIdx=processFrames
    XYProj=pCell1{fIdx};
    ZYProj=pCell2{fIdx};
    ZXProj=pCell3{fIdx};
    fRadii=radii;
    if(cumulative)
        detectionsAtFrame=detections;
        fColorIndx=colorIndx;
    else
        dIdx=find(fIdx==detFrame);
        if(~isempty(dIdx))
        detectionsAtFrame=detections(dIdx);
        fColorIndx=colorIndx{dIdx};
        if(iscell(radii))
          fRadii=radii{dIdx};   
        end
        else
          fColorIndx=[];
          fRadii=[];
          detectionsAtFrame=Detections();
        end
    end
    printVectorFile=[];
    if(~isempty(p.printVectorFilePattern))
        mkdirRobust(fileparts(p.printVectorFilePattern));
        printVectorFile=sprintfPath(p.printVectorFilePattern,fIdx);
    end

    positionsLabel=[];
    if(p.showNumber)
      N=sum(detections.getCard());
      positionsLabel={arrayfun(@(n) num2str(n),1:N,'unif',0)};
    end
    % detectionsAtFrame.getAllStruct()
    % detectionsAtFrame.zCoord(:,1)=detectionsAtFrame.zCoord(:,1)/0.378;
    [overlayXY,overlayZY,overlayZX,keptIdx(fIdx)]=overlayProjDetections( ... 
        XYProj,ZYProj,ZXProj,xBound,yBound,zBound, ...
        detectionsAtFrame,acolormap,fColorIndx, ...
        'positionsLabel',positionsLabel,'printVectorFile',printVectorFile,varargin{:},'radius',fRadii);

    oCell1{fIdx}=overlayXY;
    oCell2{fIdx}=overlayZY;
    oCell3{fIdx}=overlayZX;
end

if(~isempty(p.process))
  for fIdx=processFrames
      p.process.saveFrame(1,fIdx,oCell1{fIdx},oCell2{fIdx},oCell3{fIdx});
  end
end

if(~isempty(p.process)&&p.saveVideo) 
    ProjAnimation(p.process,'ortho').saveVideo([p.process.getOutputDir()  '.avi']);
end

