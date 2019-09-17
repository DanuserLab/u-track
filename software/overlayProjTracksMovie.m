function projImages=overlayProjTracksMovie(processProj,varargin)
  ip = inputParser;
  ip.CaseSensitive = false;
  ip.KeepUnmatched = true;
  ip.addRequired('processProj');
  ip.addParameter('tracksOrProcess',[]);
  ip.addParameter('tracks',[]);
  ip.addParameter('process',[]);
  ip.addParameter('processFrames',[]);
  ip.addParameter('showNumber',false);
  ip.addParameter('show',true);
  ip.addParameter('dragonTail',[]);
  ip.addParameter('colormap',[]);
  ip.addParameter('colorIndx',[]);
  ip.addParameter('colorLabel',[]);
  ip.addParameter('minMaxLabel',[]);
  ip.addParameter('detectionFrameIdx',[]);  % Useful for decimation: specify the frame associated to each detection 
  ip.addParameter('saveVideo',false);
  ip.addParameter('useGraph',false);
  ip.addParameter('name','tracks');
  ip.parse(processProj,varargin{:});
  p=ip.Results;

  if(isempty(p.tracks))
    tracks=p.tracksOrProcess;
  else
    tracks=p.tracks;
  end
  
  if(isa(tracks,'Process'))
      tracksFinal=tracks.loadChannelOutput(1);
      tracks=TracksHandle(tracksFinal);
  end
      
  if(isempty(tracks))
      return;
  end
  
  processFrames=p.processFrames;
  if(isempty(processFrames))
      processFrames=processProj.getStartFrame():processProj.getEndFrame();
  end

  detFrame=p.detectionFrameIdx;
  if(isempty(detFrame))
    % detFrame=1:max(processFrames);
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
    detFrame=processFrames;
  end


  projDataIdx=5; 
  ref=[];
  % try % Handle Project1D/ProjDyn with different implementation (current and then deprecated)
    ref=get(processProj,'ref');
    if(~isempty(ref)&&~isempty(p.detectionFrameIdx))
      refDecim=copy(ref);
      refDecim=refDecim.selectFrame(detFrame);
      refDecim.frame=1:numel(refDecim.frame);
      ref=refDecim;
    end
    if(~isempty(ref))
    tracks=ref.applyBase(tracks,'');
    end
    projData=processProj;
  % catch
  %   try % Handle Project1D/ProjDyn different outFilePaths_spec (need to be defined through a class...)
  %     projData=load(processProj.outFilePaths_{projDataIdx},'minXBorder', 'maxXBorder','minYBorder','maxYBorder','minZBorder','maxZBorder','frameNb');
  %   catch
  %     projDataIdx=4;
  %   end
  %   projData=load(processProj.outFilePaths_{projDataIdx},'minXBorder', 'maxXBorder','minYBorder','maxYBorder','minZBorder','maxZBorder','frameNb');
  % end


%% create projection process saving independant projection location
% if(~isempty(p.process))
% end
colorIndx=p.colorIndx;

if(isempty(colorIndx))
    colorIndx=1:numel(tracks);
end

myColormap=p.colormap;
if(isempty(myColormap))
    myColormap=255*jet(256);
end

colormapSize=size(myColormap,1);
if(~isempty(p.colorLabel))
  if(iscell(p.colorLabel))
    allLabel=cellfun(@(c) c(:),p.colorLabel,'unif',0);
    allLabel=vertcat(allLabel{:});
  else
    allLabel=p.colorLabel;  
  end

   if(isempty(p.minMaxLabel))
     minLabel=min(allLabel);
     maxLabel=max(allLabel);
   else
     minLabel=p.minMaxLabel(1);
     maxLabel=p.minMaxLabel(2);    
   end
   if(iscell(p.colorLabel))
     colorIndx=cellfun(@(a) ceil((colormapSize-1)*mat2gray(a,[minLabel,maxLabel]))+1,p.colorLabel,'unif',0);
  else
     colorIndx=ceil((colormapSize-1)*mat2gray(allLabel,[minLabel,maxLabel]))+1;
  end
end
 assert(numel(colorIndx)==numel(tracks));

saveInProcess=~isempty(p.process);


[xBound,yBound,zBound]=projData.getBoundingBox();

% [det,~,trackIndices]=Detections().getTracksCoord(tracks);
% detStruct=det.getAllStruct();
% detStruct.trackID=vertcat(trackIndices{:});

trackFrame=arrayfun(@(t) t.f,tracks,'unif',0);

frameNb=projData.frameNb;


% A attempt at a dirty trick to avoid copying p.process across parallel worker
imLoader= @(fIdx)(processProj.loadFrame(1,fIdx));
imSaved=cell(1,numel(processFrames));


parfor fIdxIdx=1:numel(processFrames)
  fIdx=processFrames(fIdxIdx);
  [XYProj,ZYProj,ZXProj,three]=feval(imLoader,fIdx);
  % dIdx=find(fIdx==detFrame);
  dIdx=fIdx;



  if(~p.useGraph)
    [tracksXY,tracksZY,tracksZX]=overlayProjTracks(XYProj,ZYProj,ZXProj, ...
      [projData.minXBorder projData.maxXBorder],...
      [projData.minYBorder projData.maxYBorder],...
      [projData.minZBorder projData.maxZBorder],...
        dIdx,tracks,myColormap,colorIndx,varargin{:});
  else
    %% All positions and edge representing track present on frame dIdx
    tracksInFrame=cellfun(@(f) any(f==dIdx),trackFrame);
    [vert,edges,frames,edgesLabel]=tracks(tracksInFrame).getGraph();


    %% position label only kept for all edges on the current view
    positionsLabel={};
    if(p.showNumber)
      positionsLabel=cell(1,size(vert,1));
      for pIdx=find(frames==dIdx)'
        tracksInFrameNum=[tracks(tracksInFrame).index];
        positionsLabel{edges(pIdx,2)}=num2str(tracksInFrameNum(edgesLabel(pIdx)));
      end
    end

    %% keep only the edge before the current time points
    tailIndx=(frames<=dIdx)&(frames>=(dIdx-p.dragonTail));
    edges=edges(tailIndx,:);

    %% sync colorIndx
    fColorIndx=colorIndx(tracksInFrame);

    if(iscell(fColorIndx))
      fColorIndx=cellfun(@(c) c(:),fColorIndx,'unif',0);
      % fColorLabel=arrayfun(@(l) l*ones(numel(fColorIndx{l}),1),1:numel(fColorIndx),'unif',0);
      fColorIndx=vertcat(fColorIndx{:});
      fColorIndx=fColorIndx(tailIndx);
    else
      fColorIndx=fColorIndx(edgesLabel(tailIndx));
    end
    assert(numel(fColorIndx)==size(edges,1));

    [tracksXY,tracksZY,tracksZX]=overlayProjGraph(XYProj,ZYProj,ZXProj, ...
      xBound,yBound,zBound,vert,edges, ... 
      myColormap,fColorIndx,'positionsLabel',positionsLabel,varargin{:});
  end


   if(saveInProcess)
    imSaved{fIdxIdx}={tracksXY,tracksZY,tracksZX};
   end
  
end

for fIdxIdx=1:numel(processFrames)
    fIdx=processFrames(fIdxIdx);
  p.process.saveFrame(1,fIdx,imSaved{fIdxIdx}{:});
end

if(~isempty(p.process)&&p.saveVideo) 
    ProjAnimation(p.process,'ortho').saveVideo([p.process.getOutputDir()  '.avi']);
end



