function [tracksXY,tracksZY,tracksZX]=overlayProjTracks(XYProj,ZYProj,ZXProj,XLimit,YLimit,ZLimit,fIdx,tracksInMask,myColormap,colorIndx,varargin)
  ip = inputParser;
  ip.CaseSensitive = false;
  ip.KeepUnmatched = true;
  ip.addOptional('cumulative',false);
  ip.addOptional('detectionBorderDisplay',0);
  ip.parse(varargin{:});
  p=ip.Results;

  minXBorder=XLimit(1);
  maxXBorder=XLimit(2);
  minYBorder=YLimit(1);
  maxYBorder=YLimit(2);
  minZBorder=ZLimit(1);
  maxZBorder=ZLimit(2);

  %% Print tracks on the projections
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

  if(isempty(colorIndx))
%    colorIndx=floor(linspace(1,255,length(tracksInMask)));
     colorIndx=ones(1,length(tracksInMask));
  end
  
  keepIdx=[];
  if(~isempty(tracksInMask))
      keepIdx=false(1,length(tracksInMask));
        for tIdx=1:length(tracksInMask)
          tr=tracksInMask(tIdx);
          coordIdx=(tr.f==fIdx);
          if(any(coordIdx))
            X=tr.x(coordIdx);Y=tr.y(coordIdx); Z=tr.z(coordIdx);
            keepIdx(tIdx)= (Z>=minZBorder-p.detectionBorderDisplay)&&(Z<=maxZBorder+p.detectionBorderDisplay) ... 
                          &&(X>=minXBorder-p.detectionBorderDisplay)&&(X<=maxXBorder+p.detectionBorderDisplay) ...
                          &&(Y>=minYBorder-p.detectionBorderDisplay)&&(Y<=maxYBorder+p.detectionBorderDisplay);
          end
      end
      tracksInMask=tracksInMask(keepIdx);
      colorIndx=colorIndx(keepIdx);
  end

  if(~isempty(tracksInMask))
    tracksXY=trackBinaryOverlay(XYProj,[minXBorder maxXBorder],[minYBorder maxYBorder],tracksInMask,fIdx,colorIndx,myColormap,'trackIDs',find(keepIdx),varargin{:});
  else
    tracksXY=XYProj;
  end

  if(~isempty(tracksInMask))
      numTrack=length(tracksInMask);
      transTracks(numTrack)=struct('x',[],'y',[],'f',[]);
      for tIdx=1:numTrack
          transTracks(tIdx).x=tracksInMask(tIdx).z;
          transTracks(tIdx).y=tracksInMask(tIdx).y;
          transTracks(tIdx).f=tracksInMask(tIdx).f;
      end

%     capturedEB3ZY=tracksInMask.copy();
%     for tIdx=1:length(capturedEB3ZY)
%       capturedEB3ZY(tIdx).x=tracksInMask(tIdx).z ;%*MD.pixelSize_/MD.pixelSizeZ_;
%     end
    tracksZY=trackBinaryOverlay(ZYProj,[minZBorder maxZBorder],[minYBorder maxYBorder],transTracks,fIdx,colorIndx,myColormap,'trackIDs',find(keepIdx),varargin{:});
  else
    tracksZY=ZYProj;
  end;


  if(~isempty(tracksInMask))
      for tIdx=1:numTrack
          transTracks(tIdx).y=tracksInMask(tIdx).x;
      end

%     capturedEB3ZX=capturedEB3ZY.copy();
%     for tIdx=1:length(capturedEB3ZX)
%       capturedEB3ZX(tIdx).y=tracksInMask(tIdx).x;
%     end
    tracksZX=trackBinaryOverlay(ZXProj,[minZBorder maxZBorder],[minXBorder maxXBorder],transTracks,fIdx,colorIndx,myColormap,'trackIDs',find(keepIdx),varargin{:});
  else
    tracksZX=ZXProj;
  end;
