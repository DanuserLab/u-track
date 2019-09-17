function [tracksXY,tracksZY,tracksZX,keepIndx]=overlayProjDetections(XYProj,ZYProj,ZXProj,XLimit,YLimit,ZLimit,detections,myColormap,colorIndx,varargin)
ip = inputParser;
ip.CaseSensitive = false;
ip.KeepUnmatched = true;
ip.addOptional('cumulative',false);
ip.addOptional('detectionBorderDisplay',0);
ip.addOptional('positionsLabel',[]);
ip.addOptional('radius',2);
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
    colorIndx=arrayfun(@(d) ones(1,size(d.zCoord,1)),detections,'unif',0);
  end

  if(~iscell(colorIndx))
    colorIndx={colorIndx};
  end
  
  radius=p.radius;
  if(~iscell(radius))
    radius={radius};
  end
  if(isempty(myColormap))
    myColormap=255*autumn(length(unique(vertcat(colorIndx{:}))));
  end

  positionsLabel=p.positionsLabel;

  if(~isempty(detections))
      keepIndx=cell(1,length(detections));
      for dIdx=1:length(detections)
          d=detections(dIdx);
          if(~isempty(d.xCoord))
              
              keepIndx{dIdx}= (d.zCoord(:,1)>=minZBorder-p.detectionBorderDisplay)&(d.zCoord(:,1)<=maxZBorder+p.detectionBorderDisplay)& ...
                              (d.xCoord(:,1)>=minXBorder-p.detectionBorderDisplay)&(d.xCoord(:,1)<=maxXBorder+p.detectionBorderDisplay)& ...
                              (d.yCoord(:,1)>=minYBorder-p.detectionBorderDisplay)&(d.yCoord(:,1)<=maxYBorder+p.detectionBorderDisplay);
          end

          if(~isempty(positionsLabel))
            positionsLabel{dIdx}=positionsLabel{dIdx}(keepIndx{dIdx});
          end
          colorIndx{dIdx}=colorIndx{dIdx}(keepIndx{dIdx});
          if(numel(radius{dIdx})>1)
                radius{dIdx}=radius{dIdx}(keepIndx{dIdx});
          end
      end
      detections=detections.copy().selectIdx(keepIndx);
  end


  if(~isempty(detections))
      % Only Keep detections within ZLimit
    tracksXY=detectionBinaryOverlay(XYProj,[minXBorder maxXBorder],[minYBorder maxYBorder],detections,colorIndx,myColormap,varargin{:}, ...
                                    'radius',radius,'positionsLabel',positionsLabel);
  else
    tracksXY=XYProj;
  end

  if(~isempty(detections))
      numDet=length(detections);
      trdetections(numDet)=struct('xCoord',[],'yCoord',[]);
      for tIdx=1:numDet
          trdetections(tIdx).xCoord=detections(tIdx).zCoord;
          trdetections(tIdx).yCoord=detections(tIdx).yCoord;
      end
    
    tracksZY=detectionBinaryOverlay(ZYProj,[minZBorder maxZBorder],[minYBorder maxYBorder],trdetections,colorIndx,myColormap,varargin{:}, ...
                                    'radius',radius,'positionsLabel',positionsLabel);
  else
    tracksZY=ZYProj;
  end;


  if(~isempty(detections))
      for tIdx=1:numDet
          trdetections(tIdx).yCoord=detections(tIdx).xCoord;
      end

%     capturedEB3ZX=capturedEB3ZY.copy();
%     for tIdx=1:length(capturedEB3ZX)
%       capturedEB3ZX(tIdx).y=tracksInMask(tIdx).x;
%     end
    tracksZX=detectionBinaryOverlay(ZXProj,[minZBorder maxZBorder],[minXBorder maxXBorder],trdetections,colorIndx,myColormap,varargin{:}, ...
                                    'radius',radius,'positionsLabel',positionsLabel);
  else
    tracksZX=ZXProj;
  end;
