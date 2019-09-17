function [tracksXY,tracksZY,tracksZX]=overlayProjGraph(XYProj,ZYProj,ZXProj,XLimit,YLimit,ZLimit,positions,edges,myColormap,colorIndx,varargin)
ip = inputParser;
ip.CaseSensitive = false;
ip.KeepUnmatched = true;
ip.addOptional('cumulative',false);
ip.addOptional('detectionBorderDisplay',0);
ip.addParameter('showROIOnly',true); % show only the edges with a vertice that belong to the ROI.
ip.addParameter('positionsLabel',{}); 
ip.parse(varargin{:});
p=ip.Results;

  minXBorder=XLimit(1);
  maxXBorder=XLimit(2);
  minYBorder=YLimit(1);
  maxYBorder=YLimit(2);
  minZBorder=ZLimit(1);
  maxZBorder=ZLimit(2);

  % %% Print tracks on the projections
  % if(isempty(colorIndx))
  %   colorIndx=arrayfun(@(d) ones(1,size(d.zCoord,1)),detections,'unif',0);
  % end
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

  % if(~iscell(colorIndx))
  %   colorIndx={colorIndx};
  % end
  
  if(isempty(myColormap))
    myColormap=255*autumn(length(unique(vertcat(colorIndx{:}))));
  end
  if(p.showROIOnly)
    if(~isempty(positions))
      Z=positions(:,3);
      Y=positions(:,2);
      X=positions(:,1);
      keepIndx= (Z>=minZBorder-p.detectionBorderDisplay)&(Z<=maxZBorder+p.detectionBorderDisplay)& ...
                (X>=minXBorder-p.detectionBorderDisplay)&(X<=maxXBorder+p.detectionBorderDisplay)& ...
                (Y>=minYBorder-p.detectionBorderDisplay)&(Y<=maxYBorder+p.detectionBorderDisplay);
      keepIndx=find(keepIndx);
      keepEdge=ismember(edges(:,1),keepIndx)|ismember(edges(:,2),keepIndx);
      colorIndx=colorIndx(keepEdge);
      edges=edges(keepEdge,:);
    end
  end


  if(~isempty(positions))
    tracksXY=graphBinaryOverlay(XYProj,[minXBorder maxXBorder],[minYBorder maxYBorder], ... 
                                positions(:,1:2),edges,colorIndx,myColormap,varargin{:},'positionsLabel',p.positionsLabel);
  else
    tracksXY=XYProj;
  end

  if(~isempty(positions))
    tracksZY=graphBinaryOverlay(ZYProj,[minZBorder maxZBorder],[minYBorder maxYBorder], ... 
                                positions(:,[3 2]),edges,colorIndx,myColormap,varargin{:},'positionsLabel',p.positionsLabel);
  else
    tracksZY=ZYProj;
  end;

  if(~isempty(positions))
    tracksZX=graphBinaryOverlay(ZXProj,[minZBorder maxZBorder],[minXBorder maxXBorder], ...
                                positions(:,[3 1]),edges,colorIndx,myColormap,varargin{:},'positionsLabel',p.positionsLabel);
  else
    tracksZX=ZXProj;
  end;
