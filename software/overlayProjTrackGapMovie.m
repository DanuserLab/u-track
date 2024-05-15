function overlayProjTrackGapMovie(processProj,varargin)
ip = inputParser;
ip.CaseSensitive = false;
ip.KeepUnmatched = true;
ip.addOptional('tracks',[]);
ip.parse(processProj,varargin{:});
p=ip.Results;

tracks=p.tracks;
gapDetections=collectTrackGaps(tracks);
overlayProjDetectionMovie(processProj,'detections',gapDetections,varargin{:},'showNumber',false); 


function gapDetections=collectTrackGaps(tracks)
  pos=nan(3,numel(tracks),max([tracks.endFrame]));
  se=[zeros(1,tracks.numTimePoints()) 1 ones(1,tracks.numTimePoints())];

  for tIdx=1:numel(tracks)
    F=tracks(tIdx).f;
    gi=isnan(tracks(tIdx).x);    
    if(any(gi))    
    copyIdx=1:tracks(tIdx).lifetime;
    copyIdx(gi)=0;
    copyIdx=imdilate(copyIdx,se);
    copyIdx=copyIdx(gi);
    pos(1,tIdx,F(gi))=tracks(tIdx).x(copyIdx);
    pos(2,tIdx,F(gi))=tracks(tIdx).y(copyIdx);
    pos(3,tIdx,F(gi))=tracks(tIdx).z(copyIdx);
    end
  end
  posCell=arrayfun(@(f) pos(:,~isnan(pos(1,:,f)),f)',1:max([tracks.endFrame]),'unif',0);
  gapDetections=Detections().initFromPosMatrices(posCell,posCell);


function gapDetections=collectTrackSOE(tracks)
    pos=nan(3,numel(tracks),max([tracks.endFrame]));
    for tIdx=1:numel(tracks)
      F=tracks(tIdx).f;
      gi=isnan(tracks(tIdx).x);        
      pos(1,tIdx,F(gi))=tracks(tIdx).x(gi);
      pos(2,tIdx,F(gi))=tracks(tIdx).y(gi);
      pos(3,tIdx,F(gi))=tracks(tIdx).z(gi);
    end
    posCell=arrayfun(@(f) pos(:,~isnan(pos(1,:,f)),f)',1:max([tracks.endFrame]),'unif',0);
    gapDetections=Detections().initFromPosMatrices(posCell,posCell);
%
% Copyright (C) 2024, Danuser Lab - UTSouthwestern 
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
