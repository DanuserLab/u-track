classdef TracksROI < DynROI
    properties (SetAccess = public, GetAccess = public)
    tracks;
    fringe;  
    end

    methods
        function obj = TracksROI(tracks,fringe,overlapping)
            if(nargin>1 )
                tracksCpy=tracks.copy();
                tr=tracks(1).copy();    

                %% Cut the tracks to minimum lifetime
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
                if((nargin<3)||overlapping)
                    for trIdx=2:numel(tracksCpy)
                        tr=tr.getOverlapping(tracksCpy(trIdx));
                    end
                    for trIdx=1:numel(tracksCpy)
                        tracksCpy(trIdx)=tracksCpy(trIdx).getOverlapping(tr);
                    end
                end

                obj.tracks=fillTrackGaps(tracksCpy);
                obj.fringe=fringe;  
            else
                obj.tracks=[];
                obj.fringe=0;
            end

            %% Default ref is the lab FoF centered around the starting position of the first track
            if(~isempty(obj.tracks))
                reftr=TracksHandle();
                reftr.startFrame=min([obj.tracks.startFrame]);
                reftr.endFrame=max([obj.tracks.endFrame]);
                reftr.x=obj.tracks(1).x(1)*ones(1,reftr.lifetime);
                reftr.y=obj.tracks(1).y(1)*ones(1,reftr.lifetime);
                reftr.z=obj.tracks(1).z(1)*ones(1,reftr.lifetime);
                ref=FrameOfRef().setOriginFromTrack(reftr).genCanonicalBase();
                obj.setDefaultRef(ref);
            end
        end

        function [idx,dist]=mapPosition(obj,positionMatrixCell)
            %% Frame Idx is ignore, we map on the full bounding box. Even if the object move.
            %% To map synchronously, map in the frame of ref.
            idx=cell(1,numel(positionMatrixCell));
            dist=cell(1,numel(positionMatrixCell));
            for pIdx=1:numel(positionMatrixCell)
                [minCoord, maxCoord]=obj.getBoundingBox([],pIdx);
                [idx{pIdx}, dist{pIdx}] = KDTreeRangeQuery(positionMatrixCell{pIdx}, ... 
                            [minCoord + maxCoord]/2,[maxCoord - minCoord]);
            end
        end

        function [mappedTracks,indices]=mapTracks(obj,tracks)
            [det,~,trackIndices]=Detections().getTracksCoord(tracks);
            [posIndices]=obj.mapPosition(det.getPosMatrix());
            indicesCell=cellfun(@(ti,pi) ti(pi),trackIndices,posIndices,'unif',0);
            indices=false(size(tracks));
            indices(unique(vertcat(indicesCell{:})))=true;
            mappedTracks=tracks(indices);
        end

        function [minCoord,maxCoord]=getBoundingBoxOptim(obj,ref,frameIdx)
            pos=Detections().getTracksCoord(obj.tracks);
            if(nargin>2)
                pos=pos(min(end,frameIdx));
            else
                frameIdx=1:numel(pos);
            end
            pos=pos.getPosMatrix();
            if(nargin>1)&&(~isempty(ref))
                parfor fIdx=1:numel(pos)
                    if(~isempty(pos{fIdx}))
                        pos{fIdx}=ref.applyBaseToPosPointCloud(pos{fIdx}(:,[1 2 3]),frameIdx(fIdx));
                    end
                end
            end
            pos=vertcat(pos{:});
            if(~isempty(pos))
                fringeWidth=obj.fringe;
                minXBorder=floor(min(pos(:,1)))-fringeWidth;
                minYBorder=floor(min(pos(:,2)))-fringeWidth;
                minZBorder=floor(min(pos(:,3)))-fringeWidth;
                maxXBorder=ceil(max(pos(:,1)))+fringeWidth;
                maxYBorder=ceil(max(pos(:,2)))+fringeWidth;
                maxZBorder=ceil(max(pos(:,3)))+fringeWidth;
            else
                [m,M]=obj.getBoundingBox(ref);
                minXBorder=m(1);
                minYBorder=m(2);
                minZBorder=m(3);
                maxXBorder=M(1); 
                maxYBorder=M(2);
                maxZBorder=M(3);
            end
            minCoord=[minXBorder minYBorder minZBorder];
            maxCoord=[maxXBorder maxYBorder maxZBorder];
        end
        
        function frame=getStartFrame(obj)
            frame=min([obj.tracks.startFrame]);
        end

        function frame=getEndFrame(obj)
            frame=max([obj.tracks.endFrame]);
        end

        function [minCoord, maxCoord]=getBoundingBox(obj,ref,frameIdx)
            if(nargin>1)&&(~isempty(ref))
                trs=ref.applyBase(obj.tracks);
            else
                trs=obj.tracks;
            end
            frameRange=[];
            if(nargin>2)
                frameRange=frameIdx;
            else
                frameRange=min([trs.startFrame]):max([trs.endFrame]);
            end
                        
            minX=[];
            minY=[];
            minZ=[];
            maxX=[]; 
            maxY=[];
            maxZ=[];
            for iP=1:length(trs)
                pIndices=ismember(trs(iP).f,frameRange);
                if(any(pIndices))
                    if(isempty(minX))
                        minX=floor((min(trs(iP).x(pIndices))));
                        minY=floor((min(trs(iP).y(pIndices))));
                        minZ=floor((min(trs(iP).z(pIndices))));
                        maxX=ceil((max(trs(iP).x(pIndices))));
                        maxY=ceil((max(trs(iP).y(pIndices))));
                        maxZ=ceil((max(trs(iP).z(pIndices))));
                    else
                        minX=floor(min(min(trs(iP).x(pIndices)),minX));
                        minY=floor(min(min(trs(iP).y(pIndices)),minY));
                        minZ=floor(min(min(trs(iP).z(pIndices)),minZ));
                        maxX=ceil(max(max(trs(iP).x(pIndices)),maxX));
                        maxY=ceil(max(max(trs(iP).y(pIndices)),maxY));
                        maxZ=ceil(max(max(trs(iP).z(pIndices)),maxZ));
                    end

                    % if(isempty(minX))
                    %     minX=floor((min(trs(iP).x(pIndices))));
                    %     minY=floor((min(trs(iP).y(pIndices))));
                    %     minZ=floor((min(trs(iP).z(pIndices))));
                    %     maxX=ceil((max(trs(iP).x(pIndices))));
                    %     maxY=ceil((max(trs(iP).y(pIndices))));
                    %     maxZ=ceil((max(trs(iP).z(pIndices))));
                    % else
                    %     minX=floor(min(min(trs(iP).x(pIndices)),minX));
                    %     minY=floor(min(min(trs(iP).y(pIndices)),minY));
                    %     minZ=floor(min(min(trs(iP).z(pIndices)),minZ));
                    %     maxX=ceil(max(max(trs(iP).x(pIndices)),maxX));
                    %     maxY=ceil(max(max(trs(iP).y(pIndices)),maxY));
                    %     maxZ=ceil(max(max(trs(iP).z(pIndices)),maxZ));
                    % end
                end
            end
            if(isempty(maxX))
                [m,M]=obj.getBoundingBox(ref);
                minXBorder=m(1);
                minYBorder=m(2);
                minZBorder=m(3);
                maxXBorder=M(1); 
                maxYBorder=M(2);
                maxZBorder=M(3);
            else
                fringeWidth=obj.fringe;
                maxXBorder=(maxX+fringeWidth);
                maxYBorder=(maxY+fringeWidth);
                maxZBorder=(maxZ+fringeWidth);
                minXBorder=(minX-fringeWidth);
                minYBorder=(minY-fringeWidth);
                minZBorder=(minZ-fringeWidth);
            end
            minCoord=[minXBorder minYBorder minZBorder];
            maxCoord=[maxXBorder maxYBorder maxZBorder];
        end

    end

end
