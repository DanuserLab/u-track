classdef FrameOfRef < handle  & matlab.mixin.Copyable
   properties
      origin;
      X; % N*3 double
      Y; % N*3 double
      Z; % N*3 double
      frame; % N*1 uint16 the frames index during which the FrameOfReference is defined
   end
   methods
      function obj=setOriginFromTrack(obj,tr)
        tr=fillTrackGaps(tr);
         obj.origin=[tr.x' tr.y' tr.z'];
         obj.frame=tr.f;
      end

      function obj=setZFromTrack(obj,tr)
         % Only the overlapping frame can be kept
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
         [F,idxTr,idxObj] = intersect(tr.f,obj.frame);
         obj.frame=F;
         obj.origin=obj.origin(idxObj,:);
         obj.Z=[tr.x(idxTr)' tr.y(idxTr)' tr.z(idxTr)']-obj.origin;
         obj.Z=obj.Z./repmat(sum(obj.Z.^2,2).^0.5,1,3);
      end

      function obj=setZFromTrackUnion(obj,tr)
         % Extrapolate value to cover the cumulated lifetime
         [F,idxObj,idxTr] = union(obj.frame,tr.f);

         origin=zeros(max(F)-min(F)+1,3);
         origin(ismember(F,obj.frame),:)=obj.origin;
         origin(F<min(obj.frame),:)=repmat(obj.origin(1,:),[sum(F<min(obj.frame)),1]);
         origin(F>max(obj.frame),:)=repmat(obj.origin(end,:),[sum(F>max(obj.frame)),1]);
         obj.origin=origin;

         M=[tr.x' tr.y' tr.z'];
         Z=zeros(max(F)-min(F)+1,3);
         Z(ismember(F,tr.f),:)=M;
         Z(F<min(tr.f),:)=repmat(M(1,:),[sum(F<min(tr.f)),1]);
         Z(F>max(tr.f),:)=repmat(M(end,:),[sum(F>max(tr.f)),1]);
         obj.Z=Z-obj.origin;
         obj.Z=obj.Z./repmat(sum(obj.Z.^2,2).^0.5,1,3);
         obj.frame=min(F):max(F);
      end

      function n=getFrameNumber(obj)
        n=numel(obj.frame);
      end
      function obj=genBaseFromZ(obj,trX)
          if(nargin<2)
            obj.X=[0*obj.Z(:,1),obj.Z(:,3),-obj.Z(:,2)];
            obj.X=obj.X./repmat(sum(obj.X.^2,2).^0.5,1,3);
          else
              [F,idxTr,idxObj] = intersect(trX.f,obj.frame);
              obj.frame=F;
              obj.origin=obj.origin(idxObj,:);
              obj.Z=obj.Z(idxObj,:);
              
              obj.X=[trX.x(idxTr)' trX.y(idxTr)' trX.z(idxTr)']-obj.origin;

              obj.X=obj.X./repmat(sum(obj.X.^2,2).^0.5,1,3);
              obj.X=obj.X-repmat(sum(obj.X.*obj.Z,2),1,3).*obj.Z;
              obj.X=obj.X./repmat(sum(obj.X.^2,2).^0.5,1,3);
          end
          obj.Y=cross(obj.Z,obj.X);
          obj
      end

      function obj=selectFrame(obj,selectFrameIdx)
         obj.frame=obj.frame(selectFrameIdx);
         obj.origin=obj.origin(selectFrameIdx,:);
         obj.X=obj.X(selectFrameIdx,:);
         obj.Y=obj.Y(selectFrameIdx,:);
         obj.Z=obj.Z(selectFrameIdx,:);
      end
      
      function obj=genCanonicalBase(obj)
          obj.X=repmat([1 0 0],[size(obj.origin,1) 1]);
          obj.Y=repmat([0 1 0],[size(obj.origin,1) 1]);
          obj.Z=repmat([0 0 1],[size(obj.origin,1) 1]);
      end
      function obj=genCanonicalRef(obj,frameNb)
          obj.origin=zeros(frameNb,3);
          obj.frame=1:frameNb;
          obj.genCanonicalBase();
      end
      function tracks=getTracksFromBaseVector(obj,direction)
        vector=[]
        switch direction
        case 'X'
          vector=obj.X;
        case 'Y'
          vector=obj.Y;
        case 'Z'
          vector=obj.Z;
        otherwise
          error('Direction must be X,Y or Z')
        end
        tracks=TracksHandle();
        tracks.x=vector(:,1)';
        tracks.y=vector(:,2)';
        tracks.z=vector(:,3)';
        tracks.startFrame=obj.frame(1);
        tracks.endFrame=obj.frame(end);
      end
      function newBaseObject=applyBase(obj,tracksOrDetections,name)
        if(nargin<3)
          name='';
        end
        if(isa(tracksOrDetections,'Tracks'))
            newBaseObject=applyBaseToTrack(obj,tracksOrDetections,name);
        else
            newBaseObject=applyBaseToDetection(obj,tracksOrDetections,name);
        end
      end

      function newBaseObject=applyInvBase(obj,tracksOrDetections,name)
        if(nargin<3)
          name='';
        end
        if(isa(tracksOrDetections,'Tracks'))
            newBaseObject=applyInvBaseToTrack(obj,tracksOrDetections,name);
        else
            newBaseObject=applyInvBaseToDetection(obj,tracksOrDetections,name);
        end
      end

      function tracksBase=applyBaseToTrack(obj,tracks,name)
          tracksBase=tracks.copy();
          try
              if(nargin>2)&&(~isempty(name))
                  tracks.addprop(name);
              end
          catch
          end;
          
%           [pos,~,tracksID]=Detections().getTracksCoord(tracks);
%           posCell=pos.getPosMatrix();
%           for fIdx=1:numel(pos)
%             posCell{fIdx}=obj.applyBaseToPosPointCloud(posCell{fIdx},fIdx);
%           end
%           pos.setPosMatrix(posCell);
%           tracksBase=pos.buildTracksFromDetection(tracksID);
          
          for trIdx=1:length(tracks)
              tr=tracks(trIdx);
              trBase=tracksBase(trIdx);

          end
          for trIdx=1:length(tracks)
              tr=tracks(trIdx);
              trBase=tracksBase(trIdx);
              % Copying EB3 track
              trBase.addprop('ref');
              trBase.ref=obj;

              % Register in original tr
              try
                  trBase.addprop('originalRef');
              catch
              end;
              if(nargin>2)&&(~isempty(name))
                  setfield(tr,name,trBase);
              end;
              trBase.originalRef=tr;
              for pIdx=1:length(tr.f)
%                       V=obj.applyBaseToPosPointCloud([tr.x(pIdx) tr.y(pIdx) tr.z(pIdx)],tr.f(pIdx));
%                       trBase.x(pIdx)=V(1);
%                       trBase.y(pIdx)=V(2);
%                       trBase.z(pIdx)=V(3);
                  f=tr.f(pIdx);
                  B=obj.getBase(f);
                  orig=obj.getOrigAtFrame(f);
                  recentered=[(tr.x(pIdx)-orig(1)) (tr.y(pIdx)-orig(2)) (tr.z(pIdx)-orig(3))];
                  v=recentered*B;
                  trBase.x(pIdx)=v(1); trBase.y(pIdx)=v(2); trBase.z(pIdx)= v(3);
               end;
          end
        end

        function tracksBase=applyInvBaseToTrack(obj,tracks,name)
            tracksBase=tracks.copy();
            try
                if(nargin>2)&&(~isempty(name))
                    tracks.addprop(name);
                end
            catch
            end;
            for trIdx=1:length(tracks)
                tr=tracks(trIdx);
                trBase=tracksBase(trIdx);
                % Copying EB3 track
                trBase.addprop('ref');
                trBase.ref=obj;

                % Register in original tr
                try
                    trBase.addprop('originalRef');
                catch
                end;
                if(nargin>2)&&(~isempty(name))
                    setfield(tr,name,trBase);
                end;
                trBase.originalRef=tr;               
                for pIdx=1:length(tr.f)

                    f=tr.f(pIdx);
                    B=obj.getBase(f);
                    orig=obj.getOrigAtFrame(f);
                    v=[tr.x(pIdx) tr.y(pIdx) tr.z(pIdx)]*B^-1;
                    trBase.x(pIdx)=v(1)+orig(1); trBase.y(pIdx)=v(2)+orig(2); trBase.z(pIdx)= v(3)+orig(3);
                end;
            end
          end

        function B= getBase(obj,f)
          pIdx=find(obj.frame==f,1);
          if(isempty(pIdx))
              if(f>max(obj.frame))   pIdx=length(obj.frame);  else   pIdx=1; end;
          end
          B=[obj.X(pIdx,:)' obj.Y(pIdx,:)' obj.Z(pIdx,:)'];
        end

        function obj= setBase(obj,B)
          if(~iscell(B))
            B={B};
          end
          nBase=numel(B);
          for f=obj.frame
            bf=min(f,nBase);
            obj.X(f,:)=B{bf}(:,1);
            obj.Y(f,:)=B{bf}(:,2);
            obj.Z(f,:)=B{bf}(:,3);
          end
        end

      function orig= getOrigAtFrame(obj,f)
          pIdx=find(obj.frame==f,1);
          if(isempty(pIdx)); if(f<obj.frame); pIdx=1; else pIdx=length(obj.frame); end; end;
          orig=obj.origin(pIdx,:);
      end

      function detectionsBase=applyBaseToDetection(obj,detections,name)
          if(isempty(detections))
              detectionsBase=[];
              return;
          end
          detectionsBase=detections.copy();
          try
              if(nargin>2)&&(~isempty(name))
                  detections.addprop(name);
              end
          catch
          end;
          % for fIdx=1:length(detections)
          %     detect=detections(fIdx);
          %     detectBase=detectionsBase(fIdx);
          %     % Copying EB3 track
          %     detectBase.ref=obj;
          %     B=obj.getBase(fIdx);
          %     O=obj.getOrigAtFrame(fIdx);
          %     % easily optimized as implementend in poleDist
          %     for pIdx=1:size(detectBase.xCoord,1)
          %         recentered=[(detect.xCoord(pIdx,1)-O(1)) (detect.yCoord(pIdx,1)-O(2)) (detect.zCoord(pIdx,1)-O(3))];
          %         v=recentered*B;
          %         detectBase.xCoord(pIdx,1)=v(1); detectBase.yCoord(pIdx,1)=v(2); detectBase.zCoord(pIdx,1)= v(3);
          %     end;
          %     if(nargin>2)&&(~isempty(name))
          %         setfield(detect,name,detectBase);
          %     end;
          % end

          pos=detections.getPosMatrix();
          err=detections.getErrMatrix();
          for fIdx=1:length(pos)
              if(~isempty(pos{fIdx}))
                pos{fIdx}=obj.applyBaseToPosPointCloud(pos{fIdx}(:,[1 2 3]),fIdx);
              end
          end
          detectionsBase.setPosMatrix(pos,err);

      end
      
      function detectionsBase=applyInvBaseToDetection(obj,detections,name)
          if(isempty(detections))
              detectionsBase=[];
              return;
          end
          detectionsBase=detections.copy();
          try
              if(nargin>2)&&(~isempty(name))
                  detections.addprop(name);
              end
          catch
          end;
          pos=detections.getPosMatrix();
          err=detections.getErrMatrix();
          for fIdx=1:length(pos)
            if(~isempty(pos{fIdx}))
              pos{fIdx}=obj.applyInvBaseToPosPointCloud(pos{fIdx}(:,[1 2 3]),fIdx);
            end
          end
          detectionsBase.setPosMatrix(pos,err);

      end

      function pos=applyBaseToPos(obj,pos,fIdx)
        B=obj.getBase(fIdx);
        O=obj.getOrigAtFrame(fIdx);

        for pIdx=1:size(pos,1)
          pos(pIdx,:)=[(pos(pIdx,1)-O(1)) (pos(pIdx,2)-O(2)) (pos(pIdx,3)-O(3))];
          pos(pIdx,:)=pos(pIdx,:)*B;
        end
      end


      function tform= getAffineTransform(obj,f)
        pIdx=find(obj.frame==f,1);
        if(isempty(pIdx))
            if(f>max(obj.frame))   pIdx=length(obj.frame);  else   pIdx=1; end;
        end
        B=[obj.X(pIdx,:)' obj.Y(pIdx,:)' obj.Z(pIdx,:)'];
        
        
        
        
        Orig=obj.getOrigAtFrame(f);
        A(1:3,1:3)=double(B);
        A(4,[1 2 3])=-Orig*B;
        A(4,4)=1;
        tform=affine3d(A);
      end

      function pos=applyBaseToPosPointCloud(obj,pos,fIdx)
        PC=pointCloud(pos);
        tform=obj.getAffineTransform(fIdx);
        PCRef = pctransform(PC,tform);    
        pos=PCRef.Location;
      end

      function pos=applyInvBaseToPosPointCloud(obj,pos,fIdx)
        PC=pointCloud(pos);
        tform=obj.getAffineTransform(fIdx);
        try
            PCRef = pctransform(PC,invert(tform));
            pos=PCRef.Location;
        catch 'vision:pointcloud:rigidTransformOnly'
            tform.T
            det(tform.T(1:3,1:3))
            [X,Y,Z]=transformPointsInverse(tform,pos(:,1),pos(:,2),pos(:,3));
            pos(:,1)=X; pos(:,2)=Y; pos(:,3)=Z;
        end
      end

    end

    % methods (Static)
    % function newPoints=applyBasePoint(B,O,points)
    %   for pIdx=1:size(points,1)
    %     recentered=[(points(pIdx,1)-O(1)) (points(pIdx,2)-O(2) (points(pIdx,3)-O(3))];
    %       v=recentered*B;
    %       newPoints(pIdx,:)=v;
    %     end;
    %   end
end
