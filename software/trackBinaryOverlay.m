function trackXY=   trackBinaryOverlay(img,XLimit,YLimit,tracks,frameIdx,colorIndex,colormap,varargin)
ip = inputParser;
ip.CaseSensitive = false;
ip.KeepUnmatched = true;
  ip.addOptional('cumulative',false);
  ip.addOptional('dragonTail',[]);
  ip.addOptional('trackIDs',1:numel(tracks));
  ip.addOptional('insertTrackID',false);
  ip.parse(varargin{:});
  p=ip.Results;
%%
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
trackXY=img;
if(length(size(img))<3)
    trackXY=repmat(img,1,1,3);
end

tracksColors=colormap;

if(~isempty(tracks))
    sampling=50;
    for trIdx=1:length(tracks)

        t=tracks(trIdx);
        tIdx=find(t.f==(frameIdx),1);
        if(~isempty(tIdx))
 %%
            if(p.cumulative)
              tIdx=length(t.f);
            end

            displayStart=1;
            if(~isempty(p.dragonTail))
                displayStart=max(1,tIdx-p.dragonTail);
            end
            X=t.x(displayStart:tIdx); Y=t.y(displayStart:tIdx);% Z=t.z(1:tIdx);

            if(iscell(colorIndex))
                nFrame=numel(t.f);
                allRGB=tracksColors(colorIndex{trIdx}(displayStart:min(tIdx,nFrame-1)),:);
            else
                allRGB=tracksColors(colorIndex(trIdx),:);
            end
            % Detect segment outside the limit.
            X=X-XLimit(1);
            Y=Y-YLimit(1);
            XRatio=size(img,2)/(XLimit(2)-XLimit(1));
            YRatio=size(img,1)/(YLimit(2)-YLimit(1));
            X=X*XRatio;
            Y=Y*YRatio;
            allInIdx=(X>0)&(Y>0)&(X<=size(img,2))&(Y<=size(img,1));
            allX=X;
            allY=Y;
            
            trackletInImage=bwconncomp(allInIdx);
            for tlIdx=1:trackletInImage.NumObjects
                inIdx=trackletInImage.PixelIdxList{tlIdx};
                if(numel(inIdx)>1)
                    X=allX(inIdx);
                    Y=allY(inIdx);
                    xSeg=max(1,round(linspaceMult(X(1:(end-1)),X(2:(end)),sampling)));
                    ySeg=max(1,round(linspaceMult(Y(1:(end-1)),Y(2:(end)),sampling)));
                    if(iscell(colorIndex))
                        RGB=allRGB(inIdx(1:end-1),:);
                        R=max(1,round(linspaceMult(RGB(:,1),RGB(:,1),sampling)));
                        G=max(1,round(linspaceMult(RGB(:,2),RGB(:,2),sampling)));
                        B=max(1,round(linspaceMult(RGB(:,3),RGB(:,3),sampling)));
                        indx=sub2ind(size(trackXY),ySeg,xSeg,ones(size(xSeg)));
                        trackXY(indx)=R;
                        indx=sub2ind(size(trackXY),ySeg,xSeg,2*ones(size(xSeg)));
                        trackXY(indx)=G;
                        indx=sub2ind(size(trackXY),ySeg,xSeg,3*ones(size(xSeg)));
                        trackXY(indx)=B;
                    else
                        indx=sub2ind(size(trackXY),ySeg,xSeg,ones(size(xSeg))); 
                        trackXY(indx)=allRGB(1);
                        indx=sub2ind(size(trackXY),ySeg,xSeg,2*ones(size(xSeg)));
                        trackXY(indx)=allRGB(2);
                        indx=sub2ind(size(trackXY),ySeg,xSeg,3*ones(size(xSeg)));
                        trackXY(indx)=allRGB(3);
                    end
                end

            end
            if(p.insertTrackID)
                 trackXY=insertText(trackXY,[X(end) Y(end)],num2str(p.trackIDs(trIdx)),'FontSize',16);
             end
        end
    end

end

function res=linspaceMult(X,Y,N)
    res=zeros(length(X),N);
    for r=1:size(res,1)
        res(r,:)=linspace(X(r),Y(r),N);
    end

