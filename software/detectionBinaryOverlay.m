function res=detectionBinaryOverlay(img,XLimit,YLimit,detections,colorIndex,colormap,varargin)
ip = inputParser;
ip.CaseSensitive = false;
ip.KeepUnmatched = true;
ip.addOptional('cumulative',false);
ip.addOptional('renderingMethod','Circle');
ip.addOptional('printVectorFile',[]);
ip.addOptional('radius',2);
ip.addOptional('positionsLabel',[]);
ip.addOptional('opacity',1);
ip.addOptional('linewidth',2);
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
res=img;
if(size(img,3)==1)
    res=repmat(img,1,1,3);
end
XRatio=size(img,2)/(XLimit(2)-XLimit(1));
YRatio=size(img,1)/(YLimit(2)-YLimit(1));
detColors=colormap;


if(~isempty(detections))
    for fIdx=1:length(detections)
        d=detections(fIdx);
        if(iscell(colorIndex))
            colIdx=colorIndex{fIdx};
        else
            colIdx=colorIndex;
        end



        if(~isempty(d)&&~(isempty(d.xCoord)))
            X=d.xCoord(:,1); Y=d.yCoord(:,1);% Z=t.z(1:tIdx);
            
            X=X-XLimit(1);
            Y=Y-YLimit(1);

            X=X*XRatio;
            Y=Y*YRatio;
            if(iscell(p.radius))
                radius=ceil(p.radius{fIdx}*XRatio);
            else
                radius=ceil(p.radius*XRatio);
            end
            if(numel(radius)==1)
                    radius= repmat(radius,numel(X),1);
            end
            radius=reshape(radius,[numel(radius) 1]);
            
            inIdx=(X>=0)&(Y>=0)&(X<=size(img,2))&(Y<=size(img,1));
            X=X(inIdx);
            Y=Y(inIdx);
            

            cIndex=colIdx(inIdx);
            radius=radius(inIdx);

            
            


            switch p.renderingMethod
                case 'MidpointCircle'
                    drawingBoard=zeros(size(img,1),size(img,2));
                    for dIdx=1:length(X)
                        drawingBoard=MidpointCircle(drawingBoard,radius,Y(dIdx),X(dIdx),cIndex(dIdx));
                    end
                    uniqueCIdx=unique(cIndex);
                    for ucIdx=1:length(uniqueCIdx)
                        cIdx=uniqueCIdx(ucIdx);
                        [I,J] = find(drawingBoard==cIdx);
                        indx=sub2ind(size(res),I,J,1*ones(size(I)));
                        res(indx)=colormap(cIdx,1);
                        indx=sub2ind(size(res),I,J,2*ones(size(I)));
                        res(indx)=colormap(cIdx,2);
                        indx=sub2ind(size(res),I,J,3*ones(size(I)));
                        res(indx)=colormap(cIdx,3);
                    end
                case {'Rectangle', 'FilledRectangle', 'Circle','FilledCircle'}
                    res = insertShape(res,p.renderingMethod,[X Y radius],'linewidth',p.linewidth,'Color',colormap(cIndex,:),'Opacity',p.opacity);
%                 case {'inserShape','Rectangle', 'FilledRectangle', 'Line', 'Polygon', 'FilledPolygon', 'Circle','FilledCircle'}
                otherwise
                    disp('unknow detection renderer.')
            end

            if(~isempty(p.positionsLabel))
                positionsLabel=p.positionsLabel{fIdx};
                
                positionsLabel=positionsLabel(inIdx);
                for lIdx=1:numel(positionsLabel)
                    if(~isempty(positionsLabel{lIdx}))    
                        res=insertText(res,[X(lIdx) Y(lIdx)],positionsLabel{lIdx});
                    end
                end
            end

            if(~isempty(p.printVectorFile))
                F=figure();
                intensityMinMax=[min(img(:)) max(img(:))];
                axes('Position',[0 0 1 1]);
                hold on;
                imshow(img);
%                 xlim(XLimit);
%                 ylim(YLimit);
                A=radius.^2;
                scatter(X,Y,A,colormap(cIndex,:))
%                 uniqueCIdx=unique(cIndex);
%                 for ucIdx=1:length(uniqueCIdx)
%                     cIdx=uniqueCIdx(ucIdx)
%                     toPlot=(cIdx==cIndex);
%                     plot(X(toPlot),Y(toPlot),'o','MarkerSize',radius,'MarkerEdgeColor',colormap(cIdx,:)/256);
%                 end
                print(F,p.printVectorFile,'-dsvg');
            end
        end
    end

end

% Draw a circle in a matrix using the integer midpoint circle algorithm
% Does not miss or repeat pixels
% Created by : Peter Bone
% Created : 19th March 2007
function i = MidpointCircle(i, radius, xc, yc, value)


xc = int16(xc);
yc = int16(yc);
keeper=(xc<(size(i,1)-radius))&& ... 
       (yc<(size(i,2)-radius))&& ...
       (xc>radius)&& ...
       (yc>radius);
if(keeper)

x = int16(0);
y = int16(radius);
d = int16(1 - radius);

i(xc, yc+y) = value;
i(xc, yc-y) = value;
i(xc+y, yc) = value;
i(xc-y, yc) = value;

while ( x < y - 1 )
    x = x + 1;
    if ( d < 0 ) 
        d = d + x + x + 1;
    else 
        y = y - 1;
        a = x - y + 1;
        d = d + a + a;
    end
    i( x+xc,  y+yc) = value;
    i( y+xc,  x+yc) = value;
    i( y+xc, -x+yc) = value;
    i( x+xc, -y+yc) = value;
    i(-x+xc, -y+yc) = value;
    i(-y+xc, -x+yc) = value;
    i(-y+xc,  x+yc) = value;
    i(-x+xc,  y+yc) = value;
end
end
