function rgbImg=graphBinaryOverlay(img,XLimit,YLimit,positions,edges,colorIndex,colormap,varargin)
%% positions define in XLimit,YLimit Nx2
%% edges is a MxP matrix to connect position in line of P positions (N>2);
%% P. Roudot 2018.
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

ip = inputParser;
ip.CaseSensitive = false;
ip.KeepUnmatched = true;
ip.addOptional('linewidth',2);
ip.addOptional('positionsLabel',{});   % position associated label (none if empty)
ip.parse(varargin{:});
p=ip.Results;
%%

XRatio=size(img,2)/(XLimit(2)-XLimit(1));
YRatio=size(img,1)/(YLimit(2)-YLimit(1));

rgbImg=img;

if(~isempty(edges))
    polygons=cell(1,size(edges,1));
    if(~isempty(positions))
        X=positions(:,1); Y=positions(:,2);
        X=X-XLimit(1);
        Y=Y-YLimit(1);
        X=X*XRatio;
        Y=Y*YRatio;
        % inIdx=(X>0)&(Y>0)&(X<=size(img,2))&(Y<=size(img,1));
        % X=X(inIdx);
        % Y=Y(inIdx);
        positions=[X Y];
    end
    for eIdx=1:size(edges,1)
        % if(numel(p.radius)>1)
        %     radius=ceil(p.radius(pIdx)*XRatio);
        % else
        %     radius=ceil(p.radius*XRatio);
        % end
        pol=zeros(1,4);
        pol(1:2:end)=positions(edges(eIdx,:),1);
        pol(2:2:end)=positions(edges(eIdx,:),2);
        polygons{eIdx}=pol;
    end
    rgbImg = insertShape(img,'Line',polygons,'linewidth',p.linewidth,'Color',colormap(colorIndex,:));

    if(~isempty(p.positionsLabel))
        for lIdx=1:numel(p.positionsLabel)
            if(~isempty(p.positionsLabel{lIdx}))    
                rgbImg=insertText(rgbImg,positions(lIdx,:),p.positionsLabel{lIdx},'BoxColor',[255 200 0]);
            end
        end
    end
    
end

