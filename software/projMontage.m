function three=projMontage(XY,XZ,YZ,ZUp,ZRight)
  if(nargin<4)
    ZUp=true;
  end
  if(nargin<5)
    ZRight=false;
  end

  if(ZRight)
    ZUp=true;
  end

  stripeSize=4;

  YZ_mount=YZ;
  YZ_row_position=size(XY,1)+stripeSize;
  if(~ZUp)
    YZ_row_position=0;
    if(ndims(XY)==3)
      YZ_mount=permute(YZ,[2 1 3]);
    else
      YZ_mount=permute(YZ,[2 1]);
    end
  end
  if(ndims(XY)==3)
        three=uint8(zeros(size(XY,1)+stripeSize+size(XZ,1),size(XY,2)+stripeSize+size(YZ_mount,2),3));
  else
        three=uint8(zeros(size(XY,1)+stripeSize+size(XZ,1),size(XY,2)+stripeSize+size(YZ_mount,2)));
  end

  three(1:size(XY,1),1:size(XY,2),:)=XY;
  three((size(XY,1)+stripeSize)+(1:size(XZ,1)),1:size(XZ,2),:)=XZ;
  three(YZ_row_position+(1:size(YZ_mount,1)),size(XY,2)+stripeSize+(1:size(YZ_mount,2)),:)=YZ_mount;

  if(ZRight)
    three=permute(three,[2 1 3]);
  end
  % threeTop =    [XY,        zeros(size(XY,1),stripeSize,size(XY,3)),                        zeros(size(XY,1),size(YZ,2),size(XY,3))];
  % threeBottom = [XZ,        0*ones(size(XZ,1),stripeSize,size(XY,3)),                       YZ];
  % three =       [threeTop;  ones(stripeSize,size(XY,2)+size(YZ,2)+stripeSize,size(XY,3));   threeBottom];
%
% Copyright (C) 2018, Danuser Lab - UTSouthwestern 
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
