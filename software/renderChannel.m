function RGBVol=renderChannel(ch1,ch2,type,varargin)
    switch(type)
        case 'grayRed'
          RGBVol=grayRedRender(ch1,ch2);
        case 'redGray'
          RGBVol=grayRedRender(ch2,ch1);
        case 'greenRed'
            RGBVol=greenRedRender(ch1,ch2);
        case 'stereo'
            RGBVol=uint8(255*sc(cat(3,ch1,ch2),'stereo'));
        case 'stereoInv'
            RGBVol=uint8(255*sc(cat(3,ch2,ch1),'stereo'));
        case 'testMix'
            img=sc(ch1,'gray').*sc(ch1,'winter')+sc(ch2,'gray').*sc(ch2,'hot');
            RGBVol=uint8(255*img);
        case 'grayRedMix'
            img=sc(ch1,'gray')+sc(ch2,'gray').*sc(ch2,'-autumn');
            RGBVol=uint8(255*img);
        case 'redGrayMix'
            img=sc(ch2,'gray')+sc(ch1,'gray').*sc(ch1,'-autumn');
            RGBVol=uint8(255*img);
        case 'cyanMagenta'
            cyan=linspaceNDim([0 0 0],[0 1 1],64)';
            mag=linspaceNDim([0 0 0],[1 0 1],64)';
            img=sc(ch1,'gray').*sc(ch1,[min(ch1(:)) max(ch1(:))],mag)+sc(ch2,'gray').*sc(ch2,[min(ch2(:)) max(ch2(:))],cyan);
            RGBVol=uint8(255*img);
        otherwise
            error('unknown channel renderer');
    end

function RGBVol=greenRedRender(greenCh,redCh)
    RGBVol=repmat(greenCh,1,1,3);
    %RGBThree(:,:,1)=max(rMaxXY,rmaxXYKin);
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
    RGBVol(:,:,1)=redCh;
    RGBVol(:,:,2)=greenCh;
    RGBVol(:,:,3)=0;


function RGBVol=grayRedRender(grayCh,redCh)
    RGBVol=repmat(grayCh,1,1,3);
    RGBVol(:,:,1)=max(grayCh,redCh);

