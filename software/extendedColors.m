function color = extendedColors(colorSwitch,drawTable,betaTable)
%extendedColors gives a more extensive set of possible colors than the normal color switch
%
%SYNOPSIS color = extendedColors(colorSwitch, drawTable, betaTable)
%
%INPUT    colorSwitch: can be any char from list 'rgbymcadefhijlnpqstuvzkw' or a
%                      number from 1 to 24 (23 without white), or a vector,
%                      where the first entry designs the color and the second
%                      entry the brighness (>0 brighter, <0 darker)
%         drawTable  : if == 1, extendedColors draws the full color table
%         betaTable  : if the table is drawn, this gives the amount of
%                      increased/decreased brightness
%
%OUTPUT   color      : RGB-vector according to the colorswitch
%
%REMARKS Thanks to the creators of arrow3 for the typing!
%
%      r  red       1               
%      g  green     2
%      b  blue      3
%      y  yellow    4            
%      m  magenta   5           
%      c  cyan      6        
%      a  apple green 7  
%      d  dark gray 8   
%      e  evergreen 9    
%      f  fuchsia   10                  
%      h  honey     11     
%      i  indigo    12                  
%      j  jade      13               
%      l  lilac     14            
%      n  nutbrown  15         
%      p  pink      16                  
%      q  kumquat   17             
%      s  sky blue  18                    
%      t  tan       19                          
%      u  umber     20                            
%      v  violet    21                           
%      z  zinc      22
%      k  black     23
%      w  white     -1
%
% c: 04-03 jonas
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

%init valid colors
valCol = 'rgbymcadefhijlnpqstuvzkw';

%----------parse input
if nargin == 0
    error('not enough input arguments!');
end

%test colorSwitch
if nargin == 1 && ~isempty(colorSwitch)
    if length(colorSwitch)>2
        % recursively call extendedColors
        nColors = length(colorSwitch);
        color = zeros(nColors,3);
        for i=1:nColors
            color(i,:) = extendedColors(colorSwitch(i));
        end
        return
    end
    
    switch ischar(colorSwitch) + 2*isnumeric(colorSwitch)
        
        case 1 %char
            
            if ~strfind(valCol,colorSwitch) | length(colorSwitch)>1
                error('unsupported colorString!');
            else
                colorNum = strfind(valCol,colorSwitch);
                beta = 0;
            end
            
        case 2 %num
            
            if length(colorSwitch) == 2
                beta = colorSwitch(2);
            else
                beta = 0;
            end
            
            %take no chances with col num - make sure that there is no 0
            if colorSwitch(1)>0
            colorNum = round(abs(rem(colorSwitch(1)-1,23)))+1;
            else
                colorNum = -1;
            end
            
    end %colorSwitch switch

elseif nargin > 1 & isempty(colorSwitch)
    colorNum = 1;
    beta = 0;
    
else
    error('insufficient input arguments');    
    
end %colorSwitch if

if nargin > 3
    error ('too many input arguments!');
end


%assign color
cn = LocalColorTable(0);
if colorNum > 0
color = cn(colorNum,:);
else
    color = cn(24,:);
end

%brighten
if beta
    color = LocalBrighten(color,beta);
end

%draw table
if nargin > 1 & ~isempty(drawTable) & drawTable ~=0
    if nargin == 2 | isempty(betaTable)
        betaTable = 0;
    end
    LocalColorTable(1,betaTable);
end
    

   
   

%==================================================================
function [cn]=LocalColorTable(n,beta)
vc='rgbymcadefhijlnpqstuvzkw';       % valid color codes
%   r               g               b               y
cn=[1.00,0.00,0.00; 0.00,1.00,0.00; 0.00,0.00,1.00; 1.00,1.00,0.00;
%   m               c               a               d
    1.00,0.00,1.00; 0.00,1.00,1.00; 0.00,0.70,0.00; 0.40,0.40,0.40;
%   e               f               h               i
    0.00,0.40,0.00; 0.90,0.00,0.40; 1.00,0.80,0.00; 0.00,0.00,0.70;
%   j               l               n               p
    0.20,0.80,0.50; 0.80,0.40,0.80; 0.50,0.20,0.00; 1.00,0.40,0.60;
%   q               s               t               u
    1.00,0.40,0.00; 0.00,0.80,1.00; 0.80,0.40,0.00; 0.70,0.00,0.00;
%   v               z               k               w
    0.60,0.00,1.00; 0.60,0.60,0.60; 0.00,0.00,0.00; 1.00,1.00,1.00;];
if n                                 % plot color table
    figure;
   name={'red','green','blue','yellow','magenta','cyan',...
      'apple green',...
      'dark gray','evergreen','fuchsia','honey',...
      'indigo','jade','lilac','nutbrown',...
      'pink','kumquat','sky blue','tan',...
      'umber','violet','zinc','black','white'};
   c=['yhtn';'gjae';'csbi';'plmv';'frqu';'wzdk'];
   clf, set(gcf,'DefaultAxesXTick',[],'DefaultAxesYTick',[],...
      'DefaultAxesXTickLabel',[],'DefaultAxesYTickLabel',[],...
      'DefaultAxesXLim',[0,0.75],'DefaultAxesYLim',[0,0.75],...
      'DefaultRectangleEdgeColor','none');
   for i=1:24
      subplot(4,6,i); j=find(vc==c(i)); title(name{j});
      dark=LocalBrighten(cn(j,:),-beta);
      light=LocalBrighten(cn(j,:),beta);
      rectangle('Position',[0,0.00,0.75,0.25],'FaceColor',dark);
      rectangle('Position',[0,0.25,0.75,0.25],'FaceColor',cn(j,:));
      rectangle('Position',[0,0.50,0.75,0.25],'FaceColor',light);
      rectangle('Position',[0,0.00,0.75,0.75],'EdgeColor','k');
      if rem(i,6)==1
         set(gca,'YTickLabel',{'dark','normal','light'},...
            'YTick',[0.125,0.375,0.625]);
         if i==19
            text(0,-0.25,['{\bf\it extendedColors}  Named Color Table  ',...
                  '( \beta = ',num2str(beta),' )']);
         end
      end
   end
end

%==========================================================================
% Brighten
function C=LocalBrighten(c,beta)
C=c.^((1-min(1-sqrt(eps),abs(beta)))^sign(beta));