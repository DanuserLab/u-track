function out=roundOddOrEven(x,oddOrEven,infOrZero)
%rounds to the next even or odd number
%
%SYNOPSIS out=roundOddOrEven(x,oddOrEven,infOrZero)
%
%INPUT    x: vector or matrix of values
%         optional:
%         oddOrEven: ['odd'],'even': if round to odd or even number [default]
%         infOrZero: 'inf',['close'],'zero': if round towards infinity, closest
%         number or zero [default]
%
%OUTPUT   out: rounded values
%
%created: 2/03, Jonas
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

%test input
switch nargin
    case 1
        oddOrEven='odd';
        infOrZero='close';
    case 2
        infOrZero='close';
        if ~strcmp(oddOrEven,'odd')&~strcmp(oddOrEven,'even')
            error('wrong input for oddOrEven');
            return
        end
    case 3
        if ~strcmp(oddOrEven,'odd')&~strcmp(oddOrEven,'even')
            error('wrong input for oddOrEven');
            return
        end
        if ~strcmp(infOrZero,'inf')&~strcmp(infOrZero,'zero')&~strcmp(infOrZero,'close')
            error('wrong input for infOrZero');
            return
        end
end

%evaluation switch
switch infOrZero
    case 'inf'
        switch oddOrEven
            case 'odd'
                %remember sign
                sig=sign(x);
                % sign(0) is 0. Round 0 toward 1
                sig(sig==0) = 1;
                %round
                x=ceil(sig.*x).*sig;
                %make odd
                out=x+(1-mod(x,2)).*sig;
            case 'even'
                %remember sign
                sig=sign(x);
                %round
                x=ceil(sig.*x).*sig;
                %make even
                out=x+(mod(x,2)).*sig;
        end
    case 'zero'
        switch oddOrEven
            case 'odd'
                %remember sign
                sig=sign(x);
                % sign(0) is 0. Round 0 toward -1
                sig(sig==0) = -1;
                %round
                x=fix(x);
                %make odd
                %first: check for zeros
                zeroIdx=find(x==0);
                if ~isempty(zeroIdx)
                    x(zeroIdx)=sig(zeroIdx);
                end
                out=x-(1-mod(x,2)).*sig;
            case 'even'
                %remember sign
                sig=sign(x);
                %round
                x=fix(x);
                %make even
                out=x-(mod(x,2)).*sig;
        end
    case 'close'
        switch oddOrEven
            case 'odd'
                %remember sign
                sig=sign(x);
                %remember modulo
                modulo=floor(mod(x.*sig,2));
                %round
                x=fix(x);
                %make odd
                out=x+(1-modulo).*sig;
            case 'even'
                %remember sign
                sig=sign(x);
                %remember modulo
                modulo=floor(mod(x.*sig,2));
                %round
                x=fix(x);
                %make even
                out=x+(modulo).*sig;
        end
end