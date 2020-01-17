function movieVar = movieInfrastructure(whatToDo,movieType,dir2saveMovie,...
    movieName,numFramesMovie,movieVar,iFrame)

switch whatToDo
    
    case 'initialize' %initialize movie

        switch movieType
            case 'mov'
                %eval(['MakeQTMovie start ', fullfile(dir2saveMovie,movieName) '.mov']);
%
% Copyright (C) 2020, Danuser Lab - UTSouthwestern 
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
                MakeQTMovie('start', [ fullfile(dir2saveMovie,movieName) '.mov']) % Edit by Tony so that file names with spaces can be used
            case {'mp4_unix','avi_unix'}
                frameDir = [dir2saveMovie filesep 'tmpFramesMovie'];
                [~,~] = mkdir(frameDir);
        end
        
    case 'addFrame' %add frame to movie
        
        switch movieType
            case 'mov'
                MakeQTMovie addfigure
                MakeQTMovie('framerate',10);
            case 'avi'
                movieVar(iFrame) = getframe(gcf);
            case {'mp4_unix','avi_unix'}
                frameDir = [dir2saveMovie filesep 'tmpFramesMovie'];
                ndigit = num2str(ceil(log10(numFramesMovie)));
                kformat = ['%.' ndigit 'd'];
                ext = '.png';
                %generate temporary EPS file
                print('-depsc2', '-loose', '-r300', [frameDir filesep 'frame.eps']);
                %save each frame as PNG in the 'tmpFramesMovie' directory
                options = '-quiet -colorspace rgb -density 150 -depth 8';
                src = [frameDir filesep 'frame.eps'];
                dest = [frameDir filesep 'frame_' num2str(iFrame, kformat) ext];
                cmd = ['convert ' options ' ' src ' ' dest];
                system(cmd);
        end
        
    case 'finalize' %finalize movie
        
        switch movieType
            case 'mov'
                MakeQTMovie finish
            case 'avi'   
                v = VideoWriter(fullfile(dir2saveMovie,movieName));
                v.FrameRate = 10;
                v.Quality = 95;
                open(v);
                writeVideo(v, movieVar);
                close(v);
            case {'mp4_unix','avi_unix'}
                frameDir = [dir2saveMovie filesep 'tmpFramesMovie'];
                ndigit = num2str(ceil(log10(numFramesMovie)));
                ext = '.png';
                cmd = ['ffmpeg -y -r 10 -i ' frameDir filesep 'frame_%0' ndigit 'd' ext ' -b 20M -r 10 ' fullfile(dir2saveMovie,movieName) '.' movieType(1:end-5)];
                system(cmd);
                rmdir(frameDir,'s')
        end
        
end