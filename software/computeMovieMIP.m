function computeMovieMIP(movieDataOrProcess, varargin)
%computeMovieMIP MovieData-Process framework wrapped printMIP function
% Andrew R. Jamieson, Aug. 2017
% Also See printMIP                                    
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
    
ip = inputParser;
ip.addRequired('MD', @(x) isa(x,'MovieData') || isa(x,'Process') && isa(x.getOwner(),'MovieData'));
ip.addOptional('paramsIn',[], @isstruct);
ip.parse(movieDataOrProcess, varargin{:});
paramsIn = ip.Results.paramsIn;

% turn a specific warning off
warning('off', 'MATLAB:imagesci:tifftagsread:expectedAsciiDataFormat');

% Gather process and MD info/parameters
[MD, thisProcess, iProc] = getOwnerAndProcess(movieDataOrProcess, 'ComputeMIPProcess', true);

if(MD.zSize_==1)
    warning('This seems to be a 2D movie, No MIP produced.');
    return;
end

% % assert(MD == movieData1);
% MD = movieData1;
%Parse input, store in parameter structure
p = parseProcessParams(thisProcess, paramsIn);

if ~isempty(p.ProcessIndex)
    % verify input process type is appropriate 
    for c = p.ChannelIndex;
        assert(~isempty(imDir(MD.processes_{p.ProcessIndex}.outFilePaths_{c})),'no process images exist')
    end
end

% ============= Configure InputPaths. ================
inFilePaths = cell(1, numel(MD.channels_));
for j = p.ChannelIndex
    if isempty(p.ProcessIndex)
        inFilePaths{1,j} = MD.getChannelPaths{j};
    else
        inFilePaths{1,j} = MD.processes_{p.ProcessIndex}.outFilePaths_{j};
    end
end
thisProcess.setInFilePaths(inFilePaths);

% ================[OUTPUT]===========================
mkClrDir(p.OutputDirectory, false);
outFilePaths = cell(5, numel(MD.channels_));
for i = p.ChannelIndex;    
    outFilePaths{1,i} = [p.OutputDirectory filesep 'ch' num2str(i) filesep 'XY'];
    outFilePaths{2,i} = [p.OutputDirectory filesep 'ch' num2str(i) filesep 'ZY'];
    outFilePaths{3,i} = [p.OutputDirectory filesep 'ch' num2str(i) filesep 'ZX'];
    outFilePaths{4,i} = [p.OutputDirectory filesep 'ch' num2str(i) filesep 'three'];
    outFilePaths{5,i} = [p.OutputDirectory filesep 'ch' num2str(i)];
    mkClrDir(outFilePaths{1,i});
    mkClrDir(outFilePaths{2,i});
    mkClrDir(outFilePaths{3,i});
    mkClrDir(outFilePaths{4,i});
end
thisProcess.setOutFilePaths(outFilePaths);

ZXRatio = MD.pixelSizeZ_/MD.pixelSize_;
minIntensityNorm = [];
maxIntensityNorm = [];

% prep for parfor linearization 
pRunSet = {}; cell(MD.nFrames_ * length(p.ChannelIndex));
i = 1;
for chIdx = p.ChannelIndex
    for frameIdx = 1:MD.nFrames_
        pRunSet{i}.chIdx = chIdx; %struct('chIdx', chIdx)
        pRunSet{i}.frameIdx = frameIdx; %struct('chIdx', chIdx)
        i = i + 1;
    end
end

parfor pfi = 1:length(pRunSet)  
% for chIdx = p.ChannelIndex
    chIdx = pRunSet{pfi}.chIdx;
    frameIdx = pRunSet{pfi}.frameIdx;
    
    savePathXY = outFilePaths{1, chIdx};
    savePathZY = outFilePaths{2, chIdx};
    savePathZX = outFilePaths{3, chIdx};
    savePathThree = outFilePaths{4, chIdx};
    savePath = outFilePaths{5, chIdx};
    if ~isdir(savePath) || ~isdir(savePathXY) || ~isdir(savePathZY) || ~isdir(savePathZX) || ~isdir([savePath filesep 'Three'])
        mkdirRobust([savePath]);
        mkdirRobust([savePathXY]);
        mkdirRobust([savePathZY]);
        mkdirRobust([savePathZX]);
        mkdirRobust([savePathThree]);
    end

    XYFilesPattern = [savePathXY filesep 'XY_frame_nb%04d.png'];
    YZFilesPattern = [savePathZY filesep 'ZY_frame_nb%04d.png'];
    XZFilesPattern = [savePathZX filesep 'ZX_frame_nb%04d.png'];
    ThreeFilesPattern = [savePathThree filesep 'Three_frame_nb%04d.png'];

    if isempty(p.ProcessIndex)
        vol = MD.getChannel(chIdx).loadStack(1);
    else
        vol = MD.processes_{p.ProcessIndex}.loadChannelOutput(chIdx,1);
    end

    minIntensityNorm = min(vol(:));
    maxIntensityNorm = max(vol(:));

%     parfor frameIdx = 1:MD.nFrames_
        if isempty(p.ProcessIndex)
            vol = MD.getChannel(chIdx).loadStack(frameIdx);
        else
            vol = MD.processes_{p.ProcessIndex}.loadChannelOutput(chIdx,frameIdx);
        end  
        [maxXY, maxZY, maxZX, three] = computeMIPs(vol, ZXRatio, minIntensityNorm, maxIntensityNorm);
        
        % save the maximum intensity projections
        imwrite(maxXY, sprintfPath(XYFilesPattern, frameIdx), 'Compression', 'none');
        imwrite(maxZY, sprintfPath(YZFilesPattern, frameIdx), 'Compression', 'none');
        imwrite(maxZX, sprintfPath(XZFilesPattern, frameIdx), 'Compression', 'none');
        imwrite(three, sprintfPath(ThreeFilesPattern, frameIdx), 'Compression', 'none');
%     end
end

parfor chIdx = p.ChannelIndex
    % savePath = outFilePaths{1,chIdx};
    ThreeFilesPattern = [outFilePaths{4, chIdx} filesep 'Three_frame_nb%04d.png'];
    threeVideo = VideoWriter([outFilePaths{5, chIdx} filesep 'threeMontage.avi']);
    % myVideo.FrameRate = 4;  % Default 30
    % myVideo.Quality = 90;    % Default 75

    open(threeVideo);
    for frameIdx = 1:MD.nFrames_
        three=imread(sprintfPath(ThreeFilesPattern, frameIdx));
        writeVideo(threeVideo,three)
    end
    close(threeVideo);
end