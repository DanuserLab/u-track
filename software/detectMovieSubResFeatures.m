function detectMovieSubResFeatures(movieDataOrProcess,varargin)
% detectMovieSubResFeatures detect sub-resolution objects in a movie
%
% detectMovieSubResFeatures 
%
% SYNOPSIS detectMovieSubResFeatures(movieData,paramsIn)
%
% INPUT   
%   movieDataOrProcess - A MovieData object describing the movie to be processed
%                      - or a Process object containing info for processing (including the movie)
%
%   paramsIn - Structure with inputs for optional parameters. The
%   parameters should be stored as fields in the structure, with the field
%   names and possible values as described below
%
% OUTPUT   
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

% Sebastien Besson, Oct 2011

%% ----------- Input ----------- %%

%Check input
ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('movieDataOrProcess', @isProcessOrMovieData);
ip.addOptional('paramsIn',[], @isstruct);
ip.parse(movieDataOrProcess,varargin{:});
paramsIn=ip.Results.paramsIn;


% Get MovieData object and Process
[movieData, subResDetProc] = getOwnerAndProcess(movieDataOrProcess,'SubResolutionProcess',true);

%Parse input, store in parameter structure
p = parseProcessParams(subResDetProc,paramsIn);

%% --------------- Initialization ---------------%%

nChan=numel(movieData.channels_);
%Set up the input directories
inFilePaths = cell(1,nChan);
for i = p.ChannelIndex
    inFilePaths{1,i} = movieData.getChannelPaths{i};
end
subResDetProc.setInFilePaths(inFilePaths);
    
% Set up the output directories
outFilePaths = cell(1,nChan);
saveResults(nChan,1)=struct();
dName = 'detections_for_channel_';
for i = p.ChannelIndex;    
    currDir = [p.OutputDirectory filesep dName num2str(i)];
    saveResults(i).dir = currDir;
    saveResults(i).filename = ['Channel_' num2str(i) '_detection_result.mat'];
    %Create string for current directory
    outFilePaths{1,i} = [saveResults(i).dir filesep saveResults(i).filename ];
    subResDetProc.setOutFilePaths(outFilePaths{1,i},i);
    mkClrDir(currDir);
end

% Get ROI mask if any.
roiMask = movieData.getROIMask;
p.detectionParam.roiMask = roiMask;

%% --------------- Sub-resolution object detection ---------------%%% 
disp('Starting detecting diffraction-limited objects...')

for i = p.ChannelIndex
    disp(['Please wait, detecting objects for channel ' num2str(i)])
    disp(inFilePaths{1,i});
    disp('Results will be saved under:')
    disp(outFilePaths{1,i});
    
    % Retrieve information about the images
    if movieData.isOmero() || movieData.isBF()
        movieParam.channel = movieData.channels_(i);
    else
        [~, base, digits4Enum, ~] = getFilenameBody(movieData.getImageFileNames{i}{1});
        digits4Enum = length(digits4Enum);
        
        movieParam.imageDir = [inFilePaths{1,i} filesep];
        movieParam.filenameBase = base;
        movieParam.digits4Enum = digits4Enum;
    end    
    movieParam.firstImageNum=p.firstImageNum;
    movieParam.lastImageNum=p.lastImageNum;

    % Call stand-alone subresolution detection function
    movieInfo = detectSubResFeatures2D_StandAlone(movieParam, p.detectionParam, saveResults(i));
end

disp('Finished detecting diffraction-limited objects...')

