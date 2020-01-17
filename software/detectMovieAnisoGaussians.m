function detectMovieAnisoGaussians(movieDataOrProcess,varargin)
% detectMovieAnisoGaussians detect objects by fitting anisotropic Gaussians
%
% SYNOPSIS detectMovieSubResFeatures(movieData,paramsIn)
%
% INPUT   
%   movieData - A MovieData object describing the movie to be processed
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

% Sebastien Besson, May 2012

%% ----------- Input ----------- %%

%Check input
ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('movieData', @isProcessOrMovieData);
ip.addOptional('paramsIn',[], @isstruct);
ip.parse(movieDataOrProcess,varargin{:});
paramsIn=ip.Results.paramsIn;

% Get MovieData object and Detection Process
[movieData, detProc] = getOwnerAndProcess(movieDataOrProcess,'AnisoGaussianDetectionProcess',true);

%Parse input, store in parameter structure
p = parseProcessParams(detProc,paramsIn);

%% --------------- Initialization ---------------%%

nChan=numel(movieData.channels_);
% Set up the input directories
inFilePaths = cell(1,nChan);
for i = p.ChannelIndex
    inFilePaths{1,i} = movieData.getChannelPaths{i};
end
detProc.setInFilePaths(inFilePaths);
    
% Set up the output directories
outFilePaths = cell(1,nChan);
saveResults(nChan,1)=struct();
for i = p.ChannelIndex;    
    saveResults(i).dir = p.OutputDirectory ;
    saveResults(i).filename = ['Channel_' num2str(i) '_detection_result.mat'];
    %Create string for current directory
    outFilePaths{1,i} = [saveResults(i).dir filesep saveResults(i).filename ];
end
mkClrDir(p.OutputDirectory);
detProc.setOutFilePaths(outFilePaths);

%% --------------- Sub-resolution object detection ---------------%%% 
disp('Starting detecting anistropic Gaussians...')

roiMask = movieData.getROIMask;

for i = p.ChannelIndex
    disp(['Please wait, detecting objects for channel ' num2str(i)])
    disp(inFilePaths{1,i});
    disp('Results will be saved under:')
    disp(outFilePaths{1,i});
    
    movieInfo(movieData.nFrames_,1)=struct('xCoord',[],'yCoord',[],...
        'amp',[],'sigmaX',[],'sigmaY',[],'theta',[],'bkg',[]);
    progressText(0,'Detecting anisotropic Gaussians');
    for j=1:movieData.nFrames_
        I=double(movieData.channels_(i).loadImage(j));
        if ~isempty(p.MaskProcessIndex)
            maskProc = movieData.processes_{p.MaskProcessIndex};
            mask = maskProc.loadChannelOutput(p.MaskChannelIndex,j);
        else
            mask = true(movieData.imSize_);
        end
        % Call stand-alone subresolution detection function
        movieInfo(j) = cometDetection(I, mask & roiMask(:,:,j), p.psfSigma,...
            'mode',p.mode,'kSigma',p.kSigma,'alpha',p.alpha,'minDist',p.minDist);
        progressText(j/movieData.nFrames_,'Detecting anisotropic Gaussians');
    end
    save(outFilePaths{1,i} ,'movieInfo');
end

disp('Finished detecting objects...')

