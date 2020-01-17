function analyzeMovieMotion(movieDataOrProcess,varargin)
% detectMovieComets detect comets in a movie using successive thresholds
%
% SYNOPSIS detectMovieComets(movieData,paramsIn)
%
% INPUT   
%   movieData - A MovieData object describing the movie to be processed
%
%   paramsIn - Structure with inputs for optional parameters. The
%   parameters should be stored as fields in the structure, with the field
%   names and possible values as described below:
%
%   Possible Parameter Structure Field Names:
%       ('FieldName' -> possible values)
%
%       ('OutputDirectory' -> character string)
%       A character string specifying the directory where to save 
%       the motion analysis results.
%
%       ('ChannelIndex' -> Positive integer scalar or vector)
%       The integer index of the channel(s) containing tracks to be analyzed.
%
%       ('probDim' -> Positive interger) 
%       Problem dimensionality. Default is 2.
% 
%       ('checkAsym' -> Boolean) 
%       1 to check for asymmetric tracks and to analyze their diffusion
%       after dimensionality reduction. Default: 0
% 
%       ('alphaValues' -> Row vector with 2 entrie) 
%       First entry is the alpha-value for MSS analysis (can take the values
%       0.2, 0.1, 0.05 and 0.01; see help of trackMSSAnalysis for most
%       up-to-date values allowed). Second entry is the alpha-value for 
%       asymmetry determination (can take  the values 0.2, 0.1, 0.05 and 0.01; 
%       see help of asymDeterm2D3D for most up-to-date values allowed).
%       Default: [0.05 0.1]. 
%
%       ('confRadMin' -> Positive integer scalar) 
%       1 to calculate the confinement radius of confined particles using 
%       the minimum positional standard deviation; OR
%       0 to calculate it using the mean positional standard deviation; OR
%       2 to approximate the confinement area by a rectangle and calculate 
%       the length of both edges.
%       Default:0
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

%
% Sebastien Besson, Feb 2012

%% Input
%Check input
ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('movieData', @isProcessOrMovieData);
ip.addOptional('paramsIn',[], @isstruct);
ip.parse(movieDataOrProcess,varargin{:});
paramsIn=ip.Results.paramsIn;

% Get the MotionAnalysisProcess and create it if it does not exist
[movieData, postProc] = getOwnerAndProcess(movieDataOrProcess,'MotionAnalysisProcess',true);

%Parse input, store in parameter structure
p = parseProcessParams(postProc,paramsIn);

%% --------------- Initialization ---------------%%

% Check detection process first
iTrackProc =movieData.getProcessIndex('TrackingProcess',1,1);

assert(~isempty(iTrackProc),['Tracking has not been run! '...
    'Please run tracking prior to post-processing!'])
trackProc = movieData.processes_{iTrackProc};

assert(all(trackProc.checkChannelOutput(p.ChannelIndex)),...
    ['Missing tracking output ! Please apply tracking before ' ...
    'running  post-processing!']);
    
% Set up the input directories (input images)
inFilePaths = cell(1,numel(movieData.channels_));
for i = p.ChannelIndex
    inFilePaths{1,i} = trackProc.outFilePaths_{1,i};
end
postProc.setInFilePaths(inFilePaths);
    
% Set up the output file
outFilePaths = cell(1,numel(movieData.channels_));
for i = p.ChannelIndex
    outFilePaths{1,i} = [p.OutputDirectory filesep 'channel_' num2str(i) '.mat'];
end
mkClrDir(p.OutputDirectory);
postProc.setOutFilePaths(outFilePaths);

%% --------------- Displacement field calculation ---------------%%% 

disp('Starting analyzing tracks motion...')

for i = p.ChannelIndex    
    tracks = trackProc.loadChannelOutput(i);
    try 
        if postProc.funParams_.driftCorrect ==1     
            [tracks] = scriptCorrectImageDrift(tracks,movieData);
        end
    catch 
        
    end
    diffAnalysisRes = trackDiffusionAnalysis1(tracks,...
        1,p.probDim,p.checkAsym,p.alphaValues,0,p.confRadMin); %#ok<NASGU>
    tracks = trackProc.loadChannelOutput(i);
    for j = 1:numel(tracks)
        tracks(j).classification = diffAnalysisRes(j).classification;
    end
    % save each projData in its own directory
    save(outFilePaths{1,i},'diffAnalysisRes', 'tracks')
end


disp('Finished analyzing motion!')
