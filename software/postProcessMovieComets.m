function postProcessMovieComets(movieData,varargin)
% Track features in a movie which has been processed by a detection method
%
% Sebastien Besson, Feb 2012
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

%% Input
%Check input
ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('movieData', @(x) isa(x,'MovieData'));
ip.addOptional('paramsIn',[], @isstruct);
ip.parse(movieData,varargin{:});
paramsIn=ip.Results.paramsIn;


%Get the indices of any previous tracking processes from this function
iProc = movieData.getProcessIndex('CometPostTrackingProcess',1,0);

%If the process doesn't exist, create it
if isempty(iProc)
    iProc = numel(movieData.processes_)+1;
    movieData.addProcess(CometPostTrackingProcess(movieData));
end
postProc = movieData.processes_{iProc};

%Parse input, store in parameter structure
p = parseProcessParams(postProc,paramsIn);

%% --------------- Initialization ---------------%%

errorMsg = @(metadata) ['Missing ' metadata '. Please in the movie''s '...
    metadata ' before  post-processing the tracks.'];
assert(~isempty(movieData.timeInterval_), errorMsg('time interval'));
assert(~isempty(movieData.pixelSize_), errorMsg('pixels size'));

% Check detection process first
iTrackProc =movieData.getProcessIndex('TrackingProcess',1,1);
assert(~isempty(iTrackProc),['Tracking has not been run! '...
    'Please run tracking prior to post-processing!'])
trackProc = movieData.processes_{iTrackProc};

assert(all(trackProc.checkChannelOutput(p.ChannelIndex)),...
    ['Missing tracking output ! Please apply tracking before ' ...
    'running  post-processing!'])

iDetProc =movieData.getProcessIndex('DetectionProcess',1,1);
assert(~isempty(iDetProc),'Please run detection first');
detProc=movieData.processes_{iDetProc};

assert(all(detProc.checkChannelOutput(p.ChannelIndex)),...
    ['Missing detection output ! Please apply detection before ' ...
    'running post-processing!']);

% Set up the input directories (input images)
inFilePaths = cell(2,numel(movieData.channels_));
for i = p.ChannelIndex
    inFilePaths{1,i} = trackProc.outFilePaths_{1,i};
    inFilePaths{2,i} = detProc.outFilePaths_{1,i};
end
postProc.setInFilePaths(inFilePaths);

% Set up the output file
outFilePaths = cell(2,numel(movieData.channels_));
for i = p.ChannelIndex
    outFilePaths{1,i} = [p.OutputDirectory filesep 'channel_' num2str(i) '.mat'];
    outFilePaths{2,i} = [p.OutputDirectory filesep 'channel_' num2str(i) '.txt'];
    outFilePaths{3,i} = [p.OutputDirectory filesep 'channel_' num2str(i) '_stats.txt'];
    outFilePaths{4,i} = [p.OutputDirectory filesep 'channel_' num2str(i) '_histograms'];
end
if isdir(p.OutputDirectory), rmdir(p.OutputDirectory, 's'); end
mkdir(p.OutputDirectory);
postProc.setOutFilePaths(outFilePaths);

%% --------------- Displacement field calculation ---------------%%%

disp('Starting post-processing...')

for i = p.ChannelIndex
    movieInfo = detProc.loadChannelOutput(i);
    tracksFinal = trackProc.loadChannelOutput(i);
    
    projData.secPerFrame = movieData.timeInterval_;
    projData.pixSizeNm = movieData.pixelSize_;
    
    % figure out which frames were used in detection
    detExists=find(arrayfun(@(x) ~isempty(x.xCoord),movieInfo));
    sF=min(detExists); eF=max(detExists);
    projData.detectionFrameRange=[sF eF];
    
    % Read tracking parameters
    gapCloseParam = trackProc.funParams_.gapCloseParam;
    costMatrices = trackProc.funParams_.costMatrices;
    projData.trackingParameters.maxGapLength=gapCloseParam.timeWindow;
    projData.trackingParameters.minTrackLen=gapCloseParam.minTrackLen;
    projData.trackingParameters.minSearchRadius=costMatrices(1,1).parameters.minSearchRadius;
    projData.trackingParameters.maxSearchRadius=costMatrices(1,1).parameters.maxSearchRadius;
    projData.trackingParameters.maxForwardAngle=costMatrices(1,2).parameters.maxFAngle;
    projData.trackingParameters.maxBackwardAngle=costMatrices(1,2).parameters.maxBAngle;
    projData.trackingParameters.backVelMultFactor=costMatrices(1,2).parameters.backVelMultFactor;
    projData.trackingParameters.fluctRadius=costMatrices(1,2).parameters.fluctRad;
    
    % Call main post-processing function
    [projData,M]=postProcessMTTracks(projData, tracksFinal, movieInfo,...
        [1 movieData.nFrames_],p.remBegEnd,...
        'fgapReclassScheme',p.fgapReclassScheme,...
        'bgapReclassScheme',p.bgapReclassScheme);
    
    % save each projData in its own directory
    save(outFilePaths{1,i},'projData')
    
    % write out speed/lifetime/displacement distributions into a text file
    dlmwrite(outFilePaths{2,i}, M,...
        'precision', 3,'delimiter', '\t','newline', 'pc');
    
    % Write stats results into a text file
    statsFile = outFilePaths{3, i};
    statsData= struct2cell(projData.stats);
    statsName = fieldnames(projData.stats);
    fid=fopen(statsFile,'w+');
    for j=1:numel(statsName)
        fprintf(fid,'%s\t%g\n',statsName{j},statsData{j});
    end
    fclose(fid);

    
    if p.makeHist==1
        plusTipMakeHistograms(M, outFilePaths{4, i});
        if movieData.isOmero() && movieData.canUpload()
            uploadHistograms(movieData, outFilePaths{4, i})
        end
        
    end
end

disp('Finished post-processing comet tracks!')

function uploadHistograms(movieData, directory)

disp('Uploading histograms to OMERO')

% Retrieve
files = dir(fullfile(directory, '*.eps'));
session =  movieData.getOmeroSession();
id = movieData.getOmeroId();
namespace = [getLCCBOmeroNamespace() '.tracking'];

for i = 1 : numel(files),
    [~, filename] = fileparts(files(i).name);
    ns = [namespace '.' filename];
    fas = getImageFileAnnotations(session, id, 'include', ns);
    
    if ~isempty(fas)
        % Read file of first found file annotation
        fa = fas(1);
        updateFileAnnotation(session, fa,...
            fullfile(directory, files(i).name));
        fprintf(1, 'Updating file annotation: %d\n', fa.getId().getValue());
    else
        fa = writeFileAnnotation(session,...
            fullfile(directory, files(i).name),...
            'description', 'Comet post-processing results', 'namespace', ns);
        linkAnnotation(session, fa, 'image', id);
        msg = 'Created file annotation %g and linked it to image %d\n';
        fprintf(1, msg, fa.getId().getValue(), id);
    end
end

