function trackMovie(processOrMovieData,varargin)
% Track features in a movie which has been processed by a detection method
%
% Sebastien Besson, 5/2011
% Updated Andrew R. Jamieson Mar 2017
%
% Copyright (C) 2025, Danuser Lab - UTSouthwestern 
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
ip.addRequired('movieData', @(x) isa(x,'Process') && isa(x.getOwner(),'MovieData') || isa(x,'MovieData'));
ip.addOptional('paramsIn',[], @isstruct);
ip.addParameter('ProcessIndex',[],@isnumeric);
ip.parse(processOrMovieData,varargin{:});
paramsIn=ip.Results.paramsIn;

% TrackingProcess default outputDirectory is owner.outputDirectory_
[movieData, trackProc] = getOwnerAndProcess(processOrMovieData,'TrackingProcess',true);

%Parse input, store in parameter structure
p = parseProcessParams(trackProc, paramsIn);

% precondition / error checking
if isa(trackProc, 'TrackingDynROIProcess')
%     If numel(buildDynROIProcId) > 1, popup window will show and let user to choose which BuildDynROIProcess. 
%     Later, added in the setting GUI for user to select a BuildDynROIProcess, so comment out.
%     buildDynROIProcId = movieData.getProcessIndex('BuildDynROIProcess');
    if isempty(p.processBuildDynROI)
        error("BuildDynROIProcess needs to be done and selected in setting before run TrackingDynROIProcess.")
    elseif ~ismember(1, p.processBuildDynROI.funParams_.ChannelIndex)
        error("Channel 1 in BuildDynROIProcess needs to be analyzed before run TrackingDynROIProcess.")
    end
end

%% --------------- Initialization ---------------%%

% Check detection process first
if isempty(p.DetProcessIndex)
    % show only relevant detection processes in the list selection dialog
    % box for TrackingDynROIProcess. edit 2021-01-06
    if isa(trackProc, 'TrackingDynROIProcess')
        p.DetProcessIndex = movieData.getProcessIndex('PointSourceDetectionProcess3DDynROI',1,1);
    else
        p.DetProcessIndex = movieData.getProcessIndex('DetectionProcess',1,1);
    end

    if isempty(p.DetProcessIndex)
        error(['Detection has not been run! '...
            'Please run detection prior to tracking!'])
    end
end
detProc = movieData.processes_{p.DetProcessIndex};

if ~detProc.checkChannelOutput(p.ChannelIndex)
    error(['Missing detection output ! Please apply detection before ' ...
        'running tracking!'])
end

% Set up the input directories (input images)
inFilePaths = cell(1,numel(movieData.channels_));
for i = p.ChannelIndex
    inFilePaths{1,i} = detProc.outFilePaths_{1,i};
end
trackProc.setInFilePaths(inFilePaths);
    
% Set up the output file
outFilePaths = cell(1,numel(movieData.channels_));
for i = p.ChannelIndex
    outFilename= ['Channel_' num2str(i) '_tracking_result'];
    outFilePaths{1,i} = [p.OutputDirectory filesep outFilename '.mat'];
    if p.saveResults.export
        outFilePaths{2,i} = [p.OutputDirectory filesep outFilename '_mat.mat'];
    end
    if movieData.is3D && p.saveResults.exportTrackabilityData
        outFilePaths{3,i} = [p.OutputDirectory filesep outFilename '_Trackability.mat'];
    end

    if movieData.is3D && (~isempty(p.processBuildDynROI))
        outFilePaths{4,i} = [p.OutputDirectory filesep outFilename '_DynROIRef.mat'];
    end
end
mkClrDir(p.OutputDirectory,0);
trackProc.setOutFilePaths(outFilePaths);





%% --------------- Displacement field calculation ---------------%%% 

disp('::::')
logMsg = @(chan) ['Tracking objects for channel ' ...
                  num2str(chan) ' under:'];
detDirs= detProc.outFilePaths_;

for i = p.ChannelIndex

    movieInfo = detProc.loadChannelOutput(i);

    %% --------------- Adding capability to select time range frames ---------------%%% 

    % use movieInfo to check start/end frames for tracking
    % then get rid of detected features in movieInfo outside this range
    if ~isempty(p.timeRange) && isequal(unique(size(p.timeRange)), [1 2])
        nFrames = length(movieInfo);
        if p.timeRange(1) <= p.timeRange(2) && p.timeRange(2) <= nFrames
            startFrame = p.timeRange(1);
        endFrame = p.timeRange(2);
    else
        startFrame = 1;
    endFrame = nFrames;
    end
    elseif ~isempty(p.timeRange)
        error('--TrackingProcess: timeRange should be [startFrame endFrame] or [] for all frames')
    end


    if ~isempty(p.timeRange)
        %% --------------- Adding capability to select time range frames ---------------%%% 
        % Initialize a new movieInfo structure array with the old movieInfo fields 
        % movieInfo fields may vary depending of the detection method 
        oldMovieInfo = movieInfo;
        clear movieInfo;
        movieFields = fieldnames(oldMovieInfo)';
        emptyFields = [movieFields; cell(size(movieFields))];
        movieInfo(nFrames,1) = struct(emptyFields{:});
        movieInfo(startFrame:endFrame,:) = oldMovieInfo(startFrame:endFrame,:);
        clear oldMovieInfo
        %% --------------- Adding capability to select time range frames ---------------%%% 
    end

    % If a dynROI specified, transform the coordinates accordingly
    if movieData.is3D && (~isempty(p.processBuildDynROI))
        if isa(p.processBuildDynROI, 'BuildDynROIProcess')
            dynROICell=p.processBuildDynROI.loadChannelOutput(1); % iChan is mandatory input for BuildDynROIProcess.loadChannelOutput, BuildDynROIProcess used in package.
        else
            dynROICell=p.processBuildDynROI.loadChannelOutput(); % iChan is not mandatory input for BuilDynROI.loadChannelOutput
        end
        dynROI=dynROICell{1};
        mappedDetections=dynROI.mapDetections(Detections(movieInfo));
        ref=dynROI.getDefaultRef();
        mappedDetectionsRef=ref.applyBase(mappedDetections);
        mappedDetectionsRef.addOffset(1000,1000,1000);
        movieInfo=mappedDetectionsRef.getStruct();
    end

    disp(logMsg(i))
    disp(detDirs{1,i});
    disp('Results will be saved under:')
    disp(outFilePaths{1,i});


    % Call function - return tracksFinal for reuse in the export
    % feature
    if movieData.is3D
        [tracksFinal,kalmanInfoLink,~,trackabilityData] = trackCloseGapsKalmanSparse(movieInfo, p.costMatrices, p.gapCloseParam,...
            p.kalmanFunctions, p.probDim, 0, p.verbose,'estimateTrackability',p.EstimateTrackability);
    else
        tracksFinal = trackCloseGapsKalmanSparse(movieInfo, p.costMatrices, p.gapCloseParam,...
        p.kalmanFunctions, p.probDim, 0, p.verbose);
    end
    



    if movieData.is3D && (~isempty(p.processBuildDynROI))
        %% replace the trajectory with original detections
        tracksFinalStripped=rmfield(tracksFinal,'tracksCoordAmpCG');
        tracksOrigDetect=TracksHandle(tracksFinalStripped,mappedDetections.getStruct()); 
        tracksFinal=tracksOrigDetect.getStruct()';

        %% Transfrom lab ref back in the DynROI Ref (useful for the GUI)
        if isa(p.processBuildDynROI, 'BuildDynROIProcess')
            dynROICell=p.processBuildDynROI.loadChannelOutput(1); % iChan is mandatory input for BuildDynROIProcess.loadChannelOutput, BuildDynROIProcess used in package.
        else
            dynROICell=p.processBuildDynROI.loadChannelOutput(); % iChan is not mandatory input for BuilDynROI.loadChannelOutput
        end
        dynROI=dynROICell{1};
        ref=dynROI.getDefaultRef();
        tracksDynROIRef=ref.applyBase(tracksOrigDetect);
        [BBmin,BBmax]=dynROI.getBoundingBox(ref);
        tracksDynROIRef.addOffset(-BBmin(1)+1,-BBmin(2)+1,-BBmin(3)+1);
        tracksFinalDynROIRef_oldFormat=tracksDynROIRef.getStruct();
        save(outFilePaths{4,i},'tracksDynROIRef','tracksFinalDynROIRef_oldFormat');

        if(p.EstimateTrackability)
            %% Project debugging information back into lab ref 
            samples=trackabilityData.samplesDetections.addOffset(-1000,-1000,-1000);
            samples=ref.applyInvBase(samples);
            trackabilityData.samplesDetections=samples;

            samples=trackabilityData.predExpectation.addOffset(-1000,-1000,-1000);
            samples=ref.applyInvBase(samples);
            trackabilityData.predExpectation=samples;

        end
    end

    if movieData.is3D && (p.EstimateTrackability)
        %% Track segment mapping
        segTrackability=cell(1,length(tracksFinal));
        for tIdx=1:length(tracksFinal)
            tr=TracksHandle(tracksFinal(tIdx));
            nonGap=~tr.gapMask();
            nonGap=nonGap(1:end-1);  % T-1 segments.
            linkIdx=find(nonGap);
            segTrackability{tIdx}=nan(size(nonGap));
            segTrackability{tIdx}(linkIdx)=arrayfun(@(pIdx) trackabilityData.trackabilityCost{tr.f(pIdx)+1}(tr.tracksFeatIndxCG(pIdx)),linkIdx);
        end
        trackabilityData.segTrackability=segTrackability;
    end
    
    save(outFilePaths{1,i},'tracksFinal');
    
    % Optional export
    if p.saveResults.export
        if ~p.gapCloseParam.mergeSplit
            [M.trackedFeatureInfo M.trackedFeatureIndx]=...
                convStruct2MatNoMS(tracksFinal);
        else
            [M.trackedFeatureInfo M.trackedFeatureIndx,M.trackStartRow,M.numSegments]=...
                convStruct2MatIgnoreMS(tracksFinal);
        end
        save(outFilePaths{2,i},'-struct','M');
        clear M;
    end

    if movieData.is3D && p.saveResults.exportTrackabilityData
        save(outFilePaths{3,i},'trackabilityData');
    end
end

disp('Finished tracking!')
