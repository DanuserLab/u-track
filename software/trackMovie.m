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
    outFilePaths{5,i} = [p.OutputDirectory filesep outFilename '.csv']; % saved tracksFinal in channel_N_tracking_result.mat to a csv format for Napari. - 2025-8
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

    % exports the tracking events (tracksFinal) into a format (csv) that is compatible with napari. - added 2025-8
    napariWriteTracks(tracksFinal, outFilePaths{5,i})

end

disp('Finished tracking!')



% export the tracking events into a format that is compatible with napari
% written by Kevin Dean, 2025-8
function napariWriteTracks(tracks, outputFile)
%EXPORTTRACKSTONAPARICSV Export u-track3D tracks to a Napari tracks CSV.
%  exportTracksToNapariCSV(tracks, outputFile) processes the given tracks structure 
%  (e.g. tracks = MovieData.processes_{1,3}.loadChannelOutput(1)) and writes a CSV 
%  file with columns: track_id, frame, z, y, x. The output is suitable for loading 
%  into Napari's Tracks layer. Gap-closed segments are interpolated and split/merge 
%  events are handled by splitting or ending track IDs appropriately.
%
% Parameters:
%    tracks      - Array of track structures from u-track3D, each with fields:
%                  .tracksFeatIndxCG, .tracksCoordAmpCG, .seqOfEvents.
%    outputFile  - Filename for the CSV output (e.g., 'tracks_napari.csv').
%
% Example:
%    tracks = MovieData.processes_{1,3}.loadChannelOutput(1);
%    exportTracksToNapariCSV(tracks, 'napari_tracks.csv');

    % Initialize an array to collect all track rows
    outputRows = [];  % will be an N x 5 matrix [track_id, frame, z, y, x]
    nextTrackID = 0;  % track_id counter (start from 0 for Napari)
    
    % Loop over each compound track in the input
    for ct = 1:length(tracks)
        trackStruct = tracks(ct);
        % Determine the global frame range spanned by this compound track.
        seq = trackStruct.seqOfEvents;
        if isempty(seq)
            continue;  % skip if no events (should not happen for valid track)
        end
        % Frames in seqOfEvents are 1-indexed absolute frames in the movie.
        firstFrame = min(seq(:,1));
        lastFrame  = max(seq(:,1));
        % Prepare shorthand variables for matrices
        coordMat = trackStruct.tracksCoordAmpCG;   % (rows = sub-tracks, cols = 8 * nFrames)
        indxMat  = trackStruct.tracksFeatIndxCG;   % (rows = sub-tracks, cols = nFrames)
        % Number of sub-tracks (rows in the matrices)
        nSubTracks = size(coordMat, 1);
        
        % Loop over each sub-track in this compound track
        for s = 1:nSubTracks
            % Skip empty sub-tracks (if any). A sub-track row is empty if it has no detections.
            if all(indxMat(s, :) == 0)
                continue;
            end
            % Determine start and end frame for this sub-track:
            relStartIdx = find(indxMat(s, :) > 0, 1, 'first');  % index in matrix where this track starts
            relEndIdx   = find(indxMat(s, :) > 0, 1, 'last');   % index where this track ends
            % Convert relative indices to absolute frame numbers:
            % relIdx = 1 corresponds to 'firstFrame'
            startFrameAbs = firstFrame + (relStartIdx - 1);
            endFrameAbs   = firstFrame + (relEndIdx   - 1);
            
            % Extract known coordinates for this track (at frames where a detection exists).
            % We will use these to interpolate gap positions.
            knownFrames = []; 
            knownPos = [];    % will be N_k x 3 (columns: x, y, z)
            for relIdx = relStartIdx:relEndIdx
                if indxMat(s, relIdx) ~= 0
                    % Absolute frame for this relative index:
                    frameAbs = firstFrame + (relIdx - 1);
                    % Extract (x, y, z) from coordMat. Each frame occupies 8 columns.
                    colStart = (relIdx - 1) * 8 + 1;
                    % Note: coordMat stores [x, y, z, amp, xStd, yStd, zStd, ampStd]
                    x = coordMat(s, colStart);
                    y = coordMat(s, colStart + 1);
                    z = coordMat(s, colStart + 2);
                    knownFrames(end+1, 1) = frameAbs;
                    knownPos(end+1, :)   = [x, y, z];
                end
            end
            
            % Interpolate to fill gaps between known detections.
            % We'll iterate through consecutive known points and linearly interpolate any gaps.
            trackRows = [];  % temporary storage for [track_id, frame, z, y, x] rows of this sub-track
            for k = 1:length(knownFrames)
                f_current = knownFrames(k);
                pos_current = knownPos(k, :);
                % Add the current known position to the track output (will interpolate after adding).
                % (We'll convert to [z, y, x] order when adding to trackRows.)
                trackRows(end+1, :) = [nextTrackID, (f_current - 1), pos_current(3), pos_current(2), pos_current(1)];
                
                if k < length(knownFrames)
                    % Look at the next known point to interpolate the interval
                    f_next = knownFrames(k+1);
                    pos_next = knownPos(k+1, :);
                    % Compute gap length (frames in between current and next)
                    gapLength = f_next - f_current - 1;
                    if gapLength > 0
                        % Interpolate linearly for each frame between f_current and f_next
                        for g = 1:gapLength
                            interpFrame = f_current + g;  % the frame to interpolate
                            t = g / (f_next - f_current); % fractional distance between current and next
                            % Linear interpolation of coordinates:
                            interpPos = pos_current + t * (pos_next - pos_current);
                            % Append the interpolated row (same track_id)
                            trackRows(end+1, :) = [nextTrackID, (interpFrame - 1), interpPos(3), interpPos(2), interpPos(1)];
                        end
                    end
                end
            end
            
            % Append this sub-track's rows to the overall output list
            outputRows = [outputRows; trackRows];
            % Increment track ID counter for the next track segment
            nextTrackID = nextTrackID + 1;
        end
    end
    
    % Sort output rows by track_id then by frame for neatness
    outputRows = sortrows(outputRows, [1 2]);
    % Write to CSV file
    writematrix(outputRows, outputFile);