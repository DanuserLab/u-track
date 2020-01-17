function overlayTracksMovieNew(tracksFinal,startend,dragtailLength,...
    saveMovie,movieName,filterSigma,classifyGaps,highlightES,showRaw,...
    imageRange,onlyTracks,classifyLft,diffAnalysisRes,intensityScale,...
    colorTracks,firstImageFile,dir2saveMovie,minLength,plotFullScreen,...
    movieType)
%OVERLAYTRACKSMOVIENEW overlays tracks obtained via trackCloseGapsKalman on movies with variable color-coding schemes
%
%SYNPOSIS overlayTracksMovieNew(tracksFinal,startend,dragtailLength,...
%    saveMovie,movieName,filterSigma,classifyGaps,highlightES,showRaw,...
%    imageRange,onlyTracks,classifyLft,diffAnalysisRes,intensityScale,...
%    colorTracks,firstImageFile,dir2saveMovieminLength,plotFullScreen,...
%    movieType)
%
%INPUT  tracksFinal   : Output of trackCloseGapsKalman.
%       startend      : Row vector indicating first and last frame to
%                       include in movie. Format: [startframe endframe].
%                       Optional. Default: [(first frame with tracks) (last frame with tracks)]
%       dragtailLength: Length of drag tail (in frames).
%                       Optional. Default: 10 frames.
%                       ** If dragtailLength = 0, then no dragtail.
%                       ** To show tracks from their beginning to their end,
%                       set dragtailLength to any value longer than the
%                       movie.
%                       ** To show tracks statically while features dance
%                       on them, use -1.
%                       ** To show tracks from their beginning to their
%                       end, and to retain tracks even after the particle
%                       disappears, use -2.
%       saveMovie     : 1 to save movie (as Quicktime), 0 otherwise.
%                       Optional. Default: 0.
%       movieName     : filename for saving movie.
%                       Optional. Default: TrackMovie (if saveMovie = 1).
%       filterSigma   : 0 to overlay on raw image, PSF sigma to overlay on
%                       image filtered with given filterSigma.
%                       Optional. Default: 0.
%       classifyGaps  : 1 to classify gaps as "good" and "bad", depending
%                       on their length relative to the legnths of the
%                       segments they connect, 0 otherwise.
%                       Optional. Default: 1.
%       highlightES   : 1 to highlight track ends and starts, 0 otherwise.
%                       Optional. Default: 1.
%       showRaw       : 1 to add raw movie to the left of the movie with
%                       tracks overlaid, 2 to add raw movie at the top of
%                       the movie with tracks overlaid, 0 otherwise.
%                       Optional. Default: 0.
%       imageRange    : Image region to make movie out of, in the form:
%                       [min pixel X, max pixel X; min pixel Y, max pixel Y].
%                       Optional. Default: Whole image.
%       onlyTracks    : 1 to show only tracks without any symbols showing
%                       detections, closed gaps, merges and splits; 0 to
%                       show symbols on top of tracks.
%                       Optional. Default: 0.
%       classifyLft   : 1 to classify objects based on (1) whether they
%                       exist throughout the whole movie, (2) whether they
%                       appear OR disappear, and (3) whether they appear
%                       AND disappear; 0 otherwise.
%                       Optional. Default: 1.
%       diffAnalysisRes:Diffusion analysis results (either output of
%                       trackDiffusionAnalysis1 or trackTransientDiffusionAnalysis2).
%                       Needed if tracks/track segments are to be
%                       colored based on their diffusion classification.
%                       With this option, classifyGaps, highlightES and
%                       classifyLft are force-set to zero, regardless of input.
%                       Optional. Default: None.
%       intensityScale: 0 to autoscale every image in the movie, 1
%                       to have a fixed scale using intensity mean and std,
%                       2 to have a fixed scale using minimum and maximum
%                       intensities.
%                       Optional. Default: 1.
%       colorTracks   : 1 to color tracks by rotating through 7 different
%                       colors, 0 otherwise. With this option,
%                       classifyGaps, highlightES and classifyLft are
%                       force-set to zero, regardless of input.
%                       Option ignored if diffAnalysisRes is supplied.
%                       Optional. Default: 0.
%       firstImageFile: Name of the first image file in the folder of
%                       images that should be overlaid. The file has to be
%                       the first image that has been analyzed even if not
%                       plotted. If file is not specified [], user will be
%                       prompted to select the first image.
%                       Optional. Default: [].
%       dir2saveMovie:  Directory where to save output movie.
%                       If not input, movie will be saved in directory where
%                       images are located.
%                       Optional. Default: [].
%       minLength     : Minimum length of tracks to be ploted.
%                       Optional. Default: 1.
%       plotFullScreen: 1 the figure will be sized to cover the whole
%                       screen. In this way the movie will be of highest
%                       possible quality. default is 0.
%       movieType     : 'mov' to make a Quicktime movie using MakeQTMovie,
%                       'avi' to make AVI movie using Matlab's movie2avi,
%                       'mp4_unix', 'avi_unix' to make an MP4 or AVI movie
%                       using ImageMagick and ffmpeg. These options work
%                       only under linux or mac.
%                       Optional. Default: 'mov'.
%
%OUTPUT the movie.
%
%REMARKS Color-coding:
%        ** Without diffusion classification, all tracks have a neutral
%        color, while objects are color coded in the following way:
%               * Detected object just after appearance: Green circle.
%               * Detected object just before disappearance: Yellow
%                 circle.
%               * Detected object in middle of trajectory that spans
%                 whole movie: White circle.
%               * Detected object in middle of trajectory that appears OR
%                 disappears within movie: Magenta circle.
%               * Detected object in middle of trajectory that appears AND
%                 disappears within movie: Red circle.
%               * Gap that is short than both segments it connects: Cyan
%                 star.
%               * Gap that is longer than at least one ofthe segments it
%                 connects: Blue star.
%               * Object before and after splitting: Green diamond.
%               * OBject before and after merging: Yellow diamond.
%           When classifyGaps = 0, all gaps are cyan.
%           When highlighES = 0, no green and yellow circles.
%           When classifyLft = 0, all objets in middle of trajectory are white.
%
%       ** With diffusion classification, all objects and gaps have neutral
%       color (merges and splits are diamonds), while tracks and track
%       segments are color-coded in the following way:
%               * Type 1: Linear + 1D confined diffusion: Orange.
%               * Type 2: Linear + 1D normal diffusion: Red.
%               * Type 3: Linear + 1D super diffusion: Green.
%               * Type 4: Linear + too short for 1D classification: Yellow.
%               * Type 5: Random/Unclassified + 2D confined diffusion: Blue.
%               * Type 6: Random/Unclassified + 2D normal diffusion: Cyan.
%               * Type 7: Random/Unclassified + 2D super diffusion: Magenta.
%               * Type 8: Random + too short for 2D classification: Purple.
%               * Type 0: Too short for any analysis: Light pink.
%
%Khuloud Jaqaman, August 2007
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

%% input - basic

%check whether correct number of input arguments was used
if nargin < 1
    disp('--overlayTracksMovieNew: Incorrect number of input arguments!');
    return
end

%ask user for images
if nargin < 16 || isempty(firstImageFile)
    [fName,dirName] = uigetfile('*.tif','specify first image in the stack - specify very first image, even if not to be plotted');
else
    if iscell(firstImageFile)
        [fpath,fname,fno,fext]=getFilenameBody(firstImageFile{1});
        dirName=[fpath,filesep];
        fName=[fname,fno,fext];
    elseif ischar(firstImageFile)
        [fpath,fname,fno,fext]=getFilenameBody(firstImageFile);
        dirName=[fpath,filesep];
        fName=[fname,fno,fext];
    end
end

%if input is valid ...
if(isa(fName,'char') && isa(dirName,'char'))
    
    %get all file names in stack
    outFileList = getFileStackNames([dirName,fName]);
    numFiles = length(outFileList);
    
    %determine     diffAnalysisRes = diffAnalysisRes(indx);which frames the files correspond to, and generate the inverse map
    %indicate missing frames with a zero
    frame2fileMap = zeros(numFiles,1);
    for iFile = 1 : numFiles
        [~,~,frameNumStr] = getFilenameBody(outFileList{iFile});
        frameNum = str2double(frameNumStr);
        frame2fileMap(frameNum) = iFile;
    end
    
    %assign as number of frames the last frame number observed
    numFrames = frameNum;
    
    %read first image to get image size
    currentImage = imread(outFileList{1});
    [isx,isy] = size(currentImage);
    
else %else, exit
    
    disp('--overlayTracksMovieNew: Bad file selection');
    return
    
end

%keep only tracks with minimum requested length
if nargin < 18 || isempty(minLength)
    minLength = 1;
end
if minLength > 1
    criteria.lifeTime.min = minLength;
    indx = chooseTracks(tracksFinal,criteria);
    tracksFinal = tracksFinal(indx,:);
end

%get first and last frames where there are tracks
allEvents = vertcat(tracksFinal.seqOfEvents);
tracksFirstFrame = min(allEvents(:,1));
tracksLastFrame = max(allEvents(:,1));

%check startend and assign default if necessary
if nargin < 2 || isempty(startend)
    startend = [tracksFirstFrame tracksLastFrame];
else
    tracksFirstFrame = min(tracksFirstFrame,startend(1));
    tracksLastFrame = max(tracksLastFrame,startend(2));
end

%keep only the frames of interest
outFileList = outFileList(frame2fileMap(startend(1)):frame2fileMap(startend(2)));
frame2fileMap = frame2fileMap(startend(1):startend(2));
indxNotZero = find(frame2fileMap~=0);
frame2fileMap(indxNotZero) = frame2fileMap(indxNotZero) - frame2fileMap(indxNotZero(1)) + 1;

%get number of frames in movie to be made
numFramesMovie = diff(startend) + 1;

%check whether an area of interest was input
if nargin < 10 || isempty(imageRange)
    imageRange = [1 isx; 1 isy];
end

%% input - additional parameters

%check dragtailLength and assign default if not necessary
if nargin < 3 || isempty(dragtailLength)
    dragtailLength = 10;
end

%check whether to save movie
if nargin < 4 || isempty(saveMovie)
    saveMovie = 0;
end

%check name for saving movie
if saveMovie && (nargin < 5 || isempty(movieName))
    movieName = 'trackMovie.mov';
end

%check whether to use filtered images
if nargin < 6 || isempty(filterSigma)
    filterSigma = 0;
end

%check whether to color-code gaps
if nargin < 7 || isempty(classifyGaps)
    classifyGaps = 1;
end

%check whether to highligh track starts and ends
if nargin < 8 || isempty(highlightES)
    highlightES = 1;
end

%check whether to put raw movie adjacent to movie with tracks overlaid
if nargin < 9 || isempty(showRaw)
    showRaw = 0;
end

%check whether to plot tracks only or also symbols
if nargin < 11 || isempty(onlyTracks)
    onlyTracks = 0;
end

%check whether to classify lifetime
if nargin < 12 || isempty(classifyLft)
    classifyLft = 1;
end

%check whether to color-code tracks based on diffusion classification
%check whether diffusion classification is for overall tracks or transient
if nargin < 13 || isempty(diffAnalysisRes)
    diffAnalysisRes = [];
    transDiffClass = 0;
else
    classifyGaps = 0;
    highlightES = 0;
    classifyLft = 0;
    if isfield(diffAnalysisRes,'segmentClass')
        transDiffClass = 1;
    else
        transDiffClass = 0;
    end
    if minLength > 1
        diffAnalysisRes = diffAnalysisRes(indx);
    end
end

%check how to scale image intensity
if nargin < 14 || isempty(intensityScale)
    intensityScale = 1;
end

%check whether to color individual tracks
if nargin < 15 || isempty(colorTracks)
    colorTracks = 0;
else
    if ~isempty(diffAnalysisRes)
        colorTracks = 0;
    end
    if colorTracks == 1
        classifyGaps = 0;
        highlightES = 0;
        classifyLft = 0;
    end
end

%check where to save resulting movie
if saveMovie && (nargin < 17 || isempty(dir2saveMovie))
    dir2saveMovie = dirName;
end

%check whether to use full screen for plotting
if nargin < 19 || isempty(plotFullScreen)
    plotFullScreen = 0;
end

%decide on movie type
if nargin < 20 || isempty(movieType)
    movieType = 'mov';
end

%define colors to loop through
colorLoop = [1 0.7 0.7; 1 0 0; 0 1 0; 0 0 1; 1 1 0; 1 0 1; 0 1 1]; %colors: 'light pink',r,g,b,y,m,c

%% store track positions, get track status and point status

%get number of tracks
numTracks = length(tracksFinal);

%get track start and end times
trackSEL = getTrackSEL(tracksFinal);

%give tracks status based on the frames they span:
%2: track exists throughout movie
%1: track exists either in first frame or in last frame
%0: track does not exist in both first frame and last frame
trackStatus  = (trackSEL(:,1) == tracksFirstFrame) + (trackSEL(:,2) == tracksLastFrame);

%give all tracks same classification if lifetime classification not
%requested
if classifyLft == 0
%     trackStatus(:) = 2;
    trackStatus(:) = 0;
end

%get number of segments making each track
numSegments = zeros(numTracks,1);
for i = 1 : numTracks
    numSegments(i) = size(tracksFinal(i).tracksCoordAmpCG,1);
end

%locate the row of the first segment of each compound track in the
%big matrices of all tracks (to be constructed in the next step)
trackStartRow = ones(numTracks,1);
for iTrack = 2 : numTracks
    trackStartRow(iTrack) = trackStartRow(iTrack-1) + numSegments(iTrack-1);
end

%find total number of segments in all tracks (i.e. number of rows in big
%matrices)
numSegmentsTracks = trackStartRow(end)+numSegments(end)-1;

%construct a matrix indicating point status in big matrices:
%-2: bad gap (gap length > either segment length on its sides)
%-1: good gap (gap length < both segment lengths on its sides)
%0 : before track start or after track end
%1 : detected feature in the middle of a track of trackStatus = 0
%2 : detected feature in the middle of a track of trackStatus = 1
%3 : detected feature in the middle of a track of trackStatus = 2
%4 : detected feature just after a birth
%5 : detected feature just before a death
%6 : detected feature just after a split
%7 : detected feature just before a merge
pointStatus = zeros(numSegmentsTracks,numFrames);

%put all tracks together in one big matrix
%put the x-coordinates in one matrix and the y-coordinates in another
%indicate the status of each point
xCoordMatAll = NaN*ones(numSegmentsTracks,numFrames);
yCoordMatAll = xCoordMatAll;
for iTrack = 1 : numTracks
    
    %get track start and end times
    startTime = trackSEL(iTrack,1);
    endTime   = trackSEL(iTrack,2);
    
    %store x-coordinates
    xCoordMatAll(trackStartRow(iTrack):trackStartRow(iTrack)+...
        numSegments(iTrack)-1,startTime:endTime) = ...
        tracksFinal(iTrack).tracksCoordAmpCG(:,1:8:end);
    
    %store y-coordinates
    yCoordMatAll(trackStartRow(iTrack):trackStartRow(iTrack)+...
        numSegments(iTrack)-1,startTime:endTime) = ...
        tracksFinal(iTrack).tracksCoordAmpCG(:,2:8:end);
    
    %assign point status for features in the middle of the track
    pointStatus(trackStartRow(iTrack):trackStartRow(iTrack)+...
        numSegments(iTrack)-1,startTime:endTime) = trackStatus(iTrack) + 1;
    
    %get sequence of events of track
    seqOfEvents = tracksFinal(iTrack).seqOfEvents;
    
    if highlightES
        
        %assign point status for features just after a birth
        points2consider = find(seqOfEvents(:,2)==1 & isnan(seqOfEvents(:,4)) ...
            & seqOfEvents(:,1)~=tracksFirstFrame)';
        for iPoint = points2consider
            pointStatus(trackStartRow(iTrack)+seqOfEvents(iPoint,3)-1,...
                seqOfEvents(iPoint,1)) = 4;
        end
        
        %assign point status for features just before a death
        points2consider = find(seqOfEvents(:,2)==2 & isnan(seqOfEvents(:,4)) ...
            & seqOfEvents(:,1)~=tracksLastFrame)';
        for iPoint = points2consider
            pointStatus(trackStartRow(iTrack)+seqOfEvents(iPoint,3)-1,...
                seqOfEvents(iPoint,1)) = 5;
        end
        
    end
    
    %assign point status for features just after and before a split
    %also, in the frame just before splitting, give the
    %splitting track the position of the track it split from
    points2consider = find(seqOfEvents(:,2)==1 & ~isnan(seqOfEvents(:,4)))';
    for iPoint = points2consider
        pointStatus(trackStartRow(iTrack)+seqOfEvents(iPoint,3)-1,...
            seqOfEvents(iPoint,1)-1:seqOfEvents(iPoint,1)) = 6;
        pointStatus(trackStartRow(iTrack)+seqOfEvents(iPoint,4)-1,...
            seqOfEvents(iPoint,1)-1:seqOfEvents(iPoint,1)) = 6;
        xCoordMatAll(trackStartRow(iTrack)+seqOfEvents(iPoint,3)-1,...
            seqOfEvents(iPoint,1)-1) = xCoordMatAll(trackStartRow(iTrack)...
            +seqOfEvents(iPoint,4)-1,seqOfEvents(iPoint,1)-1);
        yCoordMatAll(trackStartRow(iTrack)+seqOfEvents(iPoint,3)-1,...
            seqOfEvents(iPoint,1)-1) = yCoordMatAll(trackStartRow(iTrack)...
            +seqOfEvents(iPoint,4)-1,seqOfEvents(iPoint,1)-1);
    end
    
    %assign point status for features just before and after a merge
    %also, in the frame just after merging, give the
    %merging track the position of the track it merged from
    points2consider = find(seqOfEvents(:,2)==2 & ~isnan(seqOfEvents(:,4)))';
    for iPoint = points2consider
        pointStatus(trackStartRow(iTrack)+seqOfEvents(iPoint,3)-1,...
            seqOfEvents(iPoint,1)-1:seqOfEvents(iPoint,1)) = 7;
        pointStatus(trackStartRow(iTrack)+seqOfEvents(iPoint,4)-1,...
            seqOfEvents(iPoint,1)-1:seqOfEvents(iPoint,1)) = 7;
        xCoordMatAll(trackStartRow(iTrack)+seqOfEvents(iPoint,3)-1,...
            seqOfEvents(iPoint,1)) = xCoordMatAll(trackStartRow(iTrack)...
            +seqOfEvents(iPoint,4)-1,seqOfEvents(iPoint,1));
        yCoordMatAll(trackStartRow(iTrack)+seqOfEvents(iPoint,3)-1,...
            seqOfEvents(iPoint,1)) = yCoordMatAll(trackStartRow(iTrack)...
            +seqOfEvents(iPoint,4)-1,seqOfEvents(iPoint,1));
    end
    
end %(for iTrack = 1 : numTracks)

%shift coordinate to account for cropped region of interest
xCoordMatAll = xCoordMatAll - (imageRange(2,1)-1);
yCoordMatAll = yCoordMatAll - (imageRange(1,1)-1);

%find gaps in tracks
gapInfo = findTrackGaps(tracksFinal);

%for gaps, assign the position as that of the feature before the gap
%also, assign a point status of -1
for iGap = 1 : size(gapInfo,1)
    
    iTrack = gapInfo(iGap,1);
    iSegment = gapInfo(iGap,2);
    iStart = gapInfo(iGap,3);
    gapLength = gapInfo(iGap,4);
    if classifyGaps
        if gapInfo(iGap,5) <= 1 && gapInfo(iGap,6) <= 1
            gapType = -1;
        else
            gapType = -2;
        end
    else
        gapType = -1;
    end
    
    xCoordMatAll(trackStartRow(iTrack)+iSegment-1,iStart:iStart+gapLength-1) = ...
        xCoordMatAll(trackStartRow(iTrack)+iSegment-1,iStart-1); %x-coordinates
    
    yCoordMatAll(trackStartRow(iTrack)+iSegment-1,iStart:iStart+gapLength-1) = ...
        yCoordMatAll(trackStartRow(iTrack)+iSegment-1,iStart-1); %y-coordinates
    
    pointStatus(trackStartRow(iTrack)+iSegment-1,iStart:iStart+gapLength-1) = gapType; %point status
    
end

%retain in the big matrices only the frames of interest
xCoordMatAll = xCoordMatAll(:,startend(1):startend(2));
yCoordMatAll = yCoordMatAll(:,startend(1):startend(2));
pointStatus = pointStatus(:,startend(1):startend(2));

%% divide tracks based on diffusion analysis or just into groups to be colored separately

%if tracks are to be individually colored
if colorTracks
    
    %divide tracks among 9 matrices that will get their own colors
    
    %x-coordinates ...
    [xCoordMatAll0,xCoordMatAll1,xCoordMatAll2,xCoordMatAll3,xCoordMatAll4,...
        xCoordMatAll5,xCoordMatAll6,xCoordMatAll7,xCoordMatAll8] = deal(NaN(size(xCoordMatAll)));
    xCoordMatAll0(1:9:end,:) = xCoordMatAll(1:9:end,:);
    xCoordMatAll1(2:9:end,:) = xCoordMatAll(2:9:end,:);
    xCoordMatAll2(3:9:end,:) = xCoordMatAll(3:9:end,:);
    xCoordMatAll3(4:9:end,:) = xCoordMatAll(4:9:end,:);
    xCoordMatAll4(5:9:end,:) = xCoordMatAll(5:9:end,:);
    xCoordMatAll5(6:9:end,:) = xCoordMatAll(6:9:end,:);
    xCoordMatAll6(7:9:end,:) = xCoordMatAll(7:9:end,:);
    xCoordMatAll7(8:9:end,:) = xCoordMatAll(8:9:end,:);
    xCoordMatAll8(9:9:end,:) = xCoordMatAll(9:9:end,:);
    
    %y-coordinates ...
    [yCoordMatAll0,yCoordMatAll1,yCoordMatAll2,yCoordMatAll3,yCoordMatAll4,...
        yCoordMatAll5,yCoordMatAll6,yCoordMatAll7,yCoordMatAll8] = deal(NaN(size(yCoordMatAll)));
    yCoordMatAll0(1:9:end,:) = yCoordMatAll(1:9:end,:);
    yCoordMatAll1(2:9:end,:) = yCoordMatAll(2:9:end,:);
    yCoordMatAll2(3:9:end,:) = yCoordMatAll(3:9:end,:);
    yCoordMatAll3(4:9:end,:) = yCoordMatAll(4:9:end,:);
    yCoordMatAll4(5:9:end,:) = yCoordMatAll(5:9:end,:);
    yCoordMatAll5(6:9:end,:) = yCoordMatAll(6:9:end,:);
    yCoordMatAll6(7:9:end,:) = yCoordMatAll(7:9:end,:);
    yCoordMatAll7(8:9:end,:) = yCoordMatAll(8:9:end,:);
    yCoordMatAll8(9:9:end,:) = yCoordMatAll(9:9:end,:);
    
else %otherwise
    
    %copy all coordinates into new variables
    %only the "0" variables will be used when there is no diffusion
    %classification
    [xCoordMatAll0,xCoordMatAll1,xCoordMatAll2,xCoordMatAll3,xCoordMatAll4,...
        xCoordMatAll5,xCoordMatAll6,xCoordMatAll7,xCoordMatAll8] = deal(xCoordMatAll);
    [yCoordMatAll0,yCoordMatAll1,yCoordMatAll2,yCoordMatAll3,yCoordMatAll4,...
        yCoordMatAll5,yCoordMatAll6,yCoordMatAll7,yCoordMatAll8] = deal(yCoordMatAll);
    
end

%if there are diffusion analysis results ...
if ~isempty(diffAnalysisRes)
    
    if transDiffClass %if transient diffusion classification ...
        
        %put all track segment classifications into one array (note that
        %here we're talking about merging and splitting segments)
        trackSegmentClass = vertcat(diffAnalysisRes.segmentClass);
        
        %get number of time points in movie
        numTimePoints = size(xCoordMatAll,2);
        
        %go over all track segments and classify their points
        for iTrackSegment = 1 : numSegmentsTracks
            
            %get current track's classification
            trackClassCurrent = trackSegmentClass(iTrackSegment).momentScalingSpectrum(:,1:3);
            
            %map the transient classification numbers into the same numbers
            %as the whole track classification
            trackClassCol3 = trackClassCurrent(:,3);
            trackClassCurrent(trackClassCol3==1,3) = 5; %random/unclassified + 2D confined (blue)
            trackClassCurrent(trackClassCol3==2,3) = 6; %random/unclassified + 2D normal (cyan)
            trackClassCurrent(trackClassCol3==3,3) = 7; %random/unclassified + 2D super (magenta)
            trackClassCurrent(isnan(trackClassCol3),3) = 0; %completely unclassified
            
            %extract and store segments classified as 5
            classCurrent = find(trackClassCurrent(:,3)==5); %find points classified as 5
            pointsCurrent = [];
            for iClass = classCurrent'
                pointsCurrent = [pointsCurrent (trackClassCurrent(iClass,1):...
                    trackClassCurrent(iClass,2)+1)]; %#ok<AGROW> %extend each interval by 1 point beyond its end to ensure continuity
            end
            pointsNotCurrent = setdiff(1:numTimePoints,pointsCurrent); %find points not classified as 5
            xCoordMatAll5(iTrackSegment,pointsNotCurrent) = NaN; %remove time points not classified as 5
            yCoordMatAll5(iTrackSegment,pointsNotCurrent) = NaN;
            
            %extract and store segments classified as 6
            classCurrent = find(trackClassCurrent(:,3)==6); %find points classified as 6
            pointsCurrent = [];
            for iClass = classCurrent'
                pointsCurrent = [pointsCurrent (trackClassCurrent(iClass,1):...
                    trackClassCurrent(iClass,2)+1)]; %#ok<AGROW> %extend each interval by 1 point beyond its end to ensure continuity
            end
            pointsNotCurrent = setdiff(1:numTimePoints,pointsCurrent); %find points not classified as 6
            xCoordMatAll6(iTrackSegment,pointsNotCurrent) = NaN; %remove time points not classified as 6
            yCoordMatAll6(iTrackSegment,pointsNotCurrent) = NaN;
            
            %extract and store segments classified as 7
            classCurrent = find(trackClassCurrent(:,3)==7); %find points classified as 7
            pointsCurrent = [];
            for iClass = classCurrent'
                pointsCurrent = [pointsCurrent (trackClassCurrent(iClass,1):...
                    trackClassCurrent(iClass,2)+1)]; %#ok<AGROW> %extend each interval by 1 point beyond its end to ensure continuity
            end
            pointsNotCurrent = setdiff(1:numTimePoints,pointsCurrent); %find points not classified as 7
            xCoordMatAll7(iTrackSegment,pointsNotCurrent) = NaN; %remove time points not classified as 7
            yCoordMatAll7(iTrackSegment,pointsNotCurrent) = NaN;
            
            %extract and store unclassified segments
            classCurrent = find(trackClassCurrent(:,3)==0); %find points classified as 0
            pointsCurrent = [];
            for iClass = classCurrent'
                pointsCurrent = [pointsCurrent (trackClassCurrent(iClass,1):...
                    trackClassCurrent(iClass,2)+1)]; %#ok<AGROW> %extend each interval by 1 point beyond its end to ensure continuity
            end
            pointsNotCurrent = setdiff(1:numTimePoints,pointsCurrent); %find points not classified as 0
            xCoordMatAll0(iTrackSegment,pointsNotCurrent) = NaN; %remove time points not classified as 0
            yCoordMatAll0(iTrackSegment,pointsNotCurrent) = NaN;
            
        end
        
    else %if whole track diffusion classification ...
        
        %get track segment types from diffusion analysis
        trackSegmentType = vertcat(diffAnalysisRes.classification);
        
        %assign classes
        trackClass = zeros(numSegmentsTracks,1); %initialize with the indicator for undetermined
        trackClass(trackSegmentType(:,1) == 1 & trackSegmentType(:,3) == 1) = 1; %linear + 1D confined (orange)
        trackClass(trackSegmentType(:,1) == 1 & trackSegmentType(:,3) == 2) = 2; %linear + 1D normal (red)
        trackClass(trackSegmentType(:,1) == 1 & trackSegmentType(:,3) == 3) = 3; %linear + 1D super (green)
        trackClass(trackSegmentType(:,1) == 1 & isnan(trackSegmentType(:,3))) = 4; %linear + too short (yellow)
        trackClass(trackSegmentType(:,1) ~= 1 & trackSegmentType(:,2) == 1) = 5; %random/unclassified + 2D confined (blue)
        trackClass(trackSegmentType(:,1) ~= 1 & trackSegmentType(:,2) == 2) = 6; %random/unclassified + 2D normal (cyan)
        trackClass(trackSegmentType(:,1) ~= 1 & trackSegmentType(:,2) == 3) = 7; %random/unclassified + 2D super (magenta)
        trackClass(trackSegmentType(:,1) == 0 & isnan(trackSegmentType(:,2))) = 8; %random + too short (purple)
        
        %extract the tracks/track segments of different classifications
        %x-coordinates
        xCoordMatAll0(trackClass~=0,:) = NaN;
        xCoordMatAll1(trackClass~=1,:) = NaN;
        xCoordMatAll2(trackClass~=2,:) = NaN;
        xCoordMatAll3(trackClass~=3,:) = NaN;
        xCoordMatAll4(trackClass~=4,:) = NaN;
        xCoordMatAll5(trackClass~=5,:) = NaN;
        xCoordMatAll6(trackClass~=6,:) = NaN;
        xCoordMatAll7(trackClass~=7,:) = NaN;
        xCoordMatAll8(trackClass~=8,:) = NaN;
        %y-coordinates
        yCoordMatAll0(trackClass~=0,:) = NaN;
        yCoordMatAll1(trackClass~=1,:) = NaN;
        yCoordMatAll2(trackClass~=2,:) = NaN;
        yCoordMatAll3(trackClass~=3,:) = NaN;
        yCoordMatAll4(trackClass~=4,:) = NaN;
        yCoordMatAll5(trackClass~=5,:) = NaN;
        yCoordMatAll6(trackClass~=6,:) = NaN;
        yCoordMatAll7(trackClass~=7,:) = NaN;
        yCoordMatAll8(trackClass~=8,:) = NaN;
        
    end %(if transDiffClass ... else ...)
    
end

%% make movie

%initialize movie if it is to be saved
if saveMovie
    movieVar = struct('cdata',[],'colormap',[]);
    movieVar = movieInfrastructure('initialize',movieType,dir2saveMovie,...
        movieName,numFramesMovie,movieVar,[]);
end

%go over all specified frames and find minimum and maximum intensity in all
%of them combined
switch intensityScale
    case 0
        intensityMinMax = [];
    case 1
        meanIntensity = zeros(size(xCoordMatAll,2),1);
        stdIntensity = meanIntensity;
        for iFrame = 1 : size(xCoordMatAll,2)
            if frame2fileMap(iFrame) ~= 0
                imageStack = double(imread(outFileList{frame2fileMap(iFrame)}));
                %                 imageStack = imageStack(imageRange(2,1):imageRange(2,2),...
                %                     imageRange(1,1):imageRange(1,2));
                imageStack = imageStack(imageRange(1,1):imageRange(1,2),...
                    imageRange(2,1):imageRange(2,2));
                
                meanIntensity(iFrame) = mean(imageStack(:));
                stdIntensity(iFrame) = std(imageStack(:));
            end
        end
        meanIntensity = mean(meanIntensity);
        stdIntensity = mean(stdIntensity);
        intensityMinMax = [meanIntensity-2*stdIntensity meanIntensity+6*stdIntensity];
    case 2
        minIntensity = zeros(size(xCoordMatAll,2),1);
        maxIntensity = minIntensity;
        for iFrame = 1 : size(xCoordMatAll,2)
            if frame2fileMap(iFrame) ~= 0
                imageStack = double(imread(outFileList{frame2fileMap(iFrame)}));
                %                 imageStack = imageStack(imageRange(2,1):imageRange(2,2),...
                %                     imageRange(1,1):imageRange(1,2));
                imageStack = imageStack(imageRange(1,1):imageRange(1,2),...
                    imageRange(2,1):imageRange(2,2));
                minIntensity(iFrame) = min(imageStack(:));
                maxIntensity(iFrame) = max(imageStack(:));
            end
        end
        minIntensity = min(minIntensity);
        maxIntensity = max(maxIntensity);
        intensityMinMax = [minIntensity maxIntensity];
end

%go over all specified frames
if plotFullScreen
    scrsz = get(0,'ScreenSize');
    h     = figure();
    set(h,'Position',scrsz);
else
    figure
end
for iFrame = 1 : numFramesMovie
    
    if frame2fileMap(iFrame) ~= 0 %if frame exists
        
        %read specified image
        imageStack = imread(outFileList{frame2fileMap(iFrame)});
        
        %filter images if requested
        if filterSigma
            imageStack = filterGauss2D(imageStack,filterSigma);
        end
        
    else %otherwise
        
        %make empty frame
        imageStack = zeros(isx,isy);
        
    end
    
    %crop to region of interest
    imageStack = imageStack(imageRange(1,1):imageRange(1,2),...
        imageRange(2,1):imageRange(2,2),:);
    
    %     tmp = double(imageStack(:,:,1));
    %     minTmp = min(tmp(:));
    %     maxTmp = max(tmp(:));
    %     tmp = (tmp - minTmp)/(maxTmp - minTmp);
    %     imageStack(:,:,1) = uint8(tmp*255);
    %
    %     tmp = double(imageStack(:,:,2));
    %     minTmp = min(tmp(:));
    %     maxTmp = max(tmp(:));
    %     tmp = (tmp - minTmp)/(maxTmp - minTmp);
    %     imageStack(:,:,2) = uint8(tmp*255);
    %     imageStack(:,:,3) = imageStack(:,:,2);
    
    %plot image in current frame and show frame number
    clf;
    switch showRaw
        case 1
            
            axes('Position',[0 0 0.495 1]);
            imshow(imageStack,intensityMinMax);
            %             xlim(imageRange(2,:));
            %             ylim(imageRange(1,:));
            hold on;
            %             plot([80 207],[50 50],'y:','LineWidth',0.5)
            %             plot([80 207],[177 177],'y:','LineWidth',0.5)
            %             plot([80 80],[50 177],'y:','LineWidth',0.5)
            %             plot([207 207],[50 177],'y:','LineWidth',0.5)
            textDeltaCoord = min(diff(imageRange,[],2))/20;
            %             text(imageRange(2,1)+textDeltaCoord,imageRange(1,1)+...
            %                 textDeltaCoord,num2str(iFrame+startend(1)-1),...
            %                 'Color','white','FontSize',18);
            %             text(textDeltaCoord,...
            %                 textDeltaCoord,num2str(iFrame+startend(1)-1),...
            %                 'Color','white','FontSize',18);
            text(textDeltaCoord,...
                textDeltaCoord,[num2str(((iFrame+startend(1)-1)-1)*0.1,'%4.1f') ' s'],'Color','yellow');
            axes('Position',[0.505 0 0.495 1]);
            imshow(imageStack,intensityMinMax);
            %             xlim(imageRange(2,:));
            %             ylim(imageRange(1,:));
            hold on;
            %             plot([80 207],[50 50],'y:','LineWidth',0.5)
            %             plot([80 207],[177 177],'y:','LineWidth',0.5)
            %             plot([80 80],[50 177],'y:','LineWidth',0.5)
            %             plot([207 207],[50 177],'y:','LineWidth',0.5)
            
        case 2
            axes('Position',[0 0.505 1 0.495]);
            imshow(imageStack,intensityMinMax);
            %             xlim(imageRange(2,:));
            %             ylim(imageRange(1,:));
            hold on;
            textDeltaCoord = min(diff(imageRange,[],2))/20;
            %             text(imageRange(2,1)+textDeltaCoord,imageRange(1,1)+...
            %                 textDeltaCoord,num2str(iFrame+startend(1)-1),...
            %                 'Color','white','FontSize',18);
            text(textDeltaCoord,...
                textDeltaCoord,num2str(iFrame+startend(1)-1),...
                'Color','white','FontSize',18);
            %             text(textDeltaCoord-1,...
            %                 textDeltaCoord+2,[num2str((iFrame+startend(1)-2)*0.025,'%5.3f') ' s'],...
            %                 'Color','white','FontSize',18);
            axes('Position',[0 0 1 0.495]);
            imshow(imageStack,intensityMinMax);
            %             xlim(imageRange(2,:));
            %             ylim(imageRange(1,:));
            hold on;
        otherwise
            axes('Position',[0 0 1 1]);
            imshow(imageStack,intensityMinMax);
            %             xlim(imageRange(2,:));
            %             ylim(imageRange(1,:));
            hold on;
            textDeltaCoord = min(diff(imageRange,[],2))/20;
            %             text(imageRange(2,1)+textDeltaCoord,imageRange(1,1)+...
            %                 textDeltaCoord,num2str(iFrame+startend(1)-1),...
            %                 'Color','white','FontSize',18);
            text(textDeltaCoord,...
                textDeltaCoord,num2str(iFrame+startend(1)-1),...
                'Color','white','FontSize',18);
    end
    
    %get tracks to plot
    plotOrNot = 0;
    if dragtailLength >= 0 %to plot tracks dynamically
        
        if iFrame > 1 || onlyTracks
            dragTailStart = max(iFrame-dragtailLength,1);
            indx2keep = find(pointStatus(:,iFrame)~=0);
            xCoord2plot0 = (xCoordMatAll0(indx2keep,dragTailStart:iFrame))';
            yCoord2plot0 = (yCoordMatAll0(indx2keep,dragTailStart:iFrame))';
            xCoord2plot1 = (xCoordMatAll1(indx2keep,dragTailStart:iFrame))';
            yCoord2plot1 = (yCoordMatAll1(indx2keep,dragTailStart:iFrame))';
            xCoord2plot2 = (xCoordMatAll2(indx2keep,dragTailStart:iFrame))';
            yCoord2plot2 = (yCoordMatAll2(indx2keep,dragTailStart:iFrame))';
            xCoord2plot3 = (xCoordMatAll3(indx2keep,dragTailStart:iFrame))';
            yCoord2plot3 = (yCoordMatAll3(indx2keep,dragTailStart:iFrame))';
            xCoord2plot4 = (xCoordMatAll4(indx2keep,dragTailStart:iFrame))';
            yCoord2plot4 = (yCoordMatAll4(indx2keep,dragTailStart:iFrame))';
            xCoord2plot5 = (xCoordMatAll5(indx2keep,dragTailStart:iFrame))';
            yCoord2plot5 = (yCoordMatAll5(indx2keep,dragTailStart:iFrame))';
            xCoord2plot6 = (xCoordMatAll6(indx2keep,dragTailStart:iFrame))';
            yCoord2plot6 = (yCoordMatAll6(indx2keep,dragTailStart:iFrame))';
            xCoord2plot7 = (xCoordMatAll7(indx2keep,dragTailStart:iFrame))';
            yCoord2plot7 = (yCoordMatAll7(indx2keep,dragTailStart:iFrame))';
            xCoord2plot8 = (xCoordMatAll8(indx2keep,dragTailStart:iFrame))';
            yCoord2plot8 = (yCoordMatAll8(indx2keep,dragTailStart:iFrame))';
            plotOrNot = 1;
        end
        
    elseif dragtailLength == -1 %to plot tracks statically
        
        xCoord2plot0 = (xCoordMatAll0)';
        yCoord2plot0 = (yCoordMatAll0)';
        xCoord2plot1 = (xCoordMatAll1)';
        yCoord2plot1 = (yCoordMatAll1)';
        xCoord2plot2 = (xCoordMatAll2)';
        yCoord2plot2 = (yCoordMatAll2)';
        xCoord2plot3 = (xCoordMatAll3)';
        yCoord2plot3 = (yCoordMatAll3)';
        xCoord2plot4 = (xCoordMatAll4)';
        yCoord2plot4 = (yCoordMatAll4)';
        xCoord2plot5 = (xCoordMatAll5)';
        yCoord2plot5 = (yCoordMatAll5)';
        xCoord2plot6 = (xCoordMatAll6)';
        yCoord2plot6 = (yCoordMatAll6)';
        xCoord2plot7 = (xCoordMatAll7)';
        yCoord2plot7 = (yCoordMatAll7)';
        xCoord2plot8 = (xCoordMatAll8)';
        yCoord2plot8 = (yCoordMatAll8)';
        plotOrNot = 1;
        
    elseif dragtailLength == -2 %to plot tracks dynamically but keep them after they disappear
        
        if iFrame > 1 || onlyTracks
            xCoord2plot0 = xCoordMatAll0(:,1:iFrame)';
            yCoord2plot0 = yCoordMatAll0(:,1:iFrame)';
            xCoord2plot1 = xCoordMatAll1(:,1:iFrame)';
            yCoord2plot1 = yCoordMatAll1(:,1:iFrame)';
            xCoord2plot2 = xCoordMatAll2(:,1:iFrame)';
            yCoord2plot2 = yCoordMatAll2(:,1:iFrame)';
            xCoord2plot3 = xCoordMatAll3(:,1:iFrame)';
            yCoord2plot3 = yCoordMatAll3(:,1:iFrame)';
            xCoord2plot4 = xCoordMatAll4(:,1:iFrame)';
            yCoord2plot4 = yCoordMatAll4(:,1:iFrame)';
            xCoord2plot5 = xCoordMatAll5(:,1:iFrame)';
            yCoord2plot5 = yCoordMatAll5(:,1:iFrame)';
            xCoord2plot6 = xCoordMatAll6(:,1:iFrame)';
            yCoord2plot6 = yCoordMatAll6(:,1:iFrame)';
            xCoord2plot7 = xCoordMatAll7(:,1:iFrame)';
            yCoord2plot7 = yCoordMatAll7(:,1:iFrame)';
            xCoord2plot8 = xCoordMatAll8(:,1:iFrame)';
            yCoord2plot8 = yCoordMatAll8(:,1:iFrame)';
            plotOrNot = 1;
        end
        
    end
    
    %plot tracks
    if plotOrNot
        
        %plot basic tracks
        plot([xCoord2plot0(1,:); xCoord2plot0],[yCoord2plot0(1,:); yCoord2plot0],...
            'Color',[1 0.7 0.7],'LineWidth',1); %light pink; the artificial repetition of the first line is for avoiding a mess in the first frame when tracks are not color-coded individually
                         
        %color individual tracks randomly if requested
        if colorTracks == 1
            plot(xCoord2plot1,yCoord2plot1,'Color',[1 0.7 0],'LineWidth',1); %orange
            plot(xCoord2plot2,yCoord2plot2,'Color','r','LineWidth',1); %[1 0 0]
            plot(xCoord2plot3,yCoord2plot3,'Color','g','LineWidth',1); %[0 1 0]
            plot(xCoord2plot4,yCoord2plot4,'Color','y','LineWidth',1); %[1 1 0]
            plot(xCoord2plot5,yCoord2plot5,'Color','b','LineWidth',1); %[0 0 1]
            plot(xCoord2plot6,yCoord2plot6,'Color','c','LineWidth',1); %[0 1 1]
            plot(xCoord2plot7,yCoord2plot7,'Color','m','LineWidth',1); %[1 0 1]
            plot(xCoord2plot8,yCoord2plot8,'Color',[0.6 0 1],'LineWidth',1); %purple
        end
        
        %color-code dragtail based on diffusion analysis if supplied
        if ~isempty(diffAnalysisRes)
            plot(xCoord2plot5,yCoord2plot5,'Color','b','LineWidth',1); %[0 0 1]
            plot(xCoord2plot6,yCoord2plot6,'Color','c','LineWidth',1); %[0 1 1]
            plot(xCoord2plot7,yCoord2plot7,'Color','m','LineWidth',1); %[1 0 1]
            if ~transDiffClass
                plot(xCoord2plot1,yCoord2plot1,'Color',[1 0.7 0],'LineWidth',1); %orange
                plot(xCoord2plot2,yCoord2plot2,'Color','r','LineWidth',1); %[1 0 0]
                plot(xCoord2plot3,yCoord2plot3,'Color','g','LineWidth',1); %[0 1 0]
                plot(xCoord2plot4,yCoord2plot4,'Color','y','LineWidth',1); %[1 1 0]
                plot(xCoord2plot8,yCoord2plot8,'Color',[0.6 0 1],'LineWidth',1); %purple
            end
        end
        
    end
    
    %plot points (features + gaps + merges + splits)
    if ~onlyTracks
        
% % % % % % % % % % % % % % % % % % % % % %         %blue stars: bad gaps
% % % % % % % % % % % % % % % % % % % % % %         points2plot = find(pointStatus(:,iFrame)==-2);
% % % % % % % % % % % % % % % % % % % % % %         plot(xCoordMatAll(points2plot,iFrame),yCoordMatAll(points2plot,iFrame),'b*','MarkerSize',6);
% % % % % % % % % % % % % % % % % % % % % %         %         plot(xCoordMatAll(points2plot,iFrame),yCoordMatAll(points2plot,iFrame),'b*','MarkerSize',2);
% % % % % % % % % % % % % % % % % % % % % %         
% % % % % % % % % % % % % % % % % % % % % %         %cyan stars: good gaps
        points2plot = find(pointStatus(:,iFrame)==-1);
        if isempty(diffAnalysisRes)
            %             plot(xCoordMatAll(points2plot,iFrame),yCoordMatAll(points2plot,iFrame),'c*','MarkerSize',6);
            plot(xCoordMatAll(points2plot,iFrame),yCoordMatAll(points2plot,iFrame),'co','MarkerSize',6,'LineWidth',2);
            %             plot(xCoordMatAll(points2plot,iFrame),yCoordMatAll(points2plot,iFrame),'c*','MarkerSize',10,'LineWidth',2);
        else
            plot(xCoordMatAll(points2plot,iFrame),yCoordMatAll(points2plot,iFrame),'w*','MarkerSize',6);
        end
% % % % % % % % % % % % % % % % % % % % % %         
        %red circles: detected feature in the middle of track with status 0
% % % % % % % % % % % % % % %          points2plot = find(pointStatus(:,iFrame)==1);
% % % % % % % % % % % % % % %         plot(xCoordMatAll(points2plot,iFrame),yCoordMatAll(points2plot,iFrame),'ro','MarkerSize',8);
        %         plot(xCoordMatAll(points2plot,iFrame),yCoordMatAll(points2plot,iFrame),'wo','MarkerSize',2);
% % % % % % % % % % % % % % % % % % % % % % %         
% % % % % % % % % % % % % % % % % % % % % % %         %magenta circles: detected feature in the middle of track with status 1
% % % % % % % % % % % % % % % % % % % % % % %         points2plot = find(pointStatus(:,iFrame)==2);
% % % % % % % % % % % % % % % % % % % % % % %         plot(xCoordMatAll(points2plot,iFrame),yCoordMatAll(points2plot,iFrame),'mo','MarkerSize',6);
% % % % % % % % % % % % % % % % % % % % % % %         %         plot(xCoordMatAll(points2plot,iFrame),yCoordMatAll(points2plot,iFrame),'mo','MarkerSize',2);
% % % % % % % % % % % % % % % % % % % % % % %         %         plot(xCoordMatAll(points2plot,iFrame),yCoordMatAll(points2plot,iFrame),'mo','MarkerSize',10,'LineWidth',2);
% % % % % % % % % % % % % % % % % % % % % % %         
% % % % % % % % % % % % % % % % % % % % % % %         %white circles: detected feature in the middle of track with status 2
% % % % % % % % % % % % % % % % % % % % % % %         points2plot = find(pointStatus(:,iFrame)==3);
% % % % % % % % % % % % % % % % % % % % % % %         plot(xCoordMatAll(points2plot,iFrame),yCoordMatAll(points2plot,iFrame),'wo','MarkerSize',6);
% % % % % % % % % % % % % % % % % % % % % % %         %         plot(xCoordMatAll(points2plot,iFrame),yCoordMatAll(points2plot,iFrame),'wo','MarkerSize',5);
% % % % % % % % % % % % % % % % % % % % % % %         %         plot(xCoordMatAll(points2plot,iFrame),yCoordMatAll(points2plot,iFrame),'wo','MarkerSize',2);
% % % % % % % % % % % % % % % % % % % % % % %         
% % % % % % % % % % % % % % % % % % % % % % %         %green circles: detected feature just after birth
% % % % % % % % % % % % % % % % % % % % % % %         points2plot = find(pointStatus(:,iFrame)==4);
% % % % % % % % % % % % % % % % % % % % % % %         plot(xCoordMatAll(points2plot,iFrame),yCoordMatAll(points2plot,iFrame),'go','MarkerSize',6);
% % % % % % % % % % % % % % % % % % % % % % %         %         plot(xCoordMatAll(points2plot,iFrame),yCoordMatAll(points2plot,iFrame),'go','MarkerSize',2);
% % % % % % % % % % % % % % % % % % % % % % %         %         plot(xCoordMatAll(points2plot,iFrame),yCoordMatAll(points2plot,iFrame),'go','MarkerSize',15,'LineWidth',2);
        
% % % % % %         %yellow circles: detected feature just before death
% % % % % %         points2plot = find(pointStatus(:,iFrame)==5);
% % % % % %         plot(xCoordMatAll(points2plot,iFrame),yCoordMatAll(points2plot,iFrame),'yo','MarkerSize',6);
% % % % % %         %         plot(xCoordMatAll(points2plot,iFrame),yCoordMatAll(points2plot,iFrame),'yo','MarkerSize',2);
        %         plot(xCoordMatAll(points2plot,iFrame),yCoordMatAll(points2plot,iFrame),'yo','MarkerSize',15,'LineWidth',2);
        
       %% SPLIT  %green diamonds: detected feature just before/after a split
        points2plot = find(pointStatus(:,iFrame)==6);
        if isempty(diffAnalysisRes)
% % % % % %                         plot(xCoordMatAll(points2plot,iFrame),yCoordMatAll(points2plot,iFrame),'gd','MarkerSize',10);
% % % % % %             plot(xCoordMatAll(points2plot,iFrame),yCoordMatAll(points2plot,iFrame),'gd','MarkerSize',4);

plot(xCoordMatAll(points2plot,iFrame),yCoordMatAll(points2plot,iFrame),'gd','MarkerSize',6,'LineWidth',2);
            %             plot(xCoordMatAll(points2plot,iFrame),yCoordMatAll(points2plot,iFrame),'gd','MarkerSize',15,'LineWidth',2);
        else
            plot(xCoordMatAll(points2plot,iFrame),yCoordMatAll(points2plot,iFrame),'wd','MarkerSize',10);
            %             plot(xCoordMatAll(points2plot,iFrame),yCoordMatAll(points2plot,iFrame),'wd','MarkerSize',15,'LineWidth',2);
        end
        
     %% MERGE   %yellow diamonds: detected feature just before/after a merge
        points2plot = find(pointStatus(:,iFrame)==7);
        if isempty(diffAnalysisRes)
            %             plot(xCoordMatAll(points2plot,iFrame),yCoordMatAll(points2plot,iFrame),'yd','MarkerSize',10);
% % % % % % % %             plot(xCoordMatAll(points2plot,iFrame),yCoordMatAll(points2plot,iFrame),'yd','MarkerSize',6);
            plot(xCoordMatAll(points2plot,iFrame),yCoordMatAll(points2plot,iFrame),'yd','MarkerSize',6,'LineWidth',2);
            %             plot(xCoordMatAll(points2plot,iFrame),yCoordMatAll(points2plot,iFrame),'yd','MarkerSize',15,'LineWidth',2);
        else
            plot(xCoordMatAll(points2plot,iFrame),yCoordMatAll(points2plot,iFrame),'wd','MarkerSize',10);
            %             plot(xCoordMatAll(points2plot,iFrame),yCoordMatAll(points2plot,iFrame),'wd','MarkerSize',15,'LineWidth',2);
        end
        
    end
    
    %add frame to movie if movie is saved
    if saveMovie
        movieVar = movieInfrastructure('addFrame',movieType,dir2saveMovie,...
            movieName,numFramesMovie,movieVar,iFrame);
    end
    
    %pause for a moment to see frame
    pause(0.1);
    
end %(for iFrame = 1 : numFramesMovie)

%finish movie
if saveMovie
    movieInfrastructure('finalize',movieType,dir2saveMovie,...
        movieName,numFramesMovie,movieVar,[]);
end

%% ~~~ end ~~~

