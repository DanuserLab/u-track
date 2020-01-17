function overlayFeaturesMovie(movieInfo,startend,saveMovie,movieName,...
    filterSigma,showRaw,intensityScale,firstImageFile,dir2saveMovie,...
    movieType,plotFullScreen,iChannel)
%OVERLAYFEATURESMOVIE makes a movie of detected features overlaid on images
%
%SYNPOSIS overlayFeaturesMovie(movieInfo,startend,saveMovie,movieName,...
%    filterSigma,showRaw,intensityScale,firstImageFile,dir2saveMovie,...
%    movieType,plotFullScreen)
%
%INPUT  movieInfo   : Output of detectSubResFeatures2D_StandAlone.
%       startend    : Row vector indicating first and last frame to
%                     include in movie. Format: [startframe endframe].
%                     Optional. Default: [1 (maximum available frame)]
%       saveMovie   : 1 to save movie (as Quicktime), 0 otherwise.
%                     Optional. Default: 0
%       movieName   : filename for saving movie.
%                     Optional. Default: FeaturesMovie (if saveMovie = 1).
%       filterSigma : 0 to overlay on raw image, PSF sigma to overlay on image
%                     filtered with given filterSigma.
%                     Optional. Default: 0
%       showRaw     : 1 to add raw movie to the left of the movie with
%                     tracks overlaid, 2 to add raw movie at the top of
%                     the movie with tracks overlaid, 0 otherwise.
%                     Optional. Default: 0.
%       intensityScale: 0 to autoscale every image in the movie, 1
%                     to have a fixed scale using intensity mean and std, 2
%                     to have a fixed scale using minimum and maximum
%                     intensities.
%                     Optional. Default: 1.
%       firstImageFile: Name of the first image file in the folder of
%                     images that should be overlaid. The file has to be
%                     the first image that has been analyzed even if not
%                     plotted. If file is not specified [], user will be
%                     prompted to select the first image.
%                     Optional. Default: [].
%       dir2saveMovie: Directory where to save output movie.
%                     If not input, movie will be saved in directory where
%                     images are located.
%                     Optional. Default: [].
%       movieType   : 'mov' to make a Quicktime movie using MakeQTMovie,
%                     'avi' to make AVI movie using Matlab's movie2avi,
%                     'mp4_unix', 'avi_unix' to make an MP4 or AVI movie
%                     using ImageMagick and ffmpeg. These options works
%                     only under linux or mac.
%                     Optional. Default: 'mov'.
%       plotFullScreen: 1 the figure will be sized to cover the whole
%                       screen. In this way the movie will be of highest
%                       possible quality. default is 0.
%
%       iChannel: (default: channel=1) If analyzing multi-channel image,
%       select which channel detections should be overlaid on.
%
%OUTPUT the movie.
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
    disp('--overlayFeaturesMovie: Incorrect number of input arguments!');
    return
end

%check if specific channel should be used
if nargin < 12 || isempty(iChannel)
    iChannel = 1;
end


%ask user for images
if nargin < 8 || isempty(firstImageFile)
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
    
    %determine which frames the files correspond to, and generate the inverse map
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
    currentImage = imread(outFileList{1},iChannel);
    [isx,isy] = size(currentImage);
    
else %else, exit
    
    disp('--overlayFeaturesMovie: Bad file selection');
    return
    
end

%check startend and assign default if necessary
if nargin < 2 || isempty(startend)
    startend = [1 numFrames];
else
    startend(2) = min(startend(2),numFrames); %make sure that last frame does not exceed real last frame
end

%keep only the frames of interest
outFileList = outFileList(frame2fileMap(startend(1)):frame2fileMap(startend(2)));
frame2fileMap = frame2fileMap(startend(1):startend(2));
indxNotZero = find(frame2fileMap~=0);
frame2fileMap(indxNotZero) = frame2fileMap(indxNotZero) - frame2fileMap(indxNotZero(1)) + 1;

%retain only the movieInfo of the frames of interest
if isempty(movieInfo)
    movieInfo = repmat(struct('xCoord',[],'yCoord',[],'amp',[]),...
        startend(2)-startend(1)+1,1);
else
    movieInfo = movieInfo(startend(1):startend(2));
end

%get number of frames in movie to be made
numFramesMovie = diff(startend) + 1;

%get image size
imageRange = [1 isx; 1 isy];

%% input - additional parameters

%check whether to save movie
if nargin < 3 || isempty(saveMovie)
    saveMovie = 0;
end

%check name for saving movie
if saveMovie && (nargin < 4 || isempty(movieName))
    movieName = 'featuresMovie';
end

%check whether to use filtered images
if nargin < 5 || isempty(filterSigma)
    filterSigma = 0;
end

%check whether to put raw movie adjacent to movie with tracks features
if nargin < 6 || isempty(showRaw)
    showRaw = 0;
end

%check how to scale image intensity
if nargin < 7 || isempty(intensityScale)
    intensityScale = 1;
end

%check where to save resulting movie
if saveMovie && (nargin < 9 || isempty(dir2saveMovie))
    dir2saveMovie = dirName;
end

%decide on movie type
if nargin < 10 || isempty(movieType)
    movieType = 'mov';
end

%check whether to use full screen for plotting
if nargin < 11 || isempty(plotFullScreen)
    plotFullScreen = 0;
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
        meanIntensity = zeros(numFramesMovie,1);
        stdIntensity = meanIntensity;
        for iFrame = 1 : numFramesMovie
            if frame2fileMap(iFrame) ~= 0
                imageStack = double(imread(outFileList{frame2fileMap(iFrame)},iChannel));
                meanIntensity(iFrame) = mean(imageStack(:));
                stdIntensity(iFrame) = std(imageStack(:));
            end
        end
        meanIntensity = mean(meanIntensity);
        stdIntensity = mean(stdIntensity);
        intensityMinMax = [meanIntensity-2*stdIntensity meanIntensity+6*stdIntensity];
    case 2
        minIntensity = zeros(numFramesMovie,1);
        maxIntensity = minIntensity;
        for iFrame = 1 : numFramesMovie
            if frame2fileMap(iFrame) ~= 0
                imageStack = double(imread(outFileList{frame2fileMap(iFrame)},iChannel)); 
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
        imageStack = imread(outFileList{frame2fileMap(iFrame)},iChannel);
        
        %filter images if requested
        if filterSigma
            imageStack = filterGauss2D(imageStack,filterSigma);
        end
        
    else %otherwise
        
        %make empty frame
        imageStack = zeros(isx,isy);
        
    end
    
    %plot image in current frame
    clf;
    
    switch showRaw
        case 1
            axes('Position',[0 0 0.495 1]);
            imshow(imageStack,intensityMinMax);
            xlim(imageRange(2,:));
            ylim(imageRange(1,:));
            hold on;
            textDeltaCoord = min(diff(imageRange,[],2))/20;
            text(imageRange(1,1)+textDeltaCoord,imageRange(2,1)+...
                textDeltaCoord,num2str(iFrame+startend(1)-1),'Color','white');
            %             text(imageRange(1,1)+textDeltaCoord,imageRange(2,1)+...
            %                 textDeltaCoord,[num2str(((iFrame+startend(1)-1)-1)*0.025,'%7.3f') ' s'],'Color','white');
            %             plot([80 207],[50 50],'y:','LineWidth',0.5)
            %             plot([80 207],[177 177],'y:','LineWidth',0.5)
            %             plot([80 80],[50 177],'y:','LineWidth',0.5)
            %             plot([207 207],[50 177],'y:','LineWidth',0.5)
            axes('Position',[0.505 0 0.495 1]);
            imshow(imageStack,intensityMinMax);
            xlim(imageRange(2,:));
            ylim(imageRange(1,:));
            hold on;
            %             plot([80 207],[50 50],'y:','LineWidth',0.5)
            %             plot([80 207],[177 177],'y:','LineWidth',0.5)
            %             plot([80 80],[50 177],'y:','LineWidth',0.5)
            %             plot([207 207],[50 177],'y:','LineWidth',0.5)
        case 2
            axes('Position',[0 0.505 1 0.495]);
            imshow(imageStack,intensityMinMax);
            xlim(imageRange(2,:));
            ylim(imageRange(1,:));
            hold on;
            textDeltaCoord = min(diff(imageRange,[],2))/20;
            text(imageRange(1,1)+textDeltaCoord,imageRange(2,1)+...
                textDeltaCoord,num2str(iFrame+startend(1)-1),'Color','white');
            axes('Position',[0 0 1 0.495]);
            imshow(imageStack,intensityMinMax);
            xlim(imageRange(2,:));
            ylim(imageRange(1,:));
            hold on;
        otherwise
            axes('Position',[0 0 1 1]);
            imshow(imageStack,intensityMinMax);
            xlim(imageRange(2,:));
            ylim(imageRange(1,:));
            hold on;
            textDeltaCoord = min(diff(imageRange,[],2))/20;
            text(imageRange(1,1)+textDeltaCoord,imageRange(2,1)+...
                textDeltaCoord,num2str(iFrame+startend(1)-1),'Color','white');
    end
    
    %plot features
    if ~isempty(movieInfo(iFrame).xCoord)
        plot(movieInfo(iFrame).xCoord(:,1),movieInfo(iFrame).yCoord(:,1),'yo','MarkerSize',6);
    end
    
    %add frame to movie if movie is saved
    if saveMovie
        movieVar = movieInfrastructure('addFrame',movieType,dir2saveMovie,...
            movieName,numFramesMovie,movieVar,iFrame);
    end
    
    %pause for a moment to see frame
    pause(0.1);
    
end

%finish movie
if saveMovie
    movieInfrastructure('finalize',movieType,dir2saveMovie,...
        movieName,numFramesMovie,movieVar,[]);
end

%% ~~~ end ~~~

