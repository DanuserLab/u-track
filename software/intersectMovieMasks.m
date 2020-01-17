function movieData = intersectMovieMasks(movieDataOrProcess,varargin)
%INTERSECTMOVIEMASKS Spatially transforms the masks of the input movie
% 
% movieData = transformMovieMasks(movieData)
% 
% movieData = transformMovieMasks(movieData,paramsIn)
% 
% 
% This function performs a spatial transformation on the masks of the
% selected channels of the input movie.
% 
%
% Input:
% 
%   movieData - The movieData object describing the movie, as created using
%   setupMovieDataGUI.m
% 
%   paramsIn - Structure with inputs for optional parameters. The
%   parameters should be stored as fields in the structure, with the field
%   names and possible values as described below:
% 
%   Possible Parameter Structure Field Names:
%       ('FieldName' -> possible values)
%
%       ('OutputDirectory' -> character string)
%       Optional. A character string specifying the directory to save the
%       masks to. Masks for different channels will be saved as
%       sub-directories of this directory.
%       If not input, the masks will be saved to the same directory as the
%       movieData, in a sub-directory called "masks"
%
%       ('ChannelIndex'-> Positive integer scalar or vector) The integer
%       indices of the channel(s) to perform mask transformation on. This
%       index corresponds to the channel's location in the array
%       movieData.channels_. If not input, the user will be asked to select
%       from the available channels
%
%       ('SegProcessIndex' -> Positive integer scalar or vector) Optional.
%       This specifies MaskProcess(s) to use masks from by its
%       index in the array movieData.processes_; For each channel, masks
%       will be used from the last process specified which has valid masks
%       for that channel. That is if SegProcessIndex = [1 4] and both
%       processes 1 and 4 have masks for a given channel, then the masks
%       from process 4 will be used. If not input, and multiple
%       MaskProcesses are present, the user will be asked to select
%       one, unless batch mode is enabled in which case an error will be
%       generated.
%
% 
% Output:
% 
%   movieData - The modified movieData object, with the mask transform
%   logged in it, including all parameters used.
% 
%   The transformed masks will be written to a sub-directory of the
%   OutputDirectory.
% 
% 
% Sebastien Besson, Oct 2011
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

%% ----------- Input --------- %%

%Check input
ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('movieDataOrProcess', @isProcessOrMovieData);
ip.addOptional('paramsIn',[], @isstruct);
ip.parse(movieDataOrProcess,varargin{:});
paramsIn=ip.Results.paramsIn;

% Get MovieData object and Process
[movieData, maskIntProc,iProc] = getOwnerAndProcess(movieDataOrProcess,'MaskIntersectionProcess',true);

%Parse input, store in parameter structure
p = parseProcessParams(maskIntProc,paramsIn);

%% --------------- Initialization ---------------%%
if feature('ShowFigureWindows')
    wtBar = waitbar(0,'Initializing...','Name',maskIntProc.getName());
else 
    wtBar=-1;
end

if isempty(p.SegProcessIndex)    
    if ~feature('ShowFigureWindows')
        %If batch mode, just get all the seg processes
        p.SegProcessIndex = movieData.getProcessIndex('MaskProcess',Inf,0);            
        p.SegProcessIndex(p.SegProcessIndex == iProc) = [];
        if numel(p.SegProcessIndex) > 1
            error('In batch mode you must specify the SegProcessIndex if more than one MaskProcess is available!')
        end            
    else        
        %We need to exclude this function's process, and ask user if more
        %than one
        segProcList =  movieData.getProcessIndex('MaskProcess',Inf,0);
        segProcList(segProcList == iProc) = []; %Don't count this process
        iSegProc=1;
        if numel(segProcList) > 1
            procNames = cellfun(@(x)(x.getName),...
                        movieData.processes_(segProcList),'UniformOutput',false);
            iSegProc = listdlg('ListString',procNames,...
                               'SelectionMode','multiple',...
                               'ListSize',[400 400],...
                               'PromptString','Select the segmentation process(es) to use:');
            
        end
        p.SegProcessIndex = segProcList(iSegProc);        
    end
end

if isempty(p.SegProcessIndex) 
    error('This function requires that the input movie has already been segmented and that a valid MaskProcesses be specified!')
end

nChan = numel(p.ChannelIndex);
if isscalar(p.SegProcessIndex),
    p.SegProcessIndex = p.SegProcessIndex *ones(nChan,1);
else
    assert(numel(p.SegProcessIndex) == nChan);    
end


hasMasks = false(nChan,1);
segProc=cell(nChan,1);
%Check every specified process for masks
for j = 1:nChan

    %Make sure the specified process is a MaskProcess
    if ~isa(movieData.processes_{p.SegProcessIndex(j)},'MaskProcess')
        error(['The process specified by SegProcessIndex(' num2str(j) ') is not a MaskProcess!'])
    end
    
    segProc{j} = movieData.processes_{p.SegProcessIndex(j)};
    %Check which channels have masks from this process
    hasMasks(j) = segProc{j}.checkChannelOutput(p.ChannelIndex(j));        
end

%Make sure all the selected channels have foreground masks.
if ~all(hasMasks)
    error(['The mask process has not been run for all selected channels! '...
        'Please apply a valid mask process before running mask intersection!']);
end

       
% Set up the input directories
inFilePaths = cell(1,numel(movieData.channels_));
for j = 1:nChan
    iChan = p.ChannelIndex(j);
    inFilePaths{1,iChan} = segProc{j}.outFilePaths_{1,iChan};
end
maskIntProc.setInFilePaths(inFilePaths);
    
% Set up the output directory
outFilePaths = cell(1,numel(movieData.channels_));
for j = p.ChannelIndex
    outFilePaths{1,j} =  p.OutputDirectory;
end
mkClrDir(p.OutputDirectory);
maskIntProc.setOutFilePaths(outFilePaths);


%% ---------------- Mask Intersection --------------- %%

disp('Creating mask intersection...')
if ishandle(wtBar), waitbar(0,wtBar,'Creating mask intersection...'); end
   
% Reading various constants
nFrames = movieData.nFrames_;

%Create string for current directory
maskNames = arrayfun(@(i)segProc{i}.getOutMaskFileNames(p.ChannelIndex(i)),1:nChan);
inMask=@(i,frame) [segProc{i}.outFilePaths_{p.ChannelIndex(i)} filesep maskNames{i}{frame}];

%Format string for zero-padding file names
fString = ['%0' num2str(floor(log10(nFrames))+1) '.f'];
numStr = @(frame) num2str(frame,fString);
% Anonymous functions for reading input/output
outMask=@(frame) [p.OutputDirectory filesep 'intersected_mask_' numStr(frame) '.tif'];

% Redefine mask names
for j= 1:nFrames
    mask = true(movieData.imSize_);
    for i = 1:numel(p.ChannelIndex)
        mask = mask & logical(imread(inMask(i,j)));
    end
    imwrite(mask,outMask(j));

    if mod(j,5)==1 && ishandle(wtBar)
        waitbar(j/nFrames,wtBar);
    end  
end

if ishandle(wtBar), close (wtBar); end
disp('Finished intersecting masks!')
