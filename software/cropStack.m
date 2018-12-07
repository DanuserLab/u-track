function [area, output_path]=cropStack(area,output_path,first_image_path,num_images_to_crop)
% cropStack lets the user crop a region of interest from an image stack
%
% SYNOPSIS   area=cropStackk(area,output_path,first_image_path,num_images_to_crop)
%
% INPUT      area [y0 x0 y x]
%                           y0 : y coordinate of the top-left corner
%                           x0 : x coordinate of the top-left corner
%                           y  : y coordinate of the bottom-right corner
%                           x  : x coordinate of the bottom-right corner
%            Pass area=[] to manually draw a region of interest
%
%            output_path (optional) output path - if path is not passed, the
%            user will be prompted to select the path through a save
%            dialog. If a directory is specified which does not exist, 
%            it will be created (if possible).
%
%            first_image_path (optional) is the path to the first image to
%            be cropped. If path is not passed, the user will be prompted
%            to select the path through a save dialog.
%
%            num_images_to_crop (optional). 
%            num_images_to_crop=[] : will prompt you to enter the number of images to
%                                    be cropped.
%            num_images_to_crop=0  : all images will be cropped.
%            num_images_to_crop=n  : n images will be cropped. 
%
% Copyright (C) 2018, Danuser Lab - UTSouthwestern 
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

% Check input parameters
if nargin==0
   error('Not enough input arguments');
elseif nargin==1 || isempty(output_path)
   manualSave=1;
else
   manualSave=0;
end

% Creating output directory if needed
if manualSave==0 & exist(output_path)~=7
    % Directory does not exist - create it
    if output_path(2)==':' % Compatibility with Windows OS
        % Drive letter specified
        mkdir(output_path(1:3),output_path(4:end));
    else
        mkdir(output_path);
    end
    fprintf(1,'Directory %s successfully created.\n',output_path);
end

if nargin<3 || isempty(first_image_path)
    % Load First image
    [fName,dirName] = uigetfile(...
        {'*.tif;*.tiff;*.jpg;*.jpeg;*.TIF','Image Files (*.tif,*.tiff,*.jpg,*.jpeg)';
        '*.tif','TIF files (*.tif)'
        '*.tiff','TIFF files (*.tiff)'
        '*.jpg;','JPG files (*.jpg)'
        '*.jpeg;','JPEG files (*.jpeg)'
        '*.*','All Files (*.*)'},...
        'Select first image');
    if(isa(fName,'char') & isa(dirName,'char'))
        % Recover all file names from the stack
        first_image_path=[dirName,filesep,fName];
        outFileList=getFileStackNames(first_image_path);        
        % Number of files 
        n=length(outFileList);
    else
        area=[];
        output_path = [];
        return
    end
else
    outFileList=getFileStackNames(first_image_path);        
    % Number of files 
    n=length(outFileList);
end

if nargin<4 || isempty(num_images_to_crop)
    % The user decides the number of images to be cropped
    prompt={'Specify the number of images to be cropped'};
    dlg_title='User input requested';
    num_lines=1;
    def={num2str(n)};
    num_images_to_crop=fix(str2num(char(inputdlg(prompt,dlg_title,num_lines,def))));

    % Check the selected number
    if isempty(num_images_to_crop)
        area=[];
        return
    end

    if num_images_to_crop>n
        fprintf(1,'Invalid number of images specified. Using the default value (%d).\n',num_images_to_crop);
    else
        % Crop outFileList
        n=num_images_to_crop;
        outFileList=outFileList(1:n);
    end
elseif num_images_to_crop==0
    outFileList=outFileList(1:n);
else
    if num_images_to_crop>n
        fprintf(1,'Invalid number of images specified. Using the default value (%d).\n',num_images_to_crop);
    else
        % Crop outFileList
        n=num_images_to_crop;
        outFileList=outFileList(1:n);
    end
end

% Read first image
imgOne=double(imread(first_image_path));

if isempty(area)
    h=figure;
    set(h,'NumberTitle','off');
    set(h,'Name','Please draw region of interest');
    % Normalize imgOne
    imgOne=(imgOne-min(imgOne(:)))/(max(imgOne(:))-min(imgOne(:)));
    % Crop - if the user closes the window without drawing, roipoly will return an error
    try
        [imgCropped,area]=imcrop(imgOne);
    catch
        uiwait(msgbox('No polygon selected. Quitting','Error','modal'));
        area=[];
        return
    end
    % Close figure
    close(h);
    % Check selected polygon
    if area(3)==0 | area(4)==0
        uiwait(msgbox('Please chose an area, not a single pixel.','Error','error','modal'));
        area=[];
        return
    end
    % Round area (imcrop can give also non-integer boundaries)
    area=round(area);
    % Vertices
    y0=area(2); y=area(2)+area(4);
    x0=area(1); x=area(1)+area(3);
else
    % Vertices
    y0=area(1); y=area(3);
    x0=area(2); x=area(4);
end

% Check boundaries
if y0<=0, y0=1; end
if x0<=0, x0=1; end
if y>=size(imgOne,1), y=size(imgOne,1); end
if x>=size(imgOne,2), x=size(imgOne,2); end

% Select output directory
if manualSave==1
    output_path=uigetdir('','Select output directory');
    if output_path==0
        area=[];
        return
    end
end

% Initializing waitbar
h=waitbar(0,'Processing...');

% Processing files
for i=1:n
    
    % Current filename
    currentFile=char(outFileList(i));
    [fpath,fName,fno,fext]=getFilenameBody(currentFile);
   
    % Read image from disk
    img=imread(currentFile);
   
    % Cut image
    imgC=img(y0:y,x0:x);
   
    % Prepare filename with path
    filename=[output_path,filesep,'crop_',fName,fno,fext];
   
    % Write file to disk
    if strcmp(fext,'.tif') || strcmp(fext,'.tiff') || strcmp(fext,'.TIF') || strcmp(fext,'.jpg') || strcmp(fext,'.jpeg')
        imwrite(imgC,filename, 'Compression', 'none');
    elseif strcmp(fext,'.bmp')
        imwrite(imgC,filename);     
    else
       disp('Unknown image format!!');
       return
    end
        
    % Update waitbar
    waitbar(i/n,h);
    
end

% Correct area
area=[y0 x0 y x];

% Close waitbar
close(h);
