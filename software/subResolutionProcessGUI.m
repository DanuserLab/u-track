function varargout = subResolutionProcessGUI(varargin)
% SUBRESOLUTIONPROCESSGUI M-file for subResolutionProcessGUI.fig
%      SUBRESOLUTIONPROCESSGUI, by itself, creates a new SUBRESOLUTIONPROCESSGUI or raises the existing
%      singleton*.
%
%      H = SUBRESOLUTIONPROCESSGUI returns the handle to a new SUBRESOLUTIONPROCESSGUI or the handle to
%      the existing singleton*.
%
%      SUBRESOLUTIONPROCESSGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SUBRESOLUTIONPROCESSGUI.M with the given input arguments.
%
%      SUBRESOLUTIONPROCESSGUI('Property','Value',...) creates a new SUBRESOLUTIONPROCESSGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before subResolutionProcessGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to subResolutionProcessGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES
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

% Edit the above text to modify the response to help subResolutionProcessGUI

% Last Modified by GUIDE v2.5 16-Dec-2011 15:04:12

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @subResolutionProcessGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @subResolutionProcessGUI_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before subResolutionProcessGUI is made visible.
function subResolutionProcessGUI_OpeningFcn(hObject, eventdata, handles, varargin)

processGUI_OpeningFcn(hObject, eventdata, handles, varargin{:},'initChannel',1)

% Parameter Setup
userData = get(handles.figure1, 'UserData');
if(isempty(userData)), userData = struct(); end;

funParams = userData.crtProc.funParams_;

set(handles.edit_camBitDepth,'String',funParams.detectionParam.bitDepth);
set(handles.edit_psfSigma,'String',funParams.detectionParam.psfSigma)

arrayfun(@(x)set(handles.(['edit_alphaLocMax' num2str(x)]), 'String', funParams.detectionParam.alphaLocMax(x)), ...
                  1:length(funParams.detectionParam.alphaLocMax))
              
if  length(funParams.detectionParam.integWindow)>1 || funParams.detectionParam.integWindow
    set(handles.checkbox_rollingwindow, 'Value', 1)
    set(handles.text_windowSize, 'Enable', 'on')
    arrayfun(@(x)set(handles.(['edit_integWindow' num2str(x)]), 'Enable', 'on'), 1:3)
    arrayfun(@(x)set(handles.(['edit_integWindow' num2str(x)]), ...
        'String', num2str(funParams.detectionParam.integWindow(x)*2+1)), ...
                  1:length(funParams.detectionParam.integWindow))
end

% Mixture model fitting
set(handles.checkbox_mmf, 'Value', funParams.detectionParam.doMMF)
checkbox_mmf_Callback(hObject, eventdata, handles)

set(handles.edit_alphaR, 'String', num2str(funParams.detectionParam.testAlpha.alphaR))
set(handles.edit_alphaA, 'String', num2str(funParams.detectionParam.testAlpha.alphaA))
set(handles.edit_alphaD, 'String', num2str(funParams.detectionParam.testAlpha.alphaD))
set(handles.edit_alphaF, 'String', num2str(funParams.detectionParam.testAlpha.alphaF))

if funParams.detectionParam.numSigmaIter  
    set(handles.checkbox_iteration, 'Value', 1)
    set(handles.text_numSigmaIter, 'Enable', 'on')
    set(handles.edit_numSigmaIter, 'Enable', 'on', 'String',num2str(funParams.detectionParam.numSigmaIter) )
end

set(handles.checkbox_visual, 'Value', funParams.detectionParam.visual) 

if ~isempty(funParams.detectionParam.background)
    set(handles.checkbox_background, 'Value', 1)
    set(handles.edit_av_background, 'String', num2str(funParams.detectionParam.background.alphaLocMaxAbs))
    set(handles.edit_path, 'String', funParams.detectionParam.background.imageDir)  
else 
    set(handles.checkbox_background, 'Value', 0)
end
checkbox_background_Callback(hObject, eventdata, handles)

% funParams.movieParam
set(handles.edit_min, 'String', num2str(funParams.firstImageNum))
set(handles.edit_max, 'String', num2str(funParams.lastImageNum))
set(handles.text_framenum, 'String', ['(Totally ' num2str(userData.MD.nFrames_) ' frames in the movie)'])

% Choose default command line output for subResolutionProcessGUI
handles.output = hObject;

% Update user data and GUI data);
set(hObject, 'UserData', userData);
uicontrol(handles.pushbutton_done);
guidata(hObject, handles);


% --- Outputs from this function are returned to the command line.
function varargout = subResolutionProcessGUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;




% --- Executes on button press in checkbox_iteration.
function checkbox_mmf_Callback(hObject, eventdata, handles)

if get(handles.checkbox_mmf,'Value'), state='on'; else state='off'; end
set([handles.text_alphaR handles.edit_alphaR], 'Enable',state)


% --- Executes on button press in checkbox_iteration.
function checkbox_iteration_Callback(hObject, eventdata, handles)

if get(hObject, 'Value')
    state='on';
    % If no maximum number of iterations, use default number 10
    value = str2double(get(handles.edit_numSigmaIter, 'String'));
    if isnan(value), set(handles.edit_numSigmaIter, 'String', '10'); end
else
    state='off';
end

set([handles.text_numSigmaIter handles.edit_numSigmaIter], 'Enable',state)


% --- Executes on button press in checkbox_rollingwindow.
function checkbox_rollingwindow_Callback(hObject, eventdata, handles)

if get(hObject, 'Value')
    state='on';
    
    % If no input, use default number 1
    for i=1:3
        value = str2double(get(handles.(['edit_integWindow' num2str(i)]), 'String'));
        if ~isempty(get(handles.(['edit_alphaLocMax' num2str(i)]), 'String')) && isnan(value)          
            set(handles.(['edit_integWindow' num2str(i)]), 'String', '1')
        end
    end
    
else
    state='off';  
end
set([handles.text_windowSize handles.edit_integWindow1 handles.edit_integWindow2...
    handles.edit_integWindow3], 'Enable',state)


% --- Executes on button press in pushbutton_new.
function pushbutton_new_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_new (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[coor pathname] = cropStack([]);
if ~isempty(pathname) && ischar(pathname)
    set(handles.edit_path, 'String', pathname)
end


% --- Executes on button press in pushbutton_open.
function pushbutton_open_Callback(hObject, eventdata, handles)

pathname = uigetdir(pwd);
if isequal(pathname,0), return; end

set(handles.edit_path, 'String', pathname)

% --- Executes on button press in checkbox_background.
function checkbox_background_Callback(hObject, eventdata, handles)

if get(handles.checkbox_background, 'Value'), state= 'on'; else state='off'; end
set(get(handles.uipanel_background,'Children'),'Enable',state);

% --- Executes on button press in pushbutton_cancel.
function pushbutton_cancel_Callback(hObject, eventdata, handles)

userData = get(handles.figure1, 'UserData');
if(isempty(userData)), userData = struct(); end;

delete(handles.figure1);


% --- Executes on button press in pushbutton_done.
function pushbutton_done_Callback(hObject, eventdata, handles)
% Call back function of 'Apply' button

userData = get(handles.figure1, 'UserData');
if(isempty(userData)), userData = struct(); end;


if isempty(get(handles.listbox_selectedChannels, 'String'))
    errordlg('Please select at least one input channel from ''Available Channels''.','Setting Error','modal')
    return;
end
channelIndex = get (handles.listbox_selectedChannels, 'Userdata');
funParams.ChannelIndex = channelIndex;

% ------------------------- Check user input -----------------------------

missPara = [];

%% Required Fields

alphaLocMax{1} = get(handles.edit_alphaLocMax1, 'String');
alphaLocMax{2} = get(handles.edit_alphaLocMax2, 'String');
alphaLocMax{3} = get(handles.edit_alphaLocMax3, 'String');
integWindow{1} = get(handles.edit_integWindow1, 'String');
integWindow{2} = get(handles.edit_integWindow2, 'String');
integWindow{3} = get(handles.edit_integWindow3, 'String');

alphaR = get(handles.edit_alphaR, 'String');
alphaA = get(handles.edit_alphaA, 'String');
alphaD = get(handles.edit_alphaD, 'String');
alphaF = get(handles.edit_alphaF, 'String');

alphaLocMaxAbs = get(handles.edit_av_background, 'String');



%% Gaussian Standard Deviation: psfSigma
psfSigma = str2double(get(handles.edit_psfSigma, 'String'));
if isnan(psfSigma) || psfSigma < 0
    errordlg('Please provide a valid value to parameter "Gaussian Standard Deviation".','Error','modal')
    return
end

bitDepth = str2double(get(handles.edit_camBitDepth, 'String'));
if  isnan(bitDepth) || bitDepth < 0
    errordlg('Please provide a valid value to parameter "Camera Bit Depth".','Error','modal')
    return
end



%% Alpha-value for Local Maxima Detection: alphaLoc
validAlphaLocMax = ~cellfun(@isempty, alphaLocMax);
if all(~validAlphaLocMax)
    missPara = horzcat(missPara, sprintf('Alpha-value for Local Maxima Detection\n'));
    alphaLoc = 0.05; % default
    
elseif any(cellfun(@(x)isnan(str2double(x)), alphaLocMax(validAlphaLocMax))) || any(cellfun(@(x)(str2double(x) < 0), alphaLocMax(validAlphaLocMax)))
    errordlg('Please provide a valid value to parameter "Alpha-value for Local Maxima Detection".','Error','modal')
    return
    
else
    alphaLoc = cellfun(@(x)str2double(x), alphaLocMax(validAlphaLocMax), 'UniformOutput', true);
end

%% Window Size: integWin
temp = cellfun(@(x)~isempty(x), integWindow);
if get(handles.checkbox_rollingwindow, 'Value')
    if isempty(alphaLoc)
        missPara = horzcat(missPara, sprintf('Window Size\n'));
        integWin = (1-1)/2; % default
        
    elseif length(find(temp)) ~= length(alphaLoc)
        errordlg('The length of parameter "Camera Bit Depth" must be the same with parameter "Alpha-value for Local Maxima Detection".','Error','modal')
        return
        
    elseif any(cellfun(@(x)isnan(str2double(x)), integWindow(temp))) || any(cellfun(@(x)(str2double(x) < 0), integWindow(temp)))
        errordlg('Please provide a valid value to parameter "Camera Bit Depth".','Error','modal')   
        return
        
    else
        integWin = cellfun(@(x)str2double(x), integWindow(temp), 'UniformOutput', true);
        if ~all(mod(integWin, 2) == 1)
            errordlg('Parameter "Window Size" must be an odd number.','Error','modal')   
            return            
        end
        integWin = (integWin-1)/2;
    end
    
else
    integWin = 0;
end


%% AlphaR

if get(handles.checkbox_mmf, 'Value')
    if isempty( alphaR )
        missPara = horzcat(missPara, sprintf('Alpha Residual\n'));
        alphaR = 0.05; % default
        
    elseif isnan(str2double(alphaR)) || str2double(alphaR) < 0
        errordlg('Please provide a valid value to parameter "Alpha Residuals".','Error','modal')
        return   
        
    else
        alphaR = str2double(alphaR);
    end
else
    alphaR = .05;
end

%% AlphaA

if isempty( alphaA )
    missPara = horzcat(missPara, sprintf('Alpha Amplitude\n'));
    alphaA = 0.05; % default
        
elseif isnan(str2double(alphaA)) || str2double(alphaA) < 0
    errordlg('Please provide a valid value to parameter "Alpha Amplitude".','Error','modal')
    return   
        
else
    alphaA = str2double(alphaA);
end

%% AlphaD

if isempty( alphaD )
    missPara = horzcat(missPara, sprintf('Alpha Distance\n'));
    alphaD = 0.05; % default 
        
elseif isnan(str2double(alphaD)) || str2double(alphaD) < 0
    errordlg('Please provide a valid value to parameter "Alpha Distance".','Error','modal')
    return   
        
else
    alphaD = str2double(alphaD);
end

%% AlphaF

if isempty( alphaF )
    missPara = horzcat(missPara, sprintf('Alpha Final\n'));
    alphaF = 0; % default
        
elseif isnan(str2double(alphaF)) || str2double(alphaF) < 0
    errordlg('Please provide a valid value to parameter "Alpha Final".','Error','modal')
    return   
        
else
    alphaF = str2double(alphaF);
end

%% Maximum Number of Interations: numSigmaIter

if get(handles.checkbox_iteration, 'Value')
    numSigmaIter = str2double(get(handles.edit_numSigmaIter, 'String'));
    if isnan(numSigmaIter) || numSigmaIter < 0 || floor(numSigmaIter) ~= ceil(numSigmaIter)
        errordlg('Please provide a valid value to parameter "Maximum Number of Interations".','Error','modal')
        return   
    end
else
    numSigmaIter = 0;
end

%% Alpha-value for Local Maxima Detection

if get(handles.checkbox_background, 'Value')
    if isempty( alphaLocMaxAbs )
        missPara = horzcat(missPara, sprintf('Background Alpha-value \n'));
        alphaLocMaxAbs = 0.001; % default 

    elseif isnan(str2double(alphaLocMaxAbs)) || str2double(alphaLocMaxAbs) < 0 
        errordlg('Please provide a valid value to parameter "Background Alpha-value".','Error','modal')
        return   

    else
        alphaLocMaxAbs = str2double(alphaLocMaxAbs);
    end    
    
else
    alphaLocMaxAbs = [];
end

%% Background directory: bg_dir

if get(handles.checkbox_background, 'Value')
    bg_dir = get(handles.edit_path, 'String');
    if isempty( bg_dir )
        errordlg('Please specify a background image directory.','Error','modal')
        return
    end    
    if ~strcmp(bg_dir(end), filesep), bg_dir = [bg_dir filesep]; end
end

%% Frame Index: minid, maxid

minid = str2double(get(handles.edit_min, 'String'));
if isnan(minid) || minid <= 0 || floor(minid) ~= ceil(minid)
    errordlg('Please provide a valid value to parameter "Minimum Frame Index".','Error','modal')
    return
end

maxid = str2double(get(handles.edit_max, 'String'));
if isnan(maxid) || maxid <= 0 || floor(maxid) ~= ceil(maxid)
    errordlg('Please provide a valid value to parameter "Maximum Frame Index".','Error','modal')
    return
end

if minid > maxid
    errordlg('Minimum frame index is larger than maximum frame index.','Error','modal')
    return
elseif maxid > userData.MD.nFrames_
    errordlg('Frame index exceeds the number of frames.', 'Error', 'modal')
    return
end

%% TO-DO check validation of output dir and background image dir

% Check background image directory
if get(handles.checkbox_background, 'Value')
    
    fileNames = imDir(bg_dir);
    if isempty(fileNames)
        errordlg(sprintf('No valid image file found in the background image directory: %s.', bg_dir),'Error','modal')
        return           
    end

    [x1 filenameBase x3 x4] =  getFilenameBody(fileNames(1).name);
end


funParams.firstImageNum = minid;
funParams.lastImageNum = maxid;

%% funParams.detectionParam
funParams.detectionParam.psfSigma = psfSigma;
funParams.detectionParam.bitDepth = bitDepth;
funParams.detectionParam.alphaLocMax = alphaLoc;
funParams.detectionParam.integWindow = integWin;
funParams.detectionParam.doMMF = get(handles.checkbox_mmf, 'Value');
funParams.detectionParam.testAlpha.alphaR = alphaR;
funParams.detectionParam.testAlpha.alphaA = alphaA;
funParams.detectionParam.testAlpha.alphaD = alphaD;
funParams.detectionParam.testAlpha.alphaF = alphaF;
funParams.detectionParam.numSigmaIter = numSigmaIter;

funParams.detectionParam.visual =  get(handles.checkbox_visual, 'Value');

if get(handles.checkbox_background, 'Value')
    funParams.detectionParam.background.imageDir = bg_dir;
    funParams.detectionParam.background.alphaLocMaxAbs = alphaLocMaxAbs;
    funParams.detectionParam.background.filenameBase = filenameBase;
else
    funParams.detectionParam.background = [];
end

% Apply parameters
setLastImageNum =@(x) parseProcessParams(x,struct('lastImageNum',...
    min(x.funParams_.lastImageNum,x.owner_.nFrames_)));
processGUI_ApplyFcn(hObject,eventdata,handles,funParams,{setLastImageNum});


% --- Executes on key press with focus on figure1 and none of its controls.
function figure1_KeyPressFcn(hObject, eventdata, handles)

if strcmp(eventdata.Key, 'return')
    pushbutton_done_Callback(handles.pushbutton_done, [], handles);
end
