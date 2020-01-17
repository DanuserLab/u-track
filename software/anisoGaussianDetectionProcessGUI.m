function varargout = anisoGaussianDetectionProcessGUI(varargin)
% anisoGaussianDetectionProcessGUI M-file for anisoGaussianDetectionProcessGUI.fig
%      anisoGaussianDetectionProcessGUI, by itself, creates a new anisoGaussianDetectionProcessGUI or raises the existing
%      singleton*.
%
%      H = anisoGaussianDetectionProcessGUI returns the handle to a new anisoGaussianDetectionProcessGUI or the handle to
%      the existing singleton*.
%
%      anisoGaussianDetectionProcessGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in anisoGaussianDetectionProcessGUI.M with the given input arguments.
%
%      anisoGaussianDetectionProcessGUI('Property','Value',...) creates a new anisoGaussianDetectionProcessGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before anisoGaussianDetectionProcessGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to anisoGaussianDetectionProcessGUI_OpeningFcn via varargin.
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

% Edit the above text to modify the response to help anisoGaussianDetectionProcessGUI

% Last Modified by GUIDE v2.5 18-Feb-2013 12:23:43

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @anisoGaussianDetectionProcessGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @anisoGaussianDetectionProcessGUI_OutputFcn, ...
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


% --- Executes just before anisoGaussianDetectionProcessGUI is made visible.
function anisoGaussianDetectionProcessGUI_OpeningFcn(hObject, eventdata, handles, varargin)

processGUI_OpeningFcn(hObject, eventdata, handles, varargin{:},'initChannel',1);

% Set-up parameters
userData=get(handles.figure1,'UserData');
funParams = userData.crtProc.funParams_;

% Set-up parameters
userData.numParams = {'psfSigma', 'alpha', 'kSigma', 'minDist'};
for i =1 : numel(userData.numParams)
    paramName = userData.numParams{i};
    set(handles.(['edit_' paramName]), 'String', funParams.(paramName));
end

% Update GUI user data
set(handles.figure1, 'UserData', userData);
handles.output = hObject;
guidata(hObject, handles);

% --- Outputs from this function are returned to the command line.
function varargout = anisoGaussianDetectionProcessGUI_OutputFcn(~, ~, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes on button press in pushbutton_cancel.
function pushbutton_cancel_Callback(~, ~, handles)
% Delete figure
delete(handles.figure1);

% --- Executes during object deletion, before destroying properties.
function figure1_DeleteFcn(hObject, ~, handles)
% Notify the package GUI that the setting panel is closed
userData = get(handles.figure1, 'UserData');
if(isempty(userData)), userData = struct(); end;

if isfield(userData, 'helpFig') && ishandle(userData.helpFig)
   delete(userData.helpFig) 
end

set(handles.figure1, 'UserData', userData);
guidata(hObject,handles);


% --- Executes on key press with focus on pushbutton_done and none of its controls.
function pushbutton_done_KeyPressFcn(~, eventdata, handles)

if strcmp(eventdata.Key, 'return')
    pushbutton_done_Callback(handles.pushbutton_done, [], handles);
end

% --- Executes on key press with focus on figure1 and none of its controls.
function figure1_KeyPressFcn(~, eventdata, handles)

if strcmp(eventdata.Key, 'return')
    pushbutton_done_Callback(handles.pushbutton_done, [], handles);
end

% --- Executes on button press in pushbutton_done.
function pushbutton_done_Callback(hObject, eventdata, handles)

% -------- Check user input --------

if isempty(get(handles.listbox_selectedChannels, 'String'))
    errordlg('Please select at least one input channel from ''Available Channels''.','Setting Error','modal')
    return;
end

% Retrieve GUI-defined parameters
channelIndex = get(handles.listbox_selectedChannels, 'Userdata');
funParams.ChannelIndex = channelIndex;

% Retrieve detection parameters
userData = get(handles.figure1, 'UserData');
if(isempty(userData)), userData = struct(); end;
for i = 1:numel(userData.numParams)
    paramName = userData.numParams{i};
    value = str2double(get(handles.(['edit_' paramName]),'String'));
    if isnan(value) || value < 0
        errordlg(['Please enter a valid value for '...
            get(handles.(['text_' paramName]),'String') '.'],...
            'Setting Error','modal')
        return;
    end
    funParams.(paramName)=value; 
end

% Add 64-bit warning
is64bit = ~isempty(regexp(computer ,'64$', 'once'));
if ~is64bit
    warndlg(['Your Matlab version is not detected as 64-bit. Please note '....
        'the anisotropic Gaussian detection uses compiled MEX files which '...
        'are not provided for 32-bit.'],...
        'Setting Error','modal');
end

processGUI_ApplyFcn(hObject, eventdata, handles,funParams);
