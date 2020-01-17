function varargout = plusTipCostMatCloseGapsGUI(varargin)
% PLUSTIPCOSTMATCLOSEGAPSGUI M-file for plusTipCostMatCloseGapsGUI.fig
%      PLUSTIPCOSTMATCLOSEGAPSGUI, by itself, creates a new PLUSTIPCOSTMATCLOSEGAPSGUI or raises the existing
%      singleton*.
%
%      H = PLUSTIPCOSTMATCLOSEGAPSGUI returns the handle to a new PLUSTIPCOSTMATCLOSEGAPSGUI or the handle to
%      the existing singleton*.
%
%      PLUSTIPCOSTMATCLOSEGAPSGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PLUSTIPCOSTMATCLOSEGAPSGUI.M with the given input arguments.
%
%      PLUSTIPCOSTMATCLOSEGAPSGUI('Property','Value',...) creates a new PLUSTIPCOSTMATCLOSEGAPSGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before plusTipCostMatCloseGapsGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to plusTipCostMatCloseGapsGUI_OpeningFcn via varargin.
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

% Edit the above text to modify the response to help plusTipCostMatCloseGapsGUI

% Last Modified by GUIDE v2.5 09-Feb-2012 16:27:04

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @plusTipCostMatCloseGapsGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @plusTipCostMatCloseGapsGUI_OutputFcn, ...
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


% --- Executes just before plusTipCostMatCloseGapsGUI is made visible.
function plusTipCostMatCloseGapsGUI_OpeningFcn(hObject, eventdata, handles, varargin)

costMat_OpeningFcn(hObject, eventdata, handles, varargin{:})
userData = get(handles.figure1, 'UserData');
if(isempty(userData)), userData = struct(); end; % 2014b Hackaton-fix PR

parameters = userData.parameters;

% Brownian motion parameters
userData.numParams = {'maxFAngle','maxBAngle','backVelMultFactor','fluctRad'};
set(handles.checkbox_breakNonLinearTracks, 'Value', parameters.breakNonLinearTracks)
for i=1:numel(userData.numParams)
    set(handles.(['edit_' userData.numParams{i}]), 'String', parameters.(userData.numParams{i}));
end

% Update handles structure
set(handles.figure1, 'UserData', userData);
handles.output = hObject;
guidata(hObject, handles);

% UIWAIT makes plusTipCostMatCloseGapsGUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = plusTipCostMatCloseGapsGUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton_cancel.
function pushbutton_cancel_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_cancel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
delete(handles.figure1)

% --- Executes on button press in pushbutton_done.
function pushbutton_done_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_done (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

userData = get(handles.figure1, 'UserData');
if(isempty(userData)), userData = struct(); end; % 2014b Hackaton-fix PR
parameters = userData.parameters;

% Brownian motion parameters
parameters.breakNonLinearTracks = get(handles.checkbox_breakNonLinearTracks, 'Value');
isPosScalar = @(x) isscalar(x) &&~isnan(x) && x>=0;   
for i=1:numel(userData.numParams)
    value = str2double(get(handles.(['edit_' userData.numParams{i}]), 'String'));
    if ~isPosScalar(value)
        errordlg(['Please provide a valid value to parameter '...
            get(handles.(['text_' userData.numParams{i}]), 'String') '.'],'Error','modal')
        return
    end
    parameters.(userData.numParams{i})=value;
end

u = get(userData.handles_main.popupmenu_gapclosing, 'UserData');
u{userData.procID} = parameters;

set(userData.handles_main.popupmenu_gapclosing, 'UserData', u)


set(handles.figure1, 'UserData', userData);
guidata(hObject,handles);
delete(handles.figure1);

% --- Executes on key press with focus on figure1 and none of its controls.
function figure1_KeyPressFcn(hObject, eventdata, handles)

if strcmp(eventdata.Key, 'return')
    pushbutton_done_Callback(handles.pushbutton_done, [], handles);
end
