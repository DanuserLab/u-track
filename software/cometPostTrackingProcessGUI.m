function varargout = cometPostTrackingProcessGUI(varargin)
% cometPostTrackingProcessGUI M-file for cometPostTrackingProcessGUI.fig
%      cometPostTrackingProcessGUI, by itself, creates a new cometPostTrackingProcessGUI or raises the existing
%      singleton*.
%
%      H = cometPostTrackingProcessGUI returns the handle to a new cometPostTrackingProcessGUI or the handle to
%      the existing singleton*.
%
%      cometPostTrackingProcessGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in cometPostTrackingProcessGUI.M with the given input arguments.
%
%      cometPostTrackingProcessGUI('Property','Value',...) creates a new cometPostTrackingProcessGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before cometPostTrackingProcessGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to cometPostTrackingProcessGUI_OpeningFcn via varargin.
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

% Edit the above text to modify the response to help cometPostTrackingProcessGUI

% Last Modified by GUIDE v2.5 16-Jul-2012 12:30:38

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @cometPostTrackingProcessGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @cometPostTrackingProcessGUI_OutputFcn, ...
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


% --- Executes just before cometPostTrackingProcessGUI is made visible.
function cometPostTrackingProcessGUI_OpeningFcn(hObject,eventdata,handles,varargin)

processGUI_OpeningFcn(hObject, eventdata, handles, varargin{:},'initChannel',1);

userData=get(handles.figure1,'UserData');
if(isempty(userData)), userData = struct(); end; % 2014b Hackaton-fix PR
funParams = userData.crtProc.funParams_;

set(handles.checkbox_makeHist,'Value',funParams.makeHist);
set(handles.checkbox_remBegEnd,'Value',funParams.remBegEnd);

% Initialize pop-up menus for reclassification schemes
set(handles.popupmenu_fgapReclassScheme,'String',...
    CometPostTrackingProcess.getFgapReclassificationSchemes, 'Value',...
    funParams.fgapReclassScheme);
set(handles.popupmenu_bgapReclassScheme,'String',...
    CometPostTrackingProcess.getBgapReclassificationSchemes, 'Value',...
    funParams.bgapReclassScheme);

% Choose default command line output for cometPostTrackingProcessGUI
handles.output = hObject;

% Update user data and GUI data
set(hObject, 'UserData', userData);
guidata(hObject, handles);

% --- Outputs from this function are returned to the command line.
function varargout = cometPostTrackingProcessGUI_OutputFcn(~, ~, handles) 
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
if(isempty(userData)), userData = struct(); end; % 2014b Hackaton-fix PR

if ishandle(userData.helpFig), delete(userData.helpFig); end

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

% Check user input
if isempty(get(handles.listbox_selectedChannels, 'String'))
    errordlg('Please select at least one input process from ''Available Movies''.','Setting Error','modal')
    return;
end
funParams.ChannelIndex = get(handles.listbox_selectedChannels,'UserData');

funParams.makeHist=get(handles.checkbox_makeHist,'Value');
funParams.remBegEnd=get(handles.checkbox_remBegEnd,'Value');

funParams.fgapReclassScheme=get(handles.popupmenu_fgapReclassScheme,'Value');
funParams.bgapReclassScheme=get(handles.popupmenu_bgapReclassScheme,'Value');

% Set parameters
processGUI_ApplyFcn(hObject, eventdata, handles,funParams);
