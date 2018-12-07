function varargout = motionAnalysisProcessGUI(varargin)
% motionAnalysisProcessGUI M-file for motionAnalysisProcessGUI.fig
%      motionAnalysisProcessGUI, by itself, creates a new motionAnalysisProcessGUI or raises the existing
%      singleton*.
%
%      H = motionAnalysisProcessGUI returns the handle to a new motionAnalysisProcessGUI or the handle to
%      the existing singleton*.
%
%      motionAnalysisProcessGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in motionAnalysisProcessGUI.M with the given input arguments.
%
%      motionAnalysisProcessGUI('Property','Value',...) creates a new motionAnalysisProcessGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before motionAnalysisProcessGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to motionAnalysisProcessGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES
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

% Edit the above text to modify the response to help motionAnalysisProcessGUI

% Last Modified by GUIDE v2.5 09-Jan-2017 11:37:47

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @motionAnalysisProcessGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @motionAnalysisProcessGUI_OutputFcn, ...
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


% --- Executes just before motionAnalysisProcessGUI is made visible.
function motionAnalysisProcessGUI_OpeningFcn(hObject, eventdata, handles, varargin)

processGUI_OpeningFcn(hObject, eventdata, handles, varargin{:},'initChannel',1);

% Set default parameters
userData = get(handles.figure1, 'UserData');
if(isempty(userData)), userData = struct(); end;
funParams = userData.crtProc.funParams_;

set(handles.popupmenu_probDim,'String',{'2','3'},'UserData',[2 3],...
    'Value',find(funParams.probDim==[2 3]));

set(handles.checkbox_checkAsym,'Value',funParams.checkAsym);

% Set confinement radius methods
confRadMethods = MotionAnalysisProcess.getConfinementRadiusMethods;
set(handles.popupmenu_confRadMin,'String',{confRadMethods.name},...
    'Value',find(funParams.confRadMin==[confRadMethods.type]),...
    'UserData',[confRadMethods.type]);

% Set alpha values
alphaValues = MotionAnalysisProcess.getAlphaValues;

% Set to use new MSS if MSS alpha value is negative
set(handles.checkbox_NewMSSThresh,'Value',funParams.alphaValues(1) < 0);

% Display only positive alpha values
set(handles.popupmenu_alphaValueMSS,'String',num2cell(alphaValues),...
    'Value',find(abs(funParams.alphaValues(1))==alphaValues),...
    'UserData',alphaValues);
set(handles.popupmenu_alphaValueAsym,'String',num2cell(alphaValues),...
    'Value',find(funParams.alphaValues(2)==alphaValues),...
    'UserData',alphaValues);

% Update GUI user data
handles.output = hObject;
set(hObject, 'UserData', userData);
guidata(hObject, handles);

% --- Outputs from this function are returned to the command line.
function varargout = motionAnalysisProcessGUI_OutputFcn(~, ~, handles) 
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

delete(userData.helpFig(ishandle(userData.helpFig))); 

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
    errordlg('Please select at least one input channel from ''Available Channels''.','Setting Error','modal')
    return;
end
funParams.ChannelIndex = get(handles.listbox_selectedChannels, 'Userdata');

% Get frame range
props = get(handles.popupmenu_probDim, {'UserData','Value'});
funParams.probDim=props{1}(props{2});

funParams.checkAsym=get(handles.checkbox_checkAsym,'Value');

% Get alpha values
props = get(handles.popupmenu_alphaValueMSS, {'UserData','Value'});
funParams.alphaValues(1) =props{1}(props{2});
if(get(handles.checkbox_NewMSSThresh,'Value'))
    funParams.alphaValues(1) = -funParams.alphaValues(1);
end
props = get(handles.popupmenu_alphaValueAsym, {'UserData','Value'});
funParams.alphaValues(2) =props{1}(props{2});

props = get(handles.popupmenu_confRadMin, {'UserData','Value'});
funParams.confRadMin=props{1}(props{2});


processGUI_ApplyFcn(hObject, eventdata, handles,funParams);


% --- Executes on button press in checkbox_NewMSSThresh.
function checkbox_NewMSSThresh_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_NewMSSThresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_NewMSSThresh
