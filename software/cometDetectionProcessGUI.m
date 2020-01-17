function varargout = cometDetectionProcessGUI(varargin)
% cometDetectionProcessGUI M-file for cometDetectionProcessGUI.fig
%      cometDetectionProcessGUI, by itself, creates a new cometDetectionProcessGUI or raises the existing
%      singleton*.
%
%      H = cometDetectionProcessGUI returns the handle to a new cometDetectionProcessGUI or the handle to
%      the existing singleton*.
%
%      cometDetectionProcessGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in cometDetectionProcessGUI.M with the given input arguments.
%
%      cometDetectionProcessGUI('Property','Value',...) creates a new cometDetectionProcessGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before cometDetectionProcessGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to cometDetectionProcessGUI_OpeningFcn via varargin.
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

% Edit the above text to modify the response to help cometDetectionProcessGUI

% Last Modified by GUIDE v2.5 09-Feb-2012 18:03:09

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @cometDetectionProcessGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @cometDetectionProcessGUI_OutputFcn, ...
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


% --- Executes just before cometDetectionProcessGUI is made visible.
function cometDetectionProcessGUI_OpeningFcn(hObject, eventdata, handles, varargin)

processGUI_OpeningFcn(hObject, eventdata, handles, varargin{:},'initChannel',1);

% ---------------------- Channel Setup -------------------------
userData = get(handles.figure1, 'UserData');
if(isempty(userData)), userData = struct(); end; % 2014b Hackaton-fix PR
funParams = userData.crtProc.funParams_;

maskProc =  cellfun(@(x) isa(x,'MaskProcess'),userData.MD.processes_);
maskProcID=find(maskProc);
maskProcNames = cellfun(@(x) x.getName(),userData.MD.processes_(maskProc),'Unif',false);
maskProcString = vertcat('None',maskProcNames(:));
maskProcData=horzcat({[]},num2cell(maskProcID));
maskProcValue = find(cellfun(@(x) isequal(x,funParams.MaskProcessIndex),maskProcData));
if isempty(maskProcValue), maskProcValue = 1; end
set(handles.popupmenu_MaskProcess,'String',maskProcString,...
    'UserData',maskProcData,'Value',maskProcValue,'Enable','on');


set(handles.edit_firstFrame,'String',funParams.firstFrame);
set(handles.edit_lastFrame,'String',funParams.lastFrame);
userData.numParams ={'sigma1','sigma2','multFactorThresh','multFactorStepSize'};
cellfun(@(x) set(handles.(['edit_' x]),'String',funParams.(x)),...
    userData.numParams)
% Update GUI user data
handles.output = hObject;
set(hObject, 'UserData', userData);
guidata(hObject, handles);

% --- Outputs from this function are returned to the command line.
function varargout = cometDetectionProcessGUI_OutputFcn(~, ~, handles) 
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

if isfield(userData, 'helpFig') && ishandle(userData.helpFig)
   delete(userData.helpFig) 
end

if isfield(userData, 'previewFig') && ishandle(userData.previewFig)
   delete(userData.previewFig) 
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

% Check user input
if isempty(get(handles.listbox_selectedChannels, 'String'))
    errordlg('Please select at least one input channel from ''Available Channels''.','Setting Error','modal')
    return;
end
funParams.ChannelIndex = get(handles.listbox_selectedChannels, 'Userdata');

% Get frame range
userData=get(handles.figure1,'UserData');
if(isempty(userData)), userData = struct(); end; % 2014b Hackaton-fix PR
firstFrame = str2double(get(handles.edit_firstFrame,'String'));
if ~userData.crtProc.checkFrameNum(firstFrame)
    errordlg(['Please enter a valid ' get(handles.text_frameRange,'String') '.'],...
        'Setting Error','modal')
    return;
end
funParams.firstFrame=firstFrame;

% Get last range
lastFrame = str2double(get(handles.edit_lastFrame,'String'));
if ~userData.crtProc.checkFrameNum(lastFrame) || firstFrame > lastFrame
    errordlg(['Please enter a valid ' get(handles.text_frameRange,'String') '.'],...
        'Setting Error','modal')
    return;
end
funParams.lastFrame=lastFrame;

% Retrieve detection parameters
for i=1:numel(userData.numParams)  
    value = str2double(get(handles.(['edit_' userData.numParams{i}]),'String'));
    if isnan(value) || value < 0
        errordlg(['Please enter a valid value for '...
            get(handles.(['text_' userData.numParams{i}]),'String') '.'],...
            'Setting Error','modal')
        return;
    end
    funParams.(userData.numParams{i})=value; 
end

% Retrieve numerical parameters
for i=1:numel(userData.numParams)  
    value = str2double(get(handles.(['edit_' userData.numParams{i}]),'String'));
    if isnan(value) || value < 0
        errordlg(['Please enter a valid value for '...
            get(handles.(['text_' userData.numParams{i}]),'String') '.'],...
            'Setting Error','modal')
        return;
    end
    funParams.(userData.numParams{i})=value; 
end

% Retrieve mask process index and class (for propagation)
props=get(handles.popupmenu_MaskProcess,{'UserData','Value'});
funParams.MaskProcessIndex = props{1}{props{2}};
if ~isempty(funParams.MaskProcessIndex)
    maskProcessClass=class(userData.MD.processes_{funParams.MaskProcessIndex});
else
    maskProcessClass = '';
end

% Set parameters
setMaskProcess = @(x) parseProcessParams(x, struct('MaskProcessIndex',...
    x.owner_.getProcessIndex(maskProcessClass,1,false)));
settingFcn = {setMaskProcess};

% Handle multiple movies last frame setting
if  lastFrame == userData.MD.nFrames_
    setLastFrame =@(x) parseProcessParams(x,struct('lastFrame',...
        min(x.owner_.nFrames_)));
    settingFcn = {setMaskProcess, setLastFrame};
end

processGUI_ApplyFcn(hObject, eventdata, handles,funParams,settingFcn);
