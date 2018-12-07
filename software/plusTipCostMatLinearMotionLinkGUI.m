function varargout = plusTipCostMatLinearMotionLinkGUI(varargin)
% PLUSTIPCOSTMATLINEARMOTIONLINKGUI M-file for plusTipCostMatLinearMotionLinkGUI.fig
%      PLUSTIPCOSTMATLINEARMOTIONLINKGUI, by itself, creates a new PLUSTIPCOSTMATLINEARMOTIONLINKGUI or raises the existing
%      singleton*.
%
%      H = PLUSTIPCOSTMATLINEARMOTIONLINKGUI returns the handle to a new PLUSTIPCOSTMATLINEARMOTIONLINKGUI or the handle to
%      the existing singleton*.
%
%      PLUSTIPCOSTMATLINEARMOTIONLINKGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PLUSTIPCOSTMATLINEARMOTIONLINKGUI.M with the given input arguments.
%
%      PLUSTIPCOSTMATLINEARMOTIONLINKGUI('Property','Value',...) creates a new PLUSTIPCOSTMATLINEARMOTIONLINKGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before plusTipCostMatLinearMotionLinkGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to plusTipCostMatLinearMotionLinkGUI_OpeningFcn via varargin.
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

% Edit the above text to modify the response to help plusTipCostMatLinearMotionLinkGUI

% Last Modified by GUIDE v2.5 14-Feb-2013 17:00:53

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @plusTipCostMatLinearMotionLinkGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @plusTipCostMatLinearMotionLinkGUI_OutputFcn, ...
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


% --- Executes just before plusTipCostMatLinearMotionLinkGUI is made visible.
function plusTipCostMatLinearMotionLinkGUI_OpeningFcn(hObject, eventdata, handles, varargin)

costMat_OpeningFcn(hObject, eventdata, handles, varargin{:})
userData = get(handles.figure1, 'UserData');
if(isempty(userData)), userData = struct(); end; % 2014b Hackaton-fix PR
parameters = userData.parameters;

% Parameter Setup
set(handles.edit_lower, 'String', num2str(parameters.minSearchRadius))
set(handles.edit_upper, 'String', num2str(parameters.maxSearchRadius))
set(handles.edit_brownStdMult, 'String', num2str(parameters.brownStdMult))
set(handles.checkbox_useLocalDensity, 'Value', parameters.useLocalDensity)
set(handles.edit_nnWindow, 'String', num2str(parameters.nnWindow))

if isempty(parameters.diagnostics) || (length(parameters.diagnostics) == 1 && parameters.diagnostics == 0)
    set(handles.checkbox_diagnostics, 'Value', 0);
    set(get(handles.uipanel_diagnostics,'Children'),'Enable','off');
else
    set(handles.checkbox_diagnostics, 'Value', 1);
    for i = 1:min(3,length(parameters.diagnostics))
        set(handles.(['edit_diag_' num2str(i)]),...
            'String',num2str(parameters.diagnostics(i)));
    end
end

% Update handles structure
handles.output = hObject;
guidata(hObject, handles);

% UIWAIT makes plusTipCostMatLinearMotionLinkGUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = plusTipCostMatLinearMotionLinkGUI_OutputFcn(hObject, eventdata, handles) 
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

userData = get(handles.figure1, 'UserData');
if(isempty(userData)), userData = struct(); end; % 2014b Hackaton-fix PR
lower = str2double(get(handles.edit_lower, 'String'));
upper = str2double(get(handles.edit_upper, 'String'));
brownStdMult = str2double(get(handles.edit_brownStdMult, 'String'));
nnWindow = str2double(get(handles.edit_nnWindow, 'String'));
diagnostics{1} = get(handles.edit_diag_1, 'String');
diagnostics{2} = get(handles.edit_diag_2, 'String');
diagnostics{3} = get(handles.edit_diag_3, 'String');

% lower
isPosScalar = @(x) isscalar(x) && ~isnan(x) && x>=0;
if ~isPosScalar(lower)
    errordlg('Please provide a valid value to parameter "Lower Bound".','Error','modal')
    return
end

% Upper
if ~isPosScalar(upper)
    errordlg('Please provide a valid value to parameter "Upper Bound".','Error','modal')
    return
elseif upper < lower
    errordlg('"Upper Bound" should be larger than "Lower Bound".','Error','modal')
    return
end

% brownStdMult
if ~isPosScalar(brownStdMult)
    errordlg('Please provide a valid value to parameter "Multiplication Factor for Search Radius Calculation".','Error','modal')
    return
end

% nnWindow
if ~isPosScalar(nnWindow)
    errordlg('Please provide a valid value to parameter "Number of Frames for Nearest Neighbor Distance Calculation".','Error','modal')
    return
end

% Set Parameters
parameters = userData.parameters;

parameters.minSearchRadius = lower;
parameters.maxSearchRadius = upper;
parameters.brownStdMult = brownStdMult;
parameters.useLocalDensity = get(handles.checkbox_useLocalDensity, 'Value');
parameters.nnWindow = nnWindow;


% Set diagnostics parameters
if get(handles.checkbox_diagnostics, 'Value')
    
    if all(cellfun(@isempty, diagnostics))
        errordlg('Please provide 1 or more than 1 (maximum 3) "Frame Numbers to Plot Histograms" in the text boxes.','Error','modal')
        return
    end
    
    validDiagnostics = str2double(diagnostics(~cellfun(@isempty, diagnostics)));
    nFrames= userData.crtProc.owner_.nFrames_;
    if any(isnan(validDiagnostics) | validDiagnostics<2 | validDiagnostics>nFrames-1)
        errordlg('Please provide a valid value to parameter "Frame Numbers to Plot Histograms". Note: the first or last frame of a movie is invalid.','Error','modal')
        return
    end
    parameters.diagnostics =validDiagnostics;
else
    parameters.diagnostics = [];      
end


u = get(userData.handles_main.popupmenu_linking, 'UserData');
u{userData.procID} = parameters;

set(userData.handles_main.popupmenu_linking, 'UserData', u)


set(handles.figure1, 'UserData', userData);
guidata(hObject,handles);
delete(handles.figure1);

% --- Executes on button press in checkbox_diagnostics.
function checkbox_diagnostics_Callback(hObject, eventdata, handles)

if get(hObject, 'Value'),
    set(get(handles.uipanel_diagnostics,'Children'),'Enable','on');
else
    set(get(handles.uipanel_diagnostics,'Children'),'Enable','off');
end

% --- Executes on key press with focus on figure1 and none of its controls.
function figure1_KeyPressFcn(hObject, eventdata, handles)

if strcmp(eventdata.Key, 'return')
    pushbutton_done_Callback(handles.pushbutton_done, [], handles);
end
