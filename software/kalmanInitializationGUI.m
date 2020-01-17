function varargout = kalmanInitializationGUI(varargin)
% KALMANINITIALIZATIONGUI M-file for kalmanInitializationGUI.fig
%      KALMANINITIALIZATIONGUI, by itself, creates a new KALMANINITIALIZATIONGUI or raises the existing
%      singleton*.
%
%      H = KALMANINITIALIZATIONGUI returns the handle to a new KALMANINITIALIZATIONGUI or the handle to
%      the existing singleton*.
%
%      KALMANINITIALIZATIONGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in KALMANINITIALIZATIONGUI.M with the given input arguments.
%
%      KALMANINITIALIZATIONGUI('Property','Value',...) creates a new KALMANINITIALIZATIONGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before kalmanInitializationGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to kalmanInitializationGUI_OpeningFcn via varargin.
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

% Edit the above text to modify the response to help kalmanInitializationGUI

% Last Modified by GUIDE v2.5 13-Dec-2011 17:07:36

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @kalmanInitializationGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @kalmanInitializationGUI_OutputFcn, ...
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


% --- Executes just before kalmanInitializationGUI is made visible.
function kalmanInitializationGUI_OpeningFcn(hObject, eventdata, handles, varargin)


costMat_OpeningFcn(hObject, eventdata, handles, varargin{:})
userData = get(handles.figure1, 'UserData');
if(isempty(userData)), userData = struct(); end;
parameters = userData.parameters;

% Get main figure handle and process id
props = get(userData.handles_main.popupmenu_probDim, {'UserData','Value'});
userData.probDim=props{1}(props{2});

% Parameter Setup
set(handles.radiobutton_none,'Value',1);
if ~isempty(parameters)
    if isfield(parameters, 'initVelocity') && ~isempty(parameters.initVelocity)% Initial Valocity Estimate
        for i=1:userData.probDim  
            set(handles.(['edit_v_' num2str(i)]), 'String', parameters.initVelocity(i));
        end
        
        set(handles.radiobutton_initVelocity, 'Value', 1); 
    end
        
    if isfield(parameters, 'convergePoint') && ~isempty(parameters.convergePoint) % Reference Point for Initial Estimate
        for i=1:userData.probDim  
            set(handles.(['edit_' num2str(i)]), 'String', parameters.convergePoint(i));
        end
        
        set(handles.radiobutton_convergePoint, 'Value', 1);         
    end
    
    set(handles.edit_radius, 'String', num2str(parameters.searchRadiusFirstIteration))
end



handles.output = hObject;
set(handles.figure1, 'UserData', userData)
uipanel5_SelectionChangeFcn(hObject,eventdata,handles);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes kalmanInitializationGUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = kalmanInitializationGUI_OutputFcn(hObject, eventdata, handles) 
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
if(isempty(userData)), userData = struct(); end;
parameters = userData.parameters;

initVelFlag=get(handles.radiobutton_initVelocity, 'Value');
dimensions=1:userData.probDim;
if initVelFlag
    initVelocity=arrayfun(@(x) str2double(get(handles.(['edit_v_' num2str(x)]),'String')),dimensions);
    if any(isnan(initVelocity)) || any(initVelocity<0)
        errordlg(['Please provide a valid value to parameter' ...
            get(handles.radiobutton_initVelocity,'String') '.'],'Error','modal')
        return  
    end
else
    initVelocity=[];
end

convPointFlag=get(handles.radiobutton_convergePoint, 'Value');
if convPointFlag
    convergePoint=arrayfun(@(x) str2double(get(handles.(['edit_' num2str(x)]),'String')),dimensions);
    if any(isnan(convergePoint)) || any(convergePoint<0)
        errordlg(['Please provide a valid value to parameter' ...
            get(handles.radiobutton_convergePoint,'String') '.'],'Error','modal')
        return  
    end
else
    convergePoint=[];
end 

searchRadiusFlag=~isempty(get(handles.edit_radius, 'String'));
if searchRadiusFlag
    searchRadiusFirstIteration=str2double(get(handles.edit_radius, 'String'));
    if isnan(searchRadiusFirstIteration) || searchRadiusFirstIteration <0
        errordlg('Please provide a valid value to parameter "Search Radius for Iteration".','Error','modal')
        return
    end
else
    searchRadiusFirstIteration=[];
end

if ~initVelFlag && ~convPointFlag && ~searchRadiusFlag
    parameters = [];
else
    parameters.initVelocity = initVelocity;
    parameters.convergePoint = convergePoint;
    parameters.searchRadiusFirstIteration = searchRadiusFirstIteration;
end

u = get(userData.handles_main.popupmenu_kalmanFunctions, 'UserData');
u{userData.procID} = parameters;

set(userData.handles_main.popupmenu_kalmanFunctions, 'UserData', u)   

set(handles.figure1, 'UserData', userData);
guidata(hObject,handles);
delete(handles.figure1);


% --- Executes when selected object is changed in uipanel5.
function uipanel5_SelectionChangeFcn(hObject, eventdata, handles)
handles = guidata(hObject); 

userData = get(handles.figure1, 'UserData');
if(isempty(userData)), userData = struct(); end;
% Highlight the content under new radiobuttonfunction uipanel5_SelectionChangeFcn(hObject, eventdata, handles)
selectedButton = get(get(handles.uipanel5,'SelectedObject'),'Tag');
if strcmpi(selectedButton,'radiobutton_initVelocity');
    child=get(handles.uipanel_initVelocity,'Children');
    dim = cellfun(@(x)str2double(x(end)),get(child,'Tag'));    
    set(child(dim<=userData.probDim),'Enable','on');
    set(child(dim>userData.probDim),'Enable','off');
else
    set(get(handles.uipanel_initVelocity,'Children'),'Enable','off');
end

if strcmpi(selectedButton,'radiobutton_convergePoint');
    child=get(handles.uipanel_convergePoint,'Children');
    dim = cellfun(@(x)str2double(x(end)),get(child,'Tag'));    
    set(child(dim<=userData.probDim),'Enable','on');
    set(child(dim>userData.probDim),'Enable','off');
else
    set(get(handles.uipanel_convergePoint,'Children'),'Enable','off');
end

%updates the handles structure
guidata(hObject, handles);
