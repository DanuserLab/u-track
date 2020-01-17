function varargout = trackingProcessGUI(varargin)
% TRACKINGPROCESSGUI M-file for trackingProcessGUI.fig
%      TRACKINGPROCESSGUI, by itself, creates a new TRACKINGPROCESSGUI or raises the existing
%      singleton*.
%
%      H = TRACKINGPROCESSGUI returns the handle to a new TRACKINGPROCESSGUI or the handle to
%      the existing singleton*.
%
%      TRACKINGPROCESSGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in TRACKINGPROCESSGUI.M with the given input arguments.
%
%      TRACKINGPROCESSGUI('Property','Value',...) creates a new TRACKINGPROCESSGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before trackingProcessGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to trackingProcessGUI_OpeningFcn via varargin.
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

% Edit the above text to modify the response to help trackingProcessGUI

% Last Modified by GUIDE v2.5 05-Mar-2012 18:28:45

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @trackingProcessGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @trackingProcessGUI_OutputFcn, ...
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


% --- Executes just before trackingProcessGUI is made visible.
function trackingProcessGUI_OpeningFcn(hObject, eventdata, handles, varargin)

processGUI_OpeningFcn(hObject, eventdata, handles, varargin{:},'initChannel',1);

% Parameter Setup
userData = get(handles.figure1, 'UserData');
if(isempty(userData)), userData = struct(); end;

funParams = userData.crtProc.funParams_;

set(handles.popupmenu_probDim,'String',{'2','3'},'UserData',[2 3],...
    'Value',find(funParams.probDim==[2 3]));
set(handles.checkbox_verbose, 'Value', funParams.verbose)

% gapCloseParam
set(handles.edit_maxgap, 'String', num2str(funParams.gapCloseParam.timeWindow - 1))
set(handles.edit_minlength, 'String', num2str(funParams.gapCloseParam.minTrackLen))
set(handles.checkbox_histogram, 'Value', funParams.gapCloseParam.diagnostics)

set(handles.checkbox_merging, 'Value',ismember(funParams.gapCloseParam.mergeSplit,[1 2]));
set(handles.checkbox_splitting, 'Value',ismember(funParams.gapCloseParam.mergeSplit,[1 3]));
    
% Set cost matrics
defaultLinkingCostMat = userData.crtProc.getDefaultLinkingCostMatrices(userData.MD,funParams.gapCloseParam.timeWindow);
defaultGapClosingCostMat = userData.crtProc.getDefaultGapClosingCostMatrices(userData.MD,funParams.gapCloseParam.timeWindow);
userData.cost_linking = {defaultLinkingCostMat.funcName};
userData.cost_gapclosing = {defaultGapClosingCostMat.funcName};
userData.fun_cost_linking = {defaultLinkingCostMat.GUI};
userData.fun_cost_gap = {defaultGapClosingCostMat.GUI};

% Retrieve index of default cost matrices
i1 = find(strcmp(funParams.costMatrices(1).funcName, userData.cost_linking));
i2 = find(strcmp(funParams.costMatrices(2).funcName, userData.cost_gapclosing));
assert(isscalar(i1) && isscalar(i2),'User-defined: the length of matching methods must be 1.')
nLinking=numel(defaultLinkingCostMat);
nGapClosing=numel(defaultGapClosingCostMat);
u1 = cell(1,nLinking);
u2 = cell(1,nGapClosing);
u1{i1} = funParams.costMatrices(1).parameters;
u2{i2} = funParams.costMatrices(2).parameters;
for i=setdiff(1:nLinking,i1), u1{i}=defaultLinkingCostMat(i).parameters; end
for i=setdiff(1:nGapClosing,i1), u2{i}=defaultGapClosingCostMat(i).parameters; end

set(handles.popupmenu_linking, 'Value', i1, 'UserData', u1,...
    'String',{defaultLinkingCostMat.name})
set(handles.popupmenu_gapclosing, 'Value', i2, 'UserData', u2,...
    'String',{defaultGapClosingCostMat.name})

% Kalman functions
userData.kalmanFunctions = TrackingProcess.getKalmanFunctions;
nKalmanFunctions = numel(userData.kalmanFunctions);
kalmanFields = {'reserveMem','initialize','calcGain','timeReverse'};

index=true(1,nKalmanFunctions);
for i=1:numel(kalmanFields)
    index=index & strcmp(funParams.kalmanFunctions.(kalmanFields{i}),...
        {userData.kalmanFunctions.(kalmanFields{i})});
end

assert(sum(index)==1, 'Did not find a unique Kalman set.');

u2 = cell(1,nKalmanFunctions);
u2{index} = funParams.costMatrices(1).parameters.kalmanInitParam;

set(handles.popupmenu_kalmanFunctions,'UserData',u2,...
    'String', {userData.kalmanFunctions.name}, 'Value', find(index))

set(handles.checkbox_export, 'Value', funParams.saveResults.export)

% Initialize children figure handles
userData.linkingFig=-1;
userData.gapclosingFig=-1;
userData.kalmanFig=-1;

% Choose default command line output for trackingProcessGUI
handles.output = hObject;

% Update user data and GUI data
set(hObject, 'UserData', userData);

uicontrol(handles.pushbutton_done);
guidata(hObject, handles);


% --- Outputs from this function are returned to the command line.
function varargout = trackingProcessGUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton_done.
function pushbutton_done_Callback(hObject, eventdata, handles)


% to make sure parameters are properly updated first if text edit box has been updated.
uicontrol(handles.checkbox_all) % force cursor change to trigger any callbacks -- hack
pause(.1);
qobj = findobj('Tag', 'Question'); % look for question dialoge from callbacks.
if ~isempty(qobj)
    drawnow
    uistack(qobj)
    return % require user to verify, then click apply once again to complete settings.
end
 
userData = get(handles.figure1, 'UserData');
if(isempty(userData)), userData = struct(); end;


% Check User Input
if isempty(get(handles.listbox_selectedChannels, 'String'))
    errordlg('Please select at least one input channel from ''Available Channels''.','Setting Error','modal')
    return;
end

props = get(handles.popupmenu_probDim, {'UserData','Value'});
probDim=props{1}(props{2});

timeWindow = str2double(get(handles.edit_maxgap, 'String'))+1;
if isnan(timeWindow) || timeWindow < 0 || floor(timeWindow) ~= ceil(timeWindow)
    errordlg('Please provide a valid value to parameter "Maximum Gap to Close".','Error','modal')
    return
end

minTrackLen = str2double(get(handles.edit_minlength, 'String'));
if isnan(minTrackLen) || minTrackLen < 0 || floor(minTrackLen) ~= ceil(minTrackLen)
    errordlg('Please provide a valid value to parameter "Minimum Length of Track Segment from First Step to use in Second Step".','Error','modal')
    return
end

% -------- Set parameter --------
channelIndex = get (handles.listbox_selectedChannels, 'Userdata');
funParams.ChannelIndex = channelIndex;

funParams.probDim = probDim;
funParams.verbose = get(handles.checkbox_verbose, 'Value');
funParams.gapCloseParam.timeWindow = timeWindow;
funParams.gapCloseParam.minTrackLen = minTrackLen;
funParams.gapCloseParam.diagnostics = get(handles.checkbox_histogram, 'Value');

if get(handles.checkbox_merging, 'Value') && get(handles.checkbox_splitting, 'Value')
    funParams.gapCloseParam.mergeSplit = 1;
elseif get(handles.checkbox_merging, 'Value') && ~get(handles.checkbox_splitting, 'Value')
    funParams.gapCloseParam.mergeSplit = 2;
elseif ~get(handles.checkbox_merging, 'Value') && get(handles.checkbox_splitting, 'Value')
    funParams.gapCloseParam.mergeSplit = 3;
elseif ~get(handles.checkbox_merging, 'Value') && ~get(handles.checkbox_splitting, 'Value')
    funParams.gapCloseParam.mergeSplit = 0;
end

funParams.saveResults.export = get(handles.checkbox_export, 'Value');

% Cost matrices
i_linking = get(handles.popupmenu_linking, 'Value');
i_gapclosing = get(handles.popupmenu_gapclosing, 'Value');

u_linking = get(handles.popupmenu_linking, 'UserData');
u_gapclosing = get(handles.popupmenu_gapclosing, 'UserData');

if isempty( u_linking{i_linking} )
    errordlg('Plese set up the selected cost function for "Step 1: frame-to-frame linking".','Error','modal')
end

if isempty( u_gapclosing{i_gapclosing} )
    errordlg('Plese set up the selected cost function for "Step 2: gap closing, mergin and splitting".','Error','modal')
end

funParams.costMatrices(1).funcName = userData.cost_linking{i_linking};
funParams.costMatrices(1).parameters = u_linking{i_linking};
funParams.costMatrices(2).funcName = userData.cost_gapclosing{i_gapclosing};
funParams.costMatrices(2).parameters = u_gapclosing{i_gapclosing};

% Get Kalman values
iKalman = get(handles.popupmenu_kalmanFunctions, 'Value');
kalmanFields = {'reserveMem','initialize','calcGain','timeReverse'};
for i=1:numel(kalmanFields)
    funParams.kalmanFunctions.(kalmanFields{i})=userData.kalmanFunctions(iKalman).(kalmanFields{i});
end
kalmanData = get(handles.popupmenu_kalmanFunctions, 'UserData');
funParams.costMatrices(1).parameters.kalmanInitParam = kalmanData{iKalman};

% Set up parameters effected by funParams.gapCloseParam.timeWindow
if isfield(funParams.costMatrices(2).parameters,'brownStdMult'),
    funParams.costMatrices(2).parameters.brownStdMult = funParams.costMatrices(2).parameters.brownStdMult(1) * ones(funParams.gapCloseParam.timeWindow,1);
end

if isfield(funParams.costMatrices(2).parameters,'linStdMult'),
    funParams.costMatrices(2).parameters.linStdMult = funParams.costMatrices(2).parameters.linStdMult(1) * ones(funParams.gapCloseParam.timeWindow,1);
end

processGUI_ApplyFcn(hObject,eventdata,handles,funParams)


% --- Executes on button press in pushbutton_cancel.
function pushbutton_cancel_Callback(hObject, eventdata, handles)
delete(handles.figure1);

% --- Executes on button press in pushbutton_set_linking.
function pushbutton_set_linking_Callback(hObject, eventdata, handles)
%       userData.linkingFig - the handle of setting panel for linking set-up
%       userData.gapclosingFig - the handle of setting panel for gap closing set-up
%       userData.kalmanFig

userData = get(handles.figure1, 'UserData');
if(isempty(userData)), userData = struct(); end;

parent = handles.popupmenu_linking;
procID = get(parent, 'Value');
if procID > length(userData.fun_cost_linking)
    warndlg('Please select a cost function for linking step.','Error','modal')
    return
else
    settingGUI = userData.fun_cost_linking{procID};
    userData.linkingFig = settingGUI(parent, procID);
end
set(handles.figure1, 'UserData', userData);

% --- Executes on button press in pushbutton_set_gapclosing.
function pushbutton_set_gapclosing_Callback(hObject, eventdata, handles)
%       userData.linkingFig - the handle of setting panel for linking set-up
%       userData.gapclosingFig - the handle of setting panel for gap closing set-up
%       userData.kalmanFig
userData = get(handles.figure1, 'UserData');
if(isempty(userData)), userData = struct(); end;

parent = handles.popupmenu_gapclosing;
procID = get(parent, 'Value');
if procID > length(userData.fun_cost_gap)
    warndlg('Please select a cost function for gap closing step.','Error','modal')
    return
else
    settingGUI = userData.fun_cost_gap{procID};
    userData.gapclosingFig = settingGUI(parent, procID);
end
set(handles.figure1, 'UserData', userData);



function edit_maxgap_Callback(hObject, eventdata, handles)

maxgap = str2double(get(handles.edit_maxgap, 'String'));
if isnan(maxgap) || maxgap < 0 || floor(maxgap) ~= ceil(maxgap)
    errordlg('Please provide a valid value to parameter "Maximum Gap to Close".','Warning','modal')
    return;
end

timeWindow = maxgap + 1; % Retrieve the new value for the time window

% Retrieve the parameters of the linking and gap closing matrices
u_linking = get(handles.popupmenu_linking, 'UserData');
linkingID = get(handles.popupmenu_linking, 'Value');
linkingParameters = u_linking{linkingID};
u_gapclosing = get(handles.popupmenu_gapclosing, 'UserData');
gapclosingID = get(handles.popupmenu_gapclosing, 'Value');
gapclosingParameters = u_gapclosing{gapclosingID};

% Check for changes
linkingnnWindowChange=(linkingParameters.nnWindow~=timeWindow);
gapclosingnnWindowChange=isfield(gapclosingParameters,'nnWindow') && (gapclosingParameters.nnWindow~=timeWindow);
gapclosingtimeReachConfBChange=isfield(gapclosingParameters,'timeReachConfB') && (gapclosingParameters.timeReachConfB~=timeWindow);
gapclosingtimeReachConfLChange=isfield(gapclosingParameters,'timeReachConfL') &&(gapclosingParameters.timeReachConfL~=timeWindow);

if ~linkingnnWindowChange && ~gapclosingnnWindowChange && ...
        ~gapclosingtimeReachConfBChange && ~gapclosingtimeReachConfLChange,
    return;
end
% Optional: asks the user if the time window value should be propagated
% to the linking and gap closing matrics
modifyParameters=questdlg('Do you want to propagate the changes in the maximum number of gaps to close?',...
    'Parameters update','Yes','No','Yes');
if ~strcmp(modifyParameters,'Yes'), return; end


% Set new linking time window
if linkingnnWindowChange, linkingParameters.nnWindow=timeWindow; end
u_linking{linkingID} = linkingParameters;
set(handles.popupmenu_linking, 'UserData', u_linking)

% Set new gap closing time window
if gapclosingnnWindowChange, gapclosingParameters.nnWindow=timeWindow; end
if gapclosingtimeReachConfBChange, gapclosingParameters.timeReachConfB=timeWindow; end
if gapclosingtimeReachConfLChange, gapclosingParameters.timeReachConfL=timeWindow; end
u_gapclosing{gapclosingID} = gapclosingParameters;
set(handles.popupmenu_gapclosing, 'UserData', u_gapclosing)
guidata(hObject,handles);


% --- Executes on button press in pushbutton_kalman_initialize.
function pushbutton_kalman_initialize_Callback(hObject, eventdata, handles)

userData = get(handles.figure1, 'UserData');
if(isempty(userData)), userData = struct(); end;

parent = handles.popupmenu_kalmanFunctions;
id = get(parent, 'Value');

settingGUI = userData.kalmanFunctions(id).initializeGUI;
userData.kalmanFig = settingGUI(parent, id);

set(handles.figure1, 'UserData', userData);

% --- Executes on key press with focus on figure1 and none of its controls.
function figure1_KeyPressFcn(hObject, eventdata, handles)

if strcmp(eventdata.Key, 'return')
    pushbutton_done_Callback(handles.pushbutton_done, [], handles);
end


% --- Executes during object deletion, before destroying properties.
function figure1_DeleteFcn(hObject, eventdata, handles)

userData = get(handles.figure1, 'UserData');
if(isempty(userData)), userData = struct(); end;

% Delete setting panels
if ishandle(userData.linkingFig), delete(userData.linkingFig);end
if ishandle(userData.gapclosingFig), delete(userData.gapclosingFig); end
if ishandle(userData.kalmanFig), delete(userData.kalmanFig); end


% --- Executes on button press in checkbox_export.
function checkbox_export_Callback(hObject, eventdata, handles)

if get(hObject,'Value')
    exportMsg=sprintf('The output matrices resulting from this process might be very large. Be cautious if you have large movies');
    if any([get(handles.checkbox_merging, 'Value') get(handles.checkbox_splitting, 'Value')])
        exportMsg =[exportMsg sprintf('\n \nAny merging and splitting information will be lost in the exported format.')];
    end
    warndlg(exportMsg,'Warning','modal')
end


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)

delete(hObject);
