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
% Copyright (C) 2025, Danuser Lab - UTSouthwestern 
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

% Last Modified by GUIDE v2.5 25-Oct-2019 16:13:24

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

% add new param 2019-08
if userData.MD.is3D
    set(handles.checkbox_exportTrackabilityData, 'Value', funParams.saveResults.exportTrackabilityData)
else
    delete(handles.checkbox_exportTrackabilityData);
end
    
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


if isequal(userData.procConstr, @TrackingDynROIProcess)
  % Set up available Build Dyn ROI channels
  set(handles.listbox_availableDynROIChannels,'String',userData.MD.getChannelPaths(), ...
      'UserData',1:numel(userData.MD.channels_));

  DynROIChannelIndex = funParams.buildDynROIProcessChannel;

  if ~isempty(DynROIChannelIndex)
      DynROIChannelString = userData.MD.getChannelPaths(DynROIChannelIndex);
  else
      DynROIChannelString = {};
  end
  set(handles.listbox_selectedDynROIChannels,'String',DynROIChannelString,...
      'UserData',DynROIChannelIndex);


  %Setup Build Dyn ROI process list box
  DynROIProc =  cellfun(@(x) isa(x,'BuildDynROIProcess'),userData.MD.processes_);
  DynROIProcID=find(DynROIProc);
  DynROIProcNames = cellfun(@(x) x.getName(),userData.MD.processes_(DynROIProc),'Unif',false);
  DynROIProcString = vertcat('Choose later',DynROIProcNames(:));
  DynROIProcData=horzcat({[]},num2cell(DynROIProcID));
  DynROIProcValue = find(cellfun(@(x) isequal(x,funParams.processBuildDynROI),userData.MD.processes_(DynROIProc)));
  if isempty(DynROIProcValue) && isempty(DynROIProcID)
      DynROIProcValue = 1; 
  elseif isempty(DynROIProcValue) && ~isempty(DynROIProcID) % make first available DynROIProc selected&set on the GUI, even funParams.processBuildDynROI = [].
      DynROIProcValue = 2;
  else
      DynROIProcValue = DynROIProcValue+1; 
  end
  set(handles.popupmenu_BuildDynROIProcessIndex,'String',DynROIProcString,...
      'UserData',DynROIProcData,'Value',DynROIProcValue);

  % Update channels listboxes depending on the selected process
  popupmenu_BuildDynROIProcessIndex_Callback(hObject, eventdata, handles)
else
  uipanel_DynROIProc_posi = get(handles.uipanel_DynROIProc, 'Position');
  hgtDiff = uipanel_DynROIProc_posi(4) + 7;
  delete(handles.uipanel_DynROIProc);
  set(handles.uipanel_kalman,'position', (get(handles.uipanel_kalman,'position') - [0 hgtDiff 0 0]));
  set(handles.uipanel_costfunction,'position', (get(handles.uipanel_costfunction,'position') - [0 hgtDiff 0 0]));
  set(handles.uipanel_parameters,'position', (get(handles.uipanel_parameters,'position') - [0 hgtDiff 0 0]));
  set(handles.uipanel_input,'position', (get(handles.uipanel_input,'position') - [0 hgtDiff 0 0]));
  set(handles.text_processName,'position', (get(handles.text_processName,'position') - [0 hgtDiff 0 0]));
  set(handles.axes_help,'position', (get(handles.axes_help,'position') - [0 hgtDiff 0 0]));
  set(handles.text_copyright,'position', (get(handles.text_copyright,'position') - [0 hgtDiff 0 0]));
  set(handles.figure1, 'Position', (get(handles.figure1,'position') - [0 0 0 hgtDiff]));  
end


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
if userData.MD.is3D
    funParams.saveResults.exportTrackabilityData = get(handles.checkbox_exportTrackabilityData, 'Value');
end

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

% Retrieve GUI-defined parameters for Build Dyn ROI process:
if isequal(userData.procConstr, @TrackingDynROIProcess)
  %Get selected DynROI channels
  DynROIChannelProps = get(handles.listbox_selectedDynROIChannels, {'Userdata','String'});
  funParams.buildDynROIProcessChannel = DynROIChannelProps{1};
  % Retrieve Build Dyn ROI process
  props=get(handles.popupmenu_BuildDynROIProcessIndex,{'UserData','Value'});
  DynROIProcessIndex = props{1}{props{2}};
  if ~isempty(DynROIProcessIndex)
    funParams.processBuildDynROI = userData.MD.processes_{DynROIProcessIndex};
  else 
    funParams.processBuildDynROI = [];
  end
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


% --- Executes on selection change in listbox_availableDynROIChannels.
function listbox_availableDynROIChannels_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_availableDynROIChannels (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox_availableDynROIChannels contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_availableDynROIChannels


% --- Executes during object creation, after setting all properties.
function listbox_availableDynROIChannels_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox_availableDynROIChannels (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox_DynROI_all.
function checkbox_DynROI_all_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_DynROI_all (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_DynROI_all
contents1 = get(handles.listbox_availableDynROIChannels, 'String');

chanIndex1 = get(handles.listbox_availableDynROIChannels, 'Userdata');
chanIndex2 = get(handles.listbox_selectedDynROIChannels, 'Userdata');

% Return if listbox1 is empty
if isempty(contents1)
    return;
end

switch get(hObject,'Value')
    case 1
        set(handles.listbox_selectedDynROIChannels, 'String', contents1);
        chanIndex2 = chanIndex1;
    case 0
        set(handles.listbox_selectedDynROIChannels, 'String', {}, 'Value',1);
        chanIndex2 = [ ];
end
set(handles.listbox_selectedDynROIChannels, 'UserData', chanIndex2);

% --- Executes on button press in pushbutton_DynROI_select.
function pushbutton_DynROI_select_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_DynROI_select (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
contents1 = get(handles.listbox_availableDynROIChannels, 'String');
contents2 = get(handles.listbox_selectedDynROIChannels, 'String');
id = get(handles.listbox_availableDynROIChannels, 'Value');

% If channel has already been added, return;
chanIndex1 = get(handles.listbox_availableDynROIChannels, 'Userdata');
chanIndex2 = get(handles.listbox_selectedDynROIChannels, 'Userdata');

for i = id

        contents2{end+1} = contents1{i};
        
        chanIndex2 = cat(2, chanIndex2, chanIndex1(i));

end

set(handles.listbox_selectedDynROIChannels, 'String', contents2, 'Userdata', chanIndex2);


% --- Executes on button press in pushbutton_DynROI_delete.
function pushbutton_DynROI_delete_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_DynROI_delete (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Call back function of 'delete' button
contents = get(handles.listbox_selectedDynROIChannels,'String');
id = get(handles.listbox_selectedDynROIChannels,'Value');

% Return if list is empty
if isempty(contents) || isempty(id)
    return;
end

% Delete selected item
contents(id) = [ ];

% Delete userdata
chanIndex2 = get(handles.listbox_selectedDynROIChannels, 'Userdata');
chanIndex2(id) = [ ];
set(handles.listbox_selectedDynROIChannels, 'Userdata', chanIndex2);

% Point 'Value' to the second last item in the list once the 
% last item has been deleted
if (id >length(contents) && id>1)
    set(handles.listbox_selectedDynROIChannels,'Value',length(contents));
end
% Refresh listbox
set(handles.listbox_selectedDynROIChannels,'String',contents);


% --- Executes on selection change in listbox_selectedDynROIChannels.
function listbox_selectedDynROIChannels_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_selectedDynROIChannels (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox_selectedDynROIChannels contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_selectedDynROIChannels


% --- Executes during object creation, after setting all properties.
function listbox_selectedDynROIChannels_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox_selectedDynROIChannels (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_DynROI_up.
function pushbutton_DynROI_up_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_DynROI_up (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% call back of 'Up' button

id = get(handles.listbox_selectedDynROIChannels,'Value');
contents = get(handles.listbox_selectedDynROIChannels,'String');


% Return if list is empty
if isempty(contents) || isempty(id) || id == 1
    return;
end

temp = contents{id};
contents{id} = contents{id-1};
contents{id-1} = temp;

chanIndex = get(handles.listbox_selectedDynROIChannels, 'Userdata');
temp = chanIndex(id);
chanIndex(id) = chanIndex(id-1);
chanIndex(id-1) = temp;

set(handles.listbox_selectedDynROIChannels, 'String', contents, 'Userdata', chanIndex);
set(handles.listbox_selectedDynROIChannels, 'value', id-1);


% --- Executes on button press in pushbutton_DynROI_down.
function pushbutton_DynROI_down_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_DynROI_down (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

id = get(handles.listbox_selectedDynROIChannels,'Value');
contents = get(handles.listbox_selectedDynROIChannels,'String');

% Return if list is empty
if isempty(contents) || isempty(id) || id == length(contents)
    return;
end

temp = contents{id};
contents{id} = contents{id+1};
contents{id+1} = temp;

chanIndex = get(handles.listbox_selectedDynROIChannels, 'Userdata');
temp = chanIndex(id);
chanIndex(id) = chanIndex(id+1);
chanIndex(id+1) = temp;

set(handles.listbox_selectedDynROIChannels, 'string', contents, 'Userdata',chanIndex);
set(handles.listbox_selectedDynROIChannels, 'value', id+1);


% --- Executes on selection change in popupmenu_BuildDynROIProcessIndex.
function popupmenu_BuildDynROIProcessIndex_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_BuildDynROIProcessIndex (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_BuildDynROIProcessIndex contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_BuildDynROIProcessIndex

% Retrieve selected process ID
props= get(handles.popupmenu_BuildDynROIProcessIndex,{'UserData','Value'});
procID = props{1}{props{2}};

% Read process and check available channels
userData = get(handles.figure1, 'UserData');
if(isempty(userData)), userData = struct(); end;

if isempty(procID)
    allChannelIndex=1:numel(userData.MD.channels_);
else
    allChannelIndex = find(userData.MD.processes_{procID}.checkChannelOutput);
end

% Set up available channels listbox
if ~isempty(allChannelIndex)
    if isempty(procID)
        channelString = userData.MD.getChannelPaths(allChannelIndex);
    else
        channelString = userData.MD.processes_{procID}.outFilePaths_(1,allChannelIndex);
    end
else
    channelString = {};
end
set(handles.listbox_availableDynROIChannels,'String',channelString,'UserData',allChannelIndex);

% Set up selected channels listbox
channelIndex = get(handles.listbox_selectedDynROIChannels, 'UserData');
channelIndex(~ismember(channelIndex,allChannelIndex)) = [];%So that indices may repeat, and handles empty better than intersect
if ~isempty(channelIndex)
    if isempty(procID)
        channelString = userData.MD.getChannelPaths(channelIndex);
    else
        channelString = userData.MD.processes_{procID}.outFilePaths_(1,channelIndex);
    end
else
    channelString = {};
    channelIndex = [];%Because the intersect command returns a 0x1 instead of 0x0 which causes concatenation errors
end
set(handles.listbox_selectedDynROIChannels,'String',channelString,'UserData',channelIndex);


% --- Executes during object creation, after setting all properties.
function popupmenu_BuildDynROIProcessIndex_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_BuildDynROIProcessIndex (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%% other functions added automatically by GUIDE on 10/25/2019:

% --- Executes on button press in checkbox_applytoall.
function checkbox_applytoall_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_applytoall (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_applytoall


% --- Executes on selection change in popupmenu_kalmanFunctions.
function popupmenu_kalmanFunctions_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_kalmanFunctions (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_kalmanFunctions contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_kalmanFunctions


% --- Executes during object creation, after setting all properties.
function popupmenu_kalmanFunctions_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_kalmanFunctions (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu_linking.
function popupmenu_linking_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_linking (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_linking contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_linking


% --- Executes during object creation, after setting all properties.
function popupmenu_linking_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_linking (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu_gapclosing.
function popupmenu_gapclosing_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_gapclosing (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_gapclosing contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_gapclosing


% --- Executes during object creation, after setting all properties.
function popupmenu_gapclosing_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_gapclosing (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu_probDim.
function popupmenu_probDim_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_probDim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_probDim contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_probDim


% --- Executes during object creation, after setting all properties.
function popupmenu_probDim_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_probDim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function edit_maxgap_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_maxgap (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_minlength_Callback(hObject, eventdata, handles)
% hObject    handle to edit_minlength (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_minlength as text
%        str2double(get(hObject,'String')) returns contents of edit_minlength as a double


% --- Executes during object creation, after setting all properties.
function edit_minlength_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_minlength (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox_merging.
function checkbox_merging_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_merging (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_merging


% --- Executes on button press in checkbox_splitting.
function checkbox_splitting_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_splitting (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_splitting


% --- Executes on button press in checkbox_histogram.
function checkbox_histogram_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_histogram (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_histogram


% --- Executes on button press in checkbox_verbose.
function checkbox_verbose_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_verbose (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_verbose


% --- Executes on button press in checkbox_exportTrackabilityData.
function checkbox_exportTrackabilityData_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_exportTrackabilityData (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_exportTrackabilityData


% --- Executes on selection change in listbox_availableChannels.
function listbox_availableChannels_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_availableChannels (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox_availableChannels contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_availableChannels


% --- Executes during object creation, after setting all properties.
function listbox_availableChannels_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox_availableChannels (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in listbox_selectedChannels.
function listbox_selectedChannels_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_selectedChannels (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox_selectedChannels contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_selectedChannels


% --- Executes during object creation, after setting all properties.
function listbox_selectedChannels_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox_selectedChannels (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox_all.
function checkbox_all_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_all (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_all


% --- Executes on button press in pushbutton_select.
function pushbutton_select_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_select (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton_delete.
function pushbutton_delete_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_delete (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
