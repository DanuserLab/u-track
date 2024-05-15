function varargout = pointSourceDetectionProcessGUI3D(varargin)
% anisoGaussianDetectionProcessGUI M-file for anisoGaussianDetectionProcessGUI.fig
%      anisoGaussianDetectionProcessGUI, by itself, creates a new pointSourceDetectionProcessGUI3D or raises the existing
%      singleton*.
%
%      H = anisoGaussianDetectionProcessGUI returns the handle to a new pointSourceDetectionProcessGUI3D or the handle to
%      the existing singleton*.
%
%      anisoGaussianDetectionProcessGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in pointSourceDetectionProcessGUI3D.M with the given input arguments.
%
%      anisoGaussianDetectionProcessGUI('Property','Value',...) creates a new pointSourceDetectionProcessGUI3D or raises the
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
% Copyright (C) 2024, Danuser Lab - UTSouthwestern 
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

% Last Modified by GUIDE v2.5 09-Jan-2020 14:00:22

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @pointSourceDetectionProcessGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @pointSourceDetectionProcessGUI_OutputFcn, ...
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
function pointSourceDetectionProcessGUI_OpeningFcn(hObject, eventdata, handles, varargin)

processGUI_OpeningFcn(hObject, eventdata, handles, varargin{:},'initChannel',1);

% Set-up parameters
userData = get(handles.figure1,'UserData');
funParams = userData.crtProc.funParams_;
%Remove the output directory as we don't want to replicate it to other
%movies if the "apply to all movies" box is checked. Ideally we would
%explicitly only replicate the parameters we set in this GUI but this is a
%quick fix. - HLE
if isfield(funParams,'OutputDirectory')
    funParams = rmfield(funParams,'OutputDirectory');
end

set(handles.popupmenu_CurrentChannel,'UserData',funParams);

% Set Channel options
selChan = handles.listbox_selectedChannels.String;
availChan = handles.listbox_availableChannels.String;
[C ia is] = intersect(availChan, selChan);
% 
chanStr = arrayfun(@(x)(['Channel ' num2str(x)]), 1:numel(availChan),'Unif',0);
set(handles.popupmenu_CurrentChannel, 'String', chanStr);

iChan = get(handles.popupmenu_CurrentChannel,'Value');

if isempty(iChan) || iChan ~= ia(1)
    iChan = ia(1);
    set(handles.popupmenu_CurrentChannel,'Value',iChan);
end

handles.text_warningChan.Visible = 'off';


%%%%%%%%

ImProcP =  cellfun(@(x) isa(x,'ImageProcessingProcess'), userData.MD.processes_);
ImProcPID = find(ImProcP);
ImProcPNames = cellfun(@(x) x.getName(),userData.MD.processes_(ImProcP),'Unif',false);
ImProcPString = vertcat('Raw Channel Input', ImProcPNames(:));

imProcData = horzcat({[]},num2cell(ImProcPID));
imProcValue = find(cellfun(@(x) isequal(x,funParams.ChannelIndex),imProcData));
if isempty(imProcValue), imProcValue = 1; end

set(handles.edit_InputImageProcessIndex,'String', ImProcPString,...
    'UserData',imProcData,'Value',imProcValue);
%%%%%%%%%%%%%%%


% Set default channels callback function
set(handles.checkbox_all,'Callback',@(hObject,eventdata)...
    checkallChannels_Callback(hObject,eventdata,guidata(hObject)));
set(handles.pushbutton_select,'Callback',@(hObject,eventdata)...
    selectChannel_Callback(hObject,eventdata,guidata(hObject)));
set(handles.pushbutton_delete,'Callback',@(hObject,eventdata)...
    deleteChannel_Callback(hObject,eventdata,guidata(hObject)));


% Set up available mask channels
set(handles.listbox_availableMaskChannels,'String',userData.MD.getChannelPaths(), ...
    'UserData',1:numel(userData.MD.channels_));

maskChannelIndex = funParams.MaskChannelIndex;

if ~isempty(maskChannelIndex)
    maskChannelString = userData.MD.getChannelPaths(maskChannelIndex);
else
    maskChannelString = {};
end
set(handles.listbox_selectedMaskChannels,'String',maskChannelString,...
    'UserData',maskChannelIndex);


%Setup mask process list box
segProc =  cellfun(@(x) isa(x,'MaskProcess'),userData.MD.processes_);
segProcID=find(segProc);
segProcNames = cellfun(@(x) x.getName(),userData.MD.processes_(segProc),'Unif',false);
segProcString = vertcat('Choose later',segProcNames(:));
segProcData=horzcat({[]},num2cell(segProcID));
segProcValue = find(cellfun(@(x) isequal(x,funParams.MaskProcessIndex),segProcData));
if isempty(segProcValue), segProcValue = 1; end
set(handles.popupmenu_SegProcessIndex,'String',segProcString,...
    'UserData',segProcData,'Value',segProcValue);



if ~ismember(handles.listbox_availableChannels.String{iChan}, handles.listbox_selectedChannels.String)
    handles.uipanel_1.BackgroundColor = [1 .08 .18];
    handles.text_warningChan.Visible = 'off';
else
    handles.uipanel_1.BackgroundColor = [.94 .94 .94];
    handles.text_warningChan.Visible = 'on';
end

% Initialize non-channel specific parameters
% IsoCoord
set(handles.(['edit_isoCoord']), 'Value', funParams.isoCoord);

% Update channels listboxes depending on the selected process
popupmenu_SegProcessIndex_Callback(hObject, eventdata, handles)

%Update channel parameter selection dropdown
popupmenu_CurrentChannel_Callback(hObject, eventdata, handles)

% Add Scales parameter to the GUI
set(handles.edit_scales, 'String',num2str(funParams.scales))


if isequal(userData.procConstr, @PointSourceDetectionProcess3DDynROI)
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
  widthDiff = uipanel_DynROIProc_posi(3);
  delete(handles.uipanel_DynROIProc);
  set(handles.figure1, 'Position', (get(handles.figure1,'position') - [0 0 widthDiff 0]));
end


% Update GUI user data
set(handles.figure1, 'UserData', userData);
handles.output = hObject;
guidata(hObject, handles);

% --- Outputs from this function are returned to the command line.
function varargout = pointSourceDetectionProcessGUI_OutputFcn(~, ~, handles) 
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

if any(isnan(str2num(get(handles.edit_scales, 'String')))) ...
    || any(str2num(get(handles.edit_scales, 'String')) < 0) ...
    || isempty(get(handles.edit_scales, 'String'))
  errordlg('Please provide a valid input for ''Scales''.','Setting Error','modal');
  return;
end

%Save the currently set per-channel parameters
pushbutton_saveChannelParams_Callback(hObject, eventdata, handles)


% Retrieve detection parameters
funParams = get(handles.popupmenu_CurrentChannel,'UserData');

% Retrieve GUI-defined non-channel specific parameters

%Get selected image channels
channelIndex = get(handles.listbox_selectedChannels, 'Userdata');
if isempty(channelIndex)
    errordlg('Please select at least one input channel from ''Available Channels''.','Setting Error','modal')
    return;
end
funParams.ChannelIndex = channelIndex;

%Get selected mask channels
maskChannelProps = get(handles.listbox_selectedMaskChannels, {'Userdata','String'});

if ~isempty(maskChannelProps{1}) && ( numel(maskChannelProps{1}) ~= numel(channelIndex))
    errordlg('Please select either zero mask channels or the same number of mask channels as input channels.','Setting Error','modal')
    return;
end

funParams.MaskChannelIndex = maskChannelProps{1};

% Retrieve mask process index and class (for propagation)
props=get(handles.popupmenu_SegProcessIndex,{'UserData','Value'});
funParams.MaskProcessIndex = props{1}{props{2}};

funParams.UseIntersection = get(handles.edit_UseIntersection,'Value') > 0;

% Retrieve GUI-defined parameters for Build Dyn ROI process:
userData = get(handles.figure1,'UserData');

if isequal(userData.procConstr, @PointSourceDetectionProcess3DDynROI)
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

% Add Scales parameter to the GUI
funParams.scales = str2num(get(handles.edit_scales, 'String'));

% Add 64-bit warning
is64bit = ~isempty(regexp(computer ,'64$', 'once'));
if ~is64bit
    warndlg(['Your Matlab version is not detected as 64-bit. Please note '....
        'the point source detection uses compiled MEX files which '...
        'are not provided for 32-bit.'],...
        'Setting Error','modal');
end

processGUI_ApplyFcn(hObject, eventdata, handles,funParams);



function edit_Mode_Callback(hObject, eventdata, handles)
% hObject    handle to edit_Mode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_Mode as text
%        str2double(get(hObject,'String')) returns contents of edit_Mode as a double


% --- Executes during object creation, after setting all properties.
function edit_Mode_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_Mode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in edit_FitMixtures.
function edit_FitMixtures_Callback(hObject, eventdata, handles)
% hObject    handle to edit_FitMixtures (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of edit_FitMixtures



function edit_MaxMixtures_Callback(hObject, eventdata, handles)
% hObject    handle to edit_MaxMixtures (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_MaxMixtures as text
%        str2double(get(hObject,'String')) returns contents of edit_MaxMixtures as a double


% --- Executes during object creation, after setting all properties.
function edit_MaxMixtures_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_MaxMixtures (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_RedundancyRadius_Callback(hObject, eventdata, handles)
% hObject    handle to edit_RedundancyRadius (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_RedundancyRadius as text
%        str2double(get(hObject,'String')) returns contents of edit_RedundancyRadius as a double


% --- Executes during object creation, after setting all properties.
function edit_RedundancyRadius_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_RedundancyRadius (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in listbox_availableMaskChannels.
function listbox_availableMaskChannels_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_availableMaskChannels (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox_availableMaskChannels contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_availableMaskChannels


% --- Executes during object creation, after setting all properties.
function listbox_availableMaskChannels_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox_availableMaskChannels (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox_mask_all.
function checkbox_mask_all_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_mask_all (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_mask_all
contents1 = get(handles.listbox_availableMaskChannels, 'String');

chanIndex1 = get(handles.listbox_availableMaskChannels, 'Userdata');
chanIndex2 = get(handles.listbox_selectedMaskChannels, 'Userdata');

% Return if listbox1 is empty
if isempty(contents1)
    return;
end

switch get(hObject,'Value')
    case 1
        set(handles.listbox_selectedMaskChannels, 'String', contents1);
        chanIndex2 = chanIndex1;
    case 0
        set(handles.listbox_selectedMaskChannels, 'String', {}, 'Value',1);
        chanIndex2 = [ ];
end
set(handles.listbox_selectedMaskChannels, 'UserData', chanIndex2);


% --- Executes on button press in pushbutton_mask_select.
function pushbutton_mask_select_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_mask_select (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

contents1 = get(handles.listbox_availableMaskChannels, 'String');
contents2 = get(handles.listbox_selectedMaskChannels, 'String');
id = get(handles.listbox_availableMaskChannels, 'Value');

% If channel has already been added, return;
chanIndex1 = get(handles.listbox_availableMaskChannels, 'Userdata');
chanIndex2 = get(handles.listbox_selectedMaskChannels, 'Userdata');

for i = id

        contents2{end+1} = contents1{i};
        
        chanIndex2 = cat(2, chanIndex2, chanIndex1(i));

end

set(handles.listbox_selectedMaskChannels, 'String', contents2, 'Userdata', chanIndex2);

% --- Executes on button press in pushbutton_mask_delete.
function pushbutton_mask_delete_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_mask_delete (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Call back function of 'delete' button
contents = get(handles.listbox_selectedMaskChannels,'String');
id = get(handles.listbox_selectedMaskChannels,'Value');

% Return if list is empty
if isempty(contents) || isempty(id)
    return;
end

% Delete selected item
contents(id) = [ ];

% Delete userdata
chanIndex2 = get(handles.listbox_selectedMaskChannels, 'Userdata');
chanIndex2(id) = [ ];
set(handles.listbox_selectedMaskChannels, 'Userdata', chanIndex2);

% Point 'Value' to the second last item in the list once the 
% last item has been deleted
if (id >length(contents) && id>1)
    set(handles.listbox_selectedMaskChannels,'Value',length(contents));
end
% Refresh listbox
set(handles.listbox_selectedMaskChannels,'String',contents);


% --- Executes on selection change in listbox_selectedMaskChannels.
function listbox_selectedMaskChannels_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_selectedMaskChannels (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox_selectedMaskChannels contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_selectedMaskChannels


% --- Executes during object creation, after setting all properties.
function listbox_selectedMaskChannels_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox_selectedMaskChannels (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_up.
function pushbutton_up_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_up (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% call back of 'Up' button

id = get(handles.listbox_selectedMaskChannels,'Value');
contents = get(handles.listbox_selectedMaskChannels,'String');


% Return if list is empty
if isempty(contents) || isempty(id) || id == 1
    return;
end

temp = contents{id};
contents{id} = contents{id-1};
contents{id-1} = temp;

chanIndex = get(handles.listbox_selectedMaskChannels, 'Userdata');
temp = chanIndex(id);
chanIndex(id) = chanIndex(id-1);
chanIndex(id-1) = temp;

set(handles.listbox_selectedMaskChannels, 'String', contents, 'Userdata', chanIndex);
set(handles.listbox_selectedMaskChannels, 'value', id-1);


% --- Executes on button press in pushbutton_down.
function pushbutton_down_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_down (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


id = get(handles.listbox_selectedMaskChannels,'Value');
contents = get(handles.listbox_selectedMaskChannels,'String');

% Return if list is empty
if isempty(contents) || isempty(id) || id == length(contents)
    return;
end

temp = contents{id};
contents{id} = contents{id+1};
contents{id+1} = temp;

chanIndex = get(handles.listbox_selectedMaskChannels, 'Userdata');
temp = chanIndex(id);
chanIndex(id) = chanIndex(id+1);
chanIndex(id+1) = temp;

set(handles.listbox_selectedMaskChannels, 'string', contents, 'Userdata',chanIndex);
set(handles.listbox_selectedMaskChannels, 'value', id+1);


% --- Executes on selection change in popupmenu_SegProcessIndex.
function popupmenu_SegProcessIndex_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_SegProcessIndex (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_SegProcessIndex contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_SegProcessIndex
% Retrieve selected process ID
props= get(handles.popupmenu_SegProcessIndex,{'UserData','Value'});
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
set(handles.listbox_availableMaskChannels,'String',channelString,'UserData',allChannelIndex);

% Set up selected channels listbox
channelIndex = get(handles.listbox_selectedMaskChannels, 'UserData');
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
set(handles.listbox_selectedMaskChannels,'String',channelString,'UserData',channelIndex);


% --- Executes during object creation, after setting all properties.
function popupmenu_SegProcessIndex_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_SegProcessIndex (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in edit_UseIntersection.
function edit_UseIntersection_Callback(hObject, eventdata, handles)
% hObject    handle to edit_UseIntersection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of edit_UseIntersection


% --- Executes on selection change in popupmenu_CurrentChannel.
function popupmenu_CurrentChannel_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_CurrentChannel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_CurrentChannel contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_CurrentChannel
userData = get(handles.figure1,'UserData');
funParams = get(handles.popupmenu_CurrentChannel,'UserData');

iChan = get(handles.popupmenu_CurrentChannel,'Value');

if ~ismember(handles.listbox_availableChannels.String{iChan}, handles.listbox_selectedChannels.String)
    handles.uipanel_1.BackgroundColor = [.64 .08 .18];
    handles.text_warningChan.Visible = 'on';
else
    handles.uipanel_1.BackgroundColor = [.94 .94 .94];
    handles.text_warningChan.Visible = 'off';
end

% Set-up parameters
for i =1 : numel(funParams.PerChannelParams)
    paramName = funParams.PerChannelParams{i};
    if ~strcmp(paramName,'filterSigma')  && ~strcmp(paramName,'InputImageProcessIndex') && ~strcmp(paramName,'algorithmType')
        parVal = funParams.(paramName)(iChan);
        if ~islogical(funParams.(paramName))
            set(handles.(['edit_' paramName]), 'String', parVal);
        else
            set(handles.(['edit_' paramName]), 'Value', parVal);
        end
    elseif strcmp(paramName,'InputImageProcessIndex')

        selProcIndx = find(cellfun(@(x) ~(strcmp(class(x),'PointSourceDetectionProcess3D') || strcmp(class(x),'ComputeMIPProcess')), userData.MD.processes_));        
        if ~isempty(selProcIndx) && length(selProcIndx) >= 1
            ProcStr = [('Raw Channel Input')];
            ProcStr = [ProcStr arrayfun(@(x)([userData.MD.processes_{x}.name_]), selProcIndx, 'Unif',0)];
        else
            ProcStr = ['Raw Channel Input'];% arrayfun(@(x)(['NA ']),selProc,'Unif',0);
        end
        
        parVal = funParams.(paramName)(iChan);
        selProcIndx = [0 selProcIndx];
        iSelProcIn  = find(ismember(selProcIndx, parVal));
        if isempty(iSelProcIn)
            iSelProcIn = 1;
        end
        set(handles.(['edit_' paramName]), 'String', ProcStr);
        set(handles.(['edit_' paramName]), 'Value', iSelProcIn);
    elseif strcmp(paramName,'algorithmType') 
        set(handles.(['edit_' paramName]), 'String', PointSourceDetectionProcess3D.getValidAlgorithmTypes); 
        parVal = funParams.(paramName)(iChan);
        valSel  = find(ismember(PointSourceDetectionProcess3D.getValidAlgorithmTypes, parVal));
        set(handles.(['edit_' paramName]), 'Value', valSel); 
    else
        filterSigmaXY = funParams.filterSigma(1,iChan);
        filterSigmaZ = funParams.filterSigma(2,iChan);
        set(handles.('edit_filterSigmaZ'), 'String', filterSigmaZ);
        set(handles.('edit_filterSigmaXY'), 'String', filterSigmaXY);
    end
end

selType = get(handles.edit_algorithmType, 'Value'); 
algoType = PointSourceDetectionProcess3D.getValidAlgorithmTypes{selType};


set(handles.edit_isoCoord, 'enable','on');

if any(ismember(algoType,{'watershedApplegateAuto', ...
                      'watershedApplegate',...
                      'bandPassWatershed',...
                      'watershedMatlab',...
                      'markedWatershed'}))
    
    children = get(handles.uipanel_pointSource,'Children');
    set(children(strcmpi ( get (children,'Type'),'UIControl')),'enable','off')
    
    children = get(handles.uipanel_water,'Children');
    set(children(strcmpi ( get (children,'Type'),'UIControl')),'enable','on')
    
elseif any(ismember(algoType,{'pointSourceLM',...
                              'pointSource',...
                              'pointSourceAutoSigma',...
                              'pointSourceAutoSigmaFit',...
                              'pSAutoSigmaMarkedWatershed',...
                              'pointSourceAutoSigmaMixture',... 
                              'pointSourceAutoSigmaLM',...     
                              'pointSourceAutoSigmaFitSig',... 
                              'pSAutoSigmaWatershed'}))

    children = get(handles.uipanel_pointSource,'Children');
    set(children(strcmpi ( get (children,'Type'),'UIControl')),'enable','on')
    
    children = get(handles.uipanel_water,'Children');
    set(children(strcmpi ( get (children,'Type'),'UIControl')),'enable','off')

elseif any(ismember(algoType,{'multiscaleDetectionDebug'}))

    children = get(handles.uipanel_pointSource,'Children');
    set(children(strcmpi ( get (children,'Type'),'UIControl')),'enable','off');
    children = get(handles.uipanel_water,'Children');
    set(children(strcmpi ( get (children,'Type'),'UIControl')),'enable','off');

    set(handles.text_alpha, 'enable','on');
    set(handles.edit_alpha, 'enable','on');
    set(handles.text_scales, 'enable','on');
    set(handles.edit_scales, 'enable','on');
    set(handles.edit_isoCoord, 'enable','off');
                          
end





% --- Executes during object creation, after setting all properties.
function popupmenu_CurrentChannel_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_CurrentChannel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_saveChannelParams.
function pushbutton_saveChannelParams_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_saveChannelParams (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%Get settings for the current channel before switching to another
iChan = get(handles.popupmenu_CurrentChannel,'Value');

userData = get(handles.figure1,'UserData');
funParams = get(handles.popupmenu_CurrentChannel,'UserData');

for i =1 : numel(funParams.PerChannelParams)
    paramName = funParams.PerChannelParams{i};
    if ~strcmp(paramName,'filterSigma')  && ~strcmp(paramName,'InputImageProcessIndex') && ~strcmp(paramName,'algorithmType') 
        if islogical(funParams.(paramName))
            parVal = get(handles.(['edit_' paramName]), 'Value');
            funParams.(paramName)(iChan) = parVal;
        elseif iscell(funParams.(paramName))   
            parVal = get(handles.(['edit_' paramName]), 'String');
            funParams.(paramName)(iChan) = parVal;
        else
            parVal = get(handles.(['edit_' paramName]), 'String');
            funParams.(paramName)(iChan) = str2double(parVal);
        end
        
    elseif strcmp(paramName,'algorithmType')      
        selType = get(handles.(['edit_' paramName]), 'Value'); 
        parVal = PointSourceDetectionProcess3D.getValidAlgorithmTypes{selType};
        funParams.(paramName)(iChan) = {parVal};
        if isequal(parVal, 'multiscaleDetectionDebug')
            funParams.version = 'useMaxResponse';
        end
    
    elseif strcmp(paramName,'InputImageProcessIndex')
%         userData = get(handles.figure1,'UserData');
        MD = userData.MD;
        
%         selProc = 0:userData.procID-1;
selProcIndx = find(cellfun(@(x) ~(strcmp(class(x),'PointSourceDetectionProcess3D') || strcmp(class(x),'ComputeMIPProcess')), userData.MD.processes_));        
%         selProcIndx = find(cellfun(@(x) ~strcmp(class(x),'PointSourceDetectionProcess3D'), userData.MD.processes_));
        selProcIndx = [0 selProcIndx];
        iProcSel = get(handles.(['edit_' paramName]), 'Value'); 
        parVal = selProcIndx(iProcSel);
        
        if parVal ~= 0
            assert(MD.processes_{parVal}.checkChannelOutput(iChan),...
               ['No valid output for input process specified for channel ' num2str(iChan)]);
        end
        
        funParams.(paramName)(iChan) = parVal;
        
    else
        filterSigmaZ = get(handles.('edit_filterSigmaZ'), 'String');
        filterSigmaXY = get(handles.('edit_filterSigmaXY'), 'String');
        funParams.filterSigma(:,iChan) = [str2double(filterSigmaXY); str2double(filterSigmaZ)];
    end
end

% non-channel specific parameters
paramName = 'isoCoord';
parVal = get(handles.(['edit_' paramName]), 'Value');
funParams.(paramName) = parVal;

set(handles.popupmenu_CurrentChannel,'UserData',funParams);


% --- Executes on button press in edit_PreFilter.
function edit_PreFilter_Callback(hObject, eventdata, handles)
% hObject    handle to edit_PreFilter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of edit_PreFilter



function edit_ConfRadius_Callback(hObject, eventdata, handles)
% hObject    handle to edit_ConfRadius (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_ConfRadius as text
%        str2double(get(hObject,'String')) returns contents of edit_ConfRadius as a double


% --- Executes during object creation, after setting all properties.
function edit_ConfRadius_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_ConfRadius (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_WindowSize_Callback(hObject, eventdata, handles)
% hObject    handle to edit_WindowSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_WindowSize as text
%        str2double(get(hObject,'String')) returns contents of edit_WindowSize as a double


% --- Executes during object creation, after setting all properties.
function edit_WindowSize_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_WindowSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_filterSigmaXY_Callback(hObject, eventdata, handles)
% hObject    handle to edit_filterSigmaXY (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_filterSigmaXY as text
%        str2double(get(hObject,'String')) returns contents of edit_filterSigmaXY as a double


% --- Executes during object creation, after setting all properties.
function edit_filterSigmaXY_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_filterSigmaXY (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in edit_RefineMaskLoG.
function edit_RefineMaskLoG_Callback(hObject, eventdata, handles)
% hObject    handle to edit_RefineMaskLoG (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of edit_RefineMaskLoG


% --- Executes during object creation, after setting all properties.
function edit_RefineMaskLoG_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_RefineMaskLoG (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on button press in edit_RefineMaskValid.
function edit_RefineMaskValid_Callback(hObject, eventdata, handles)
% hObject    handle to edit_RefineMaskValid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of edit_RefineMaskValid



function edit_InputImageProcessIndex_Callback(hObject, eventdata, handles)
% hObject    handle to edit_InputImageProcessIndex (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_InputImageProcessIndex as text
%        str2double(get(hObject,'String')) returns contents of edit_InputImageProcessIndex as a double


% --- Executes during object creation, after setting all properties.
function edit_InputImageProcessIndex_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_InputImageProcessIndex (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in edit_algorithmType.
function edit_algorithmType_Callback(hObject, eventdata, handles)
% hObject    handle to edit_algorithmType (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns edit_algorithmType contents as cell array
%        contents{get(hObject,'Value')} returns selected item from edit_algorithmType
    selType = get(handles.edit_algorithmType, 'Value'); 
    algoType = PointSourceDetectionProcess3D.getValidAlgorithmTypes{selType};
 
    
set(handles.edit_isoCoord, 'enable','on');

if any(ismember(algoType,{'watershedApplegateAuto', ...
                      'watershedApplegate',...
                      'bandPassWatershed',...
                      'watershedMatlab',...
                      'markedWatershed'}))
    
    children = get(handles.uipanel_pointSource,'Children');
    set(children(strcmpi ( get (children,'Type'),'UIControl')),'enable','off');
    children = get(handles.uipanel_water,'Children');
    set(children(strcmpi ( get (children,'Type'),'UIControl')),'enable','on');
    
elseif any(ismember(algoType,{'pointSourceLM',...
                              'pointSource',...
                              'pointSourceAutoSigma',...
                              'pointSourceAutoSigmaFit',...
                              'pSAutoSigmaMarkedWatershed',...
                              'pointSourceAutoSigmaMixture',... 
                              'pointSourceAutoSigmaLM',...     
                              'pointSourceAutoSigmaFitSig',... 
                              'pSAutoSigmaWatershed'}))

    children = get(handles.uipanel_pointSource,'Children');
    set(children(strcmpi ( get (children,'Type'),'UIControl')),'enable','on');
    children = get(handles.uipanel_water,'Children');
    set(children(strcmpi ( get (children,'Type'),'UIControl')),'enable','off');
  
elseif any(ismember(algoType,{'multiscaleDetectionDebug'}))

    children = get(handles.uipanel_pointSource,'Children');
    set(children(strcmpi ( get (children,'Type'),'UIControl')),'enable','off');
    children = get(handles.uipanel_water,'Children');
    set(children(strcmpi ( get (children,'Type'),'UIControl')),'enable','off');

    set(handles.text_alpha, 'enable','on');
    set(handles.edit_alpha, 'enable','on');
    set(handles.text_scales, 'enable','on');
    set(handles.edit_scales, 'enable','on');
    set(handles.edit_isoCoord, 'enable','off');

end

        

    

% --- Executes during object creation, after setting all properties.
function edit_algorithmType_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_algorithmType (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_waterThresh_Callback(hObject, eventdata, handles)
% hObject    handle to edit_waterThresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_waterThresh as text
%        str2double(get(hObject,'String')) returns contents of edit_waterThresh as a double


% --- Executes during object creation, after setting all properties.
function edit_waterThresh_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_waterThresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_waterStep_Callback(hObject, eventdata, handles)
% hObject    handle to edit_waterStep (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_waterStep as text
%        str2double(get(hObject,'String')) returns contents of edit_waterStep as a double


% --- Executes during object creation, after setting all properties.
function edit_waterStep_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_waterStep (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_lowFreq_Callback(hObject, eventdata, handles)
% hObject    handle to edit_lowFreq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_lowFreq as text
%        str2double(get(hObject,'String')) returns contents of edit_lowFreq as a double


% --- Executes during object creation, after setting all properties.
function edit_lowFreq_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_lowFreq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_highFreq_Callback(hObject, eventdata, handles)
% hObject    handle to edit_highFreq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_highFreq as text
%        str2double(get(hObject,'String')) returns contents of edit_highFreq as a double


% --- Executes during object creation, after setting all properties.
function edit_highFreq_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_highFreq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function checkallChannels_Callback(hObject, eventdata, handles)

% Retrieve available channels properties
availableProps = get(handles.listbox_availableChannels, {'String','UserData'});
if isempty(availableProps{1}), return; end

% Update selected channels
if get(hObject,'Value')
    set(handles.listbox_selectedChannels, 'String', availableProps{1},...
        'UserData',availableProps{2});
else
    set(handles.listbox_selectedChannels, 'String', {}, 'UserData',[], 'Value',1);
end
% Set Channel options
% selChan = handles.listbox_selectedChannels.String;
% availChan = handles.listbox_availableChannels.String;
% [C ia is] = intersect(availChan, selChan);
% chanStr = arrayfun(@(x)(['Channel ' num2str(x)]), ia,'Unif',0);
% set(handles.popupmenu_CurrentChannel, 'String', chanStr);
popupmenu_CurrentChannel_Callback(hObject, eventdata, handles)


function deleteChannel_Callback(hObject, eventdata, handles)
% Generic callback to be exectuted when a selected channel is removed from
% the graphical settings interface

% Get selected properties and returin if empty
selectedProps = get(handles.listbox_selectedChannels, {'String','UserData','Value'});
if isempty(selectedProps{1}) || isempty(selectedProps{3}),return; end

% Delete selected item
selectedProps{1}(selectedProps{3}) = [ ];
selectedProps{2}(selectedProps{3}) = [ ];
set(handles.listbox_selectedChannels, 'String', selectedProps{1},'UserData',selectedProps{2},...
    'Value',max(1,min(selectedProps{3},numel(selectedProps{1}))));
% Set Channel options
% selChan = handles.listbox_selectedChannels.String;
% availChan = handles.listbox_availableChannels.String;
% [C ia is] = intersect(availChan, selChan);
% chanStr = arrayfun(@(x)(['Channel ' num2str(x)]), ia,'Unif',0);
% set(handles.popupmenu_CurrentChannel, 'String', chanStr);
popupmenu_CurrentChannel_Callback(hObject, eventdata, handles)

function selectChannel_Callback(hObject, eventdata, handles)

% Retrieve  channels properties
availableProps = get(handles.listbox_availableChannels, {'String','UserData','Value'});
selectedProps = get(handles.listbox_selectedChannels, {'String','UserData'});

% Find new elements and set them to the selected listbox
newID = availableProps{3}(~ismember(availableProps{1}(availableProps{3}),selectedProps{1}));
selectedChannels = horzcat(selectedProps{1}',availableProps{1}(newID)');
selectedData = horzcat(selectedProps{2}, availableProps{2}(newID));
set(handles.listbox_selectedChannels, 'String', selectedChannels, 'UserData', selectedData);

popupmenu_CurrentChannel_Callback(hObject, eventdata, handles)
% % Set Channel options
% selChan = handles.listbox_selectedChannels.String;
% availChan = handles.listbox_availableChannels.String;
% [C ia is] = intersect(availChan, selChan);
% chanStr = arrayfun(@(x)(['Channel ' num2str(x)]), ia,'Unif',0);
% set(handles.popupmenu_CurrentChannel, 'String', chanStr);


% --- Executes on button press in edit_isoCoord.
function edit_isoCoord_Callback(hObject, eventdata, handles)
% hObject    handle to edit_isoCoord (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of edit_isoCoord


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



function edit_alpha_Callback(hObject, eventdata, handles)
% hObject    handle to edit_alpha (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_alpha as text
%        str2double(get(hObject,'String')) returns contents of edit_alpha as a double


% --- Executes during object creation, after setting all properties.
function edit_alpha_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_alpha (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_filterSigmaZ_Callback(hObject, eventdata, handles)
% hObject    handle to edit_filterSigmaZ (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_filterSigmaZ as text
%        str2double(get(hObject,'String')) returns contents of edit_filterSigmaZ as a double


% --- Executes during object creation, after setting all properties.
function edit_filterSigmaZ_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_filterSigmaZ (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_scales_Callback(hObject, eventdata, handles)
% hObject    handle to edit_scales (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_scales as text
%        str2double(get(hObject,'String')) returns contents of edit_scales as a double


% --- Executes during object creation, after setting all properties.
function edit_scales_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_scales (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
