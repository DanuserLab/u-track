function varargout = pointSourceDetectionProcessGUI(varargin)
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

% Last Modified by GUIDE v2.5 01-Oct-2013 19:23:29

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
userData=get(handles.figure1,'UserData');
funParams = userData.crtProc.funParams_;
%Remove the output directory as we don't want to replicate it to other
%movies if the "apply to all movies" box is checked. Ideally we would
%explicitly only replicate the parameters we set in this GUI but this is a
%quick fix. - HLE
if isfield(funParams,'OutputDirectory')
    funParams = rmfield(funParams,'OutputDirectory');
end

set(handles.popupmenu_CurrentChannel,'UserData',funParams);

iChan = get(handles.popupmenu_CurrentChannel,'Value');
if isempty(iChan)
    iChan = 1;
    set(handles.popupmenu_CurrentChannel,'Value',1);
end

% %Specify parameters
% userData.numParams = {'filterSigma', 'alpha','Mode','MaxMixtures','RedundancyRadius'};
% %Set the intersection checkbox 
% userData.checkParams ={'FitMixtures'};
% 

% for i =1 : numel(userData.numParams)
%     paramName = userData.numParams{i};        
%     if any(strcmp(paramName,funParams.PerChannelParams))
%         parVal = funParams.(paramName)(iChan);
%     else
%         parVal = funParams.(paramName);
%     end
%     set(handles.(['edit_' paramName]), 'String',parVal);
% end
% 
% 
% 
% set(handles.edit_UseIntersection,,'Value',funParams.UseIntersection)
% 
% for i =1 : numel(userData.checkParams)
%     paramName = userData.checkParams{i};        
%     if any(strcmp(paramName,funParams.PerChannelParams))
%         parVal = funParams.(paramName)(iChan);
%     else
%         parVal = funParams.(paramName);
%     end
%     set(handles.(['edit_' paramName]), 'Value',parVal);
% end

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

% Update channels listboxes depending on the selected process
popupmenu_SegProcessIndex_Callback(hObject, eventdata, handles)

%Update channel parameter selection dropdown
popupmenu_CurrentChannel_Callback(hObject, eventdata, handles)


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
userData=get(handles.figure1,'UserData');
funParams = get(handles.popupmenu_CurrentChannel,'UserData');

selChan = 1:numel(userData.MD.channels_);%For now just let them set parameters for all channels.
%selChan = get(handles.listbox_selectedChannels,'UserData');
chanStr = arrayfun(@(x)(['Channel ' num2str(x)]),selChan,'Unif',0);
set(handles.popupmenu_CurrentChannel,'String',chanStr);
iChan = get(handles.popupmenu_CurrentChannel,'Value');
%set(handles.popupmenu_CurrentChannel,'UserData',iChan);



% Set-up parameters
for i =1 : numel(funParams.PerChannelParams)
    paramName = funParams.PerChannelParams{i};
    parVal = funParams.(paramName)(iChan);
    if ~islogical(funParams.(paramName))
        set(handles.(['edit_' paramName]), 'String',parVal);
    else
        set(handles.(['edit_' paramName]), 'Value',parVal);
    end
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

%userData=get(handles.figure1,'UserData');
funParams = get(handles.popupmenu_CurrentChannel,'UserData');

for i =1 : numel(funParams.PerChannelParams)
    paramName = funParams.PerChannelParams{i};
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
        
end

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
