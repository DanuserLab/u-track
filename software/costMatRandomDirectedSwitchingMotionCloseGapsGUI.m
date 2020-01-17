function varargout = costMatRandomDirectedSwitchingMotionCloseGapsGUI(varargin)
% COSTMATRANDOMDIRECTEDSWITCHINGMOTIONCLOSEGAPSGUI M-file for costMatRandomDirectedSwitchingMotionCloseGapsGUI.fig
%      COSTMATRANDOMDIRECTEDSWITCHINGMOTIONCLOSEGAPSGUI, by itself, creates a new COSTMATRANDOMDIRECTEDSWITCHINGMOTIONCLOSEGAPSGUI or raises the existing
%      singleton*.
%
%      H = COSTMATRANDOMDIRECTEDSWITCHINGMOTIONCLOSEGAPSGUI returns the handle to a new COSTMATRANDOMDIRECTEDSWITCHINGMOTIONCLOSEGAPSGUI or the handle to
%      the existing singleton*.
%
%      COSTMATRANDOMDIRECTEDSWITCHINGMOTIONCLOSEGAPSGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in COSTMATRANDOMDIRECTEDSWITCHINGMOTIONCLOSEGAPSGUI.M with the given input arguments.
%
%      COSTMATRANDOMDIRECTEDSWITCHINGMOTIONCLOSEGAPSGUI('Property','Value',...) creates a new COSTMATRANDOMDIRECTEDSWITCHINGMOTIONCLOSEGAPSGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before costMatRandomDirectedSwitchingMotionCloseGapsGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to costMatRandomDirectedSwitchingMotionCloseGapsGUI_OpeningFcn via varargin.
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

% Edit the above text to modify the response to help costMatRandomDirectedSwitchingMotionCloseGapsGUI

% Last Modified by GUIDE v2.5 19-Dec-2016 15:36:25

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @costMatRandomDirectedSwitchingMotionCloseGapsGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @costMatRandomDirectedSwitchingMotionCloseGapsGUI_OutputFcn, ...
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


% --- Executes just before costMatRandomDirectedSwitchingMotionCloseGapsGUI is made visible.
function costMatRandomDirectedSwitchingMotionCloseGapsGUI_OpeningFcn(hObject, eventdata, handles, varargin)


costMat_OpeningFcn(hObject, eventdata, handles, varargin{:})
userData = get(handles.figure1, 'UserData');
if(isempty(userData)), userData = struct(); end;
parameters = userData.parameters;

% Brownian motion parameters
set(handles.edit_lower, 'String', num2str(parameters.minSearchRadius))
set(handles.edit_upper, 'String', num2str(parameters.maxSearchRadius))
set(handles.edit_brownStdMult, 'String', num2str(parameters.brownStdMult(1)))
set(handles.checkbox_useLocalDensity, 'Value', parameters.useLocalDensity)
set(handles.edit_nnWindow, 'String', num2str(parameters.nnWindow))
set(handles.edit_before, 'String', num2str(parameters.brownScaling(1)))
set(handles.edit_after, 'String', num2str(parameters.brownScaling(2)))
set(handles.edit_gapLengthTransitionB, 'String', num2str(parameters.timeReachConfB-1))
set(handles.edit_gapPenalty, 'String', num2str(parameters.gapPenalty));

% Directed Motion parameters
if parameters.linearMotion
    set(get(handles.uipanel_linearMotion,'Children'),'Enable','on');
    
    set(handles.edit_lenForClassify, 'String', num2str(parameters.lenForClassify))
    set(handles.edit_linStdMult, 'String', num2str(parameters.linStdMult(1)))
    set(handles.edit_before_2, 'String', num2str(parameters.linScaling(1)))
    set(handles.edit_after_2, 'String', num2str(parameters.linScaling(2)))
    set(handles.edit_gapLengthTransitionL, 'String', num2str(parameters.timeReachConfL-1))  
    set(handles.edit_maxAngleVV, 'String', num2str(parameters.maxAngleVV))    
else
    set(get(handles.uipanel_linearMotion,'Children'),'Enable','off');
end
set(handles.checkbox_linearMotion, 'Value', logical(parameters.linearMotion),'Enable','off');
set(handles.checkbox_immediateDirectionReversal, 'Value', parameters.linearMotion==2,'Enable','off');

% Merging/splitting parameters
mergeSplitComponents = findobj(handles.uipanel_mergeSplit,'-not','Type','uipanel');
if get(userData.handles_main.checkbox_merging,'Value') || ...
    get(userData.handles_main.checkbox_splitting,'Value')
    set(mergeSplitComponents,'Enable','on');
    if isempty(parameters.ampRatioLimit) || ...
            (length(parameters.ampRatioLimit) ==1 && parameters.ampRatioLimit == 0)
        set(handles.checkbox_ampRatioLimit, 'Value', 0)
        set(get(handles.uipanel_ampRatioLimit,'Children'),'Enable','off');
    else
        set(get(handles.uipanel_ampRatioLimit,'Children'),'Enable','on');
        set(handles.edit_min, 'String', num2str(parameters.ampRatioLimit(1)))
        set(handles.edit_max, 'String', num2str(parameters.ampRatioLimit(2)))
    end
    set(handles.edit_resLimit, 'String', num2str(parameters.resLimit))
    if ~isfield(parameters, 'gapExcludeMS')
        parameters.gapExcludeMS = 0;
    end
    set(handles.checkbox10_gapExcludeMS, 'Value', parameters.gapExcludeMS);
else
    set(mergeSplitComponents,'Enable','off');
    
end

% Birth and Death Strategy
if ~isfield(parameters, 'strategyBD')
    parameters.strategyBD = 0;
end

set(handles.edit_manual_percentile,'String','');

switch parameters.strategyBD
    case -1
        set(handles.edit_manual_percentile,'String','');
        set(handles.edit_manual_percentile,'Enable','off');
        set(handles.Derived_BDButton,'Value', 1.0);
        set(handles.perctile_Button, 'Value', 0.0);
        set(handles.perctAutomatic, 'Value', 0.0);
    
    case 0
        set(handles.edit_manual_percentile,'String','');
        set(handles.edit_manual_percentile,'Enable','off');
        set(handles.Derived_BDButton,'Value', 0.0);
        set(handles.perctile_Button, 'Value', 0.0);
        set(handles.perctAutomatic, 'Value', 1.0);
    
    otherwise
        set(handles.edit_manual_percentile,'Enable','on');
        set(handles.edit_manual_percentile, 'String', num2str(parameters.strategyBD));
        set(handles.Derived_BDButton,'Value', 0.0);
        set(handles.perctile_Button, 'Value', 1.0);
        set(handles.perctAutomatic, 'Value', 0.0);
end
    

% Update handles structure
handles.output = hObject;
guidata(hObject, handles);

% UIWAIT makes costMatRandomDirectedSwitchingMotionCloseGapsGUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = costMatRandomDirectedSwitchingMotionCloseGapsGUI_OutputFcn(hObject, eventdata, handles) 
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
if(isempty(userData)), userData = struct(); end;
parameters = userData.parameters;

% Brownian motion parameters
minSearchRadius = str2double(get(handles.edit_lower, 'String'));
maxSearchRadius = str2double(get(handles.edit_upper, 'String'));
brownStdMult = str2double(get(handles.edit_brownStdMult, 'String'));
nnWindow = str2double(get(handles.edit_nnWindow, 'String'));
brownScaling_1 = str2double(get(handles.edit_before, 'String'));
brownScaling_2 = str2double(get(handles.edit_after, 'String'));
gapLengthTransitionB = str2double(get(handles.edit_gapLengthTransitionB, 'String'));

isPosScalar = @(x) isscalar(x) &&~isnan(x) && x>=0;   
% lower
if  ~isPosScalar(minSearchRadius)
    errordlg('Please provide a valid value to parameter "Lower Bound".','Error','modal')
    return
end

% Upper
if ~isPosScalar(maxSearchRadius)
    errordlg('Please provide a valid value to parameter "Upper Bound".','Error','modal')
    return
elseif maxSearchRadius < minSearchRadius
    errordlg('"Upper Bound" should be larger than "Lower Bound".','Error','modal')
    return
end

% brownStdMult
if ~isPosScalar(brownStdMult)
    errordlg('Please provide a valid value to parameter "Multiplication Factor for Search Radius Calculation".','Error','modal')
    return
end
brownStdMult = brownStdMult*ones(userData.crtProc.funParams_.gapCloseParam.timeWindow,1);

% nnWindow
if ~isposint(nnWindow)
    errordlg('Please provide a valid value to parameter "Number of Frames for Nearest Neighbor Distance Calculation".','Error','modal')
    return
end

% brownScaling
if ~isPosScalar(brownScaling_1)
    errordlg('Please provide a valid value to parameter "Scaling Power in Fast Expansion Phase".','Error','modal')
    returnparameters.gapExcludeMS = get(handles.checkbox10_gapExcludeMS, 'Value');
end

% brownScaling
if ~isPosScalar(brownScaling_2)
    errordlg('Please provide a valid value to parameter "Scaling Power in Slow Expansion Phase".','Error','modal')
    return
end

brownScaling = [brownScaling_1 brownScaling_2];

% gapLengthTransitionB
if ~isPosScalar(gapLengthTransitionB)
    errordlg('Please provide a valid value to parameter "Gap length to transition from Fast to Slow Expansion".','Error','modal')
    return
end

% gapPenalty
gapPenalty = get(handles.edit_gapPenalty, 'String');
if isempty(gapPenalty)
    gapPenalty = [];
else
    gapPenalty = str2double(gapPenalty);
    if ~isPosScalar(gapPenalty)
        errordlg('Please provide a valid value to parameter "Time to Reach Confinement".','Error','modal')
        return
    end
end

parameters.minSearchRadius = minSearchRadius;
parameters.maxSearchRadius = maxSearchRadius;
parameters.brownStdMult = brownStdMult;
parameters.useLocalDensity = get(handles.checkbox_useLocalDensity, 'Value');
parameters.nnWindow = nnWindow;
parameters.brownScaling = brownScaling;
parameters.timeReachConfB = gapLengthTransitionB+1;
parameters.gapPenalty = gapPenalty;

% Merging/splitting parameters
if get(userData.handles_main.checkbox_merging,'Value') || ...
        get(userData.handles_main.checkbox_splitting,'Value')
    
    if ~get(handles.checkbox_ampRatioLimit, 'Value')
        ampRatioLimit = [];
    else
        ampRatioLimit_1 = str2double(get(handles.edit_min, 'String'));
        ampRatioLimit_2 = str2double(get(handles.edit_max, 'String'));
        
        % ampRatioLimit_1
        if ~isPosScalar(ampRatioLimit_1)
            errordlg('Please provide a valid value to parameter "Min Allowed".','Error','modal')
            return
        end
        
        % ampRatioLimit_2
        if ~isPosScalar(ampRatioLimit_2)
            errordlg('Please provide a valid value to parameter "Max Allowed".','Error','modal')
            return
        end
        
        if ampRatioLimit_2 <= ampRatioLimit_1
            errordlg('"Max Allowed" should be larger than "Min Allowed".','Error','modal')
            return
        end
        
        ampRatioLimit = [ampRatioLimit_1 ampRatioLimit_2];
    end
    
    % resLimit
    resLimit = get(handles.edit_resLimit, 'String');
    if isempty( resLimit )
        resLimit = [];
    else
        resLimit = str2double(resLimit);
        if ~isPosScalar(resLimit)
            errordlg('Please provide a valid value to parameter "Time to Reach Confinement".','Error','modal')
            return
        end
    end
    
    parameters.gapExcludeMS = get(handles.checkbox10_gapExcludeMS, 'Value');
    parameters.ampRatioLimit = ampRatioLimit;
    parameters.resLimit = resLimit;
end


% Birth and Death Strategy
%%% ADD error checking
BDstateSelection = get(get(handles.uibutton_BDstrat, 'SelectedObject'),'Tag');
%  handles.uibutton_BDstrat.SelectedObject.String % After 2014b
switch BDstateSelection
    case 'Derived_BDButton'
        parameters.strategyBD = -1;
    case 'perctAutomatic'
        parameters.strategyBD = 0;
    case 'perctile_Button'
        prctBD = get(handles.edit_manual_percentile, 'String');
        if isempty(prctBD)
            errordlg('Please provide a valid value to parameter "Percentile of the gap closing and merging and splitting cost distribution to use for calculating the birth and death".','Error','modal');
            return
        else
            prctBD = str2double(prctBD);
            if ~isPosScalar(prctBD) || (prctBD >= 100) || (prctBD < 0)
                errordlg('Please provide a valid value to parameter "Percentile of the gap closing and merging and splitting cost distribution to use for calculating the birth and death".','Error','modal');
                return
            else
                parameters.strategyBD = prctBD;
            end
        end

    otherwise
        error('bad tag....');
end


% Linear motion parameters
if parameters.linearMotion
    
    lenForClassify = str2double(get(handles.edit_lenForClassify, 'String'));
    linStdMult = str2double(get(handles.edit_linStdMult, 'String'));
    linScaling_1 = str2double(get(handles.edit_before_2, 'String'));
    linScaling_2 = str2double(get(handles.edit_after_2, 'String'));
    gapLengthTransitionL = str2double(get(handles.edit_gapLengthTransitionL, 'String'));
    maxAngleVV = str2double(get(handles.edit_maxAngleVV, 'String'));
    
    % lenForClassify
    if ~isPosScalar(lenForClassify)
        errordlg('Please provide a valid value to parameter "Minimum Track Segment Length to Classify it as Linear or Random".','Error','modal')
        return
    end
    
    % linStdMult
    if ~isPosScalar(linStdMult)
        errordlg('Please provide a valid value to parameter "Multiplication Factor for Linear Search Radius Calculation".','Error','modal')
        return
    end
    linStdMult = linStdMult*ones(userData.crtProc.funParams_.gapCloseParam.timeWindow,1);
    
    % linScaling_1
    if ~isPosScalar(linScaling_1)
        errordlg('Please provide a valid value to parameter "Scaling Power in Fast Expansion Phase".','Error','modal')
        return
    end
    
    % linScaling_1
    if ~isPosScalar(linScaling_2)
        errordlg('Please provide a valid value to parameter "Scaling Power in Slow Expansion Phase".','Error','modal')
        return
    end
    
    linScaling = [linScaling_1 linScaling_2];
    
    % gapLengthTransitionL
    if ~isPosScalar(gapLengthTransitionL)
        errordlg('Please provide a valid value to parameter "Gap length to transition from Fast to Slow Expansion".','Error','modal')
        return
    end
    
    % maxAngleVV
    if ~isPosScalar(maxAngleVV)
        errordlg('Please provide a valid value to parameter "Maximum Angle Between Linear Track Segments".','Error','modal')
        return
    end
    
    parameters.lenForClassify = lenForClassify;
    parameters.linStdMult = linStdMult;
    parameters.linScaling = linScaling;
    parameters.timeReachConfL = gapLengthTransitionL+1;
    parameters.maxAngleVV = maxAngleVV;
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

% --- Executes on button press in checkbox_ampRatioLimit.
function checkbox_ampRatioLimit_Callback(hObject, eventdata, handles)

if get(hObject, 'Value')
    set(get(handles.uipanel_ampRatioLimit,'Children'),'Enable','on');
else
    set(get(handles.uipanel_ampRatioLimit,'Children'),'Enable','off');  
end


% --- Executes on button press in checkbox10_gapExcludeMS.
function checkbox10_gapExcludeMS_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox10_gapExcludeMS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
disp('gapExcludeMS check toggled');
% Hint: get(hObject,'Value') returns toggle state of checkbox10_gapExcludeMS


% --- Executes on button press in perctile_Button.
function perctile_Button_Callback(hObject, eventdata, handles)
% hObject    handle to perctile_Button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.edit_manual_percentile, 'Enable', 'on');
% set(handles.edit_manual_percentile, 'String', num2str(parameters.strategyBD));

% Hint: get(hObject,'Value') returns toggle state of perctile_Button
