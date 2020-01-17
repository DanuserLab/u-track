function varargout = plusTipGroupAnalysisGUI(varargin)
% plusTipGroupAnalysisGUI M-file for plusTipGroupAnalysisGUI.fig
%      plusTipGroupAnalysisGUI, by itself, creates a new plusTipGroupAnalysisGUI or raises the existing
%      singleton*.
%
%      H = plusTipGroupAnalysisGUI returns the handle to a new plusTipGroupAnalysisGUI or the handle to
%      the existing singleton*.
%
%      plusTipGroupAnalysisGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in plusTipGroupAnalysisGUI.M with the given input arguments.
%
%      plusTipGroupAnalysisGUI('Property','Value',...) creates a new plusTipGroupAnalysisGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before plusTipGroupAnalysisGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to plusTipGroupAnalysisGUI_OpeningFcn via varargin.
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

% Edit the above text to modify the response to help plusTipGroupAnalysisGUI

% Last Modified by GUIDE v2.5 14-Apr-2014 13:02:10

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @plusTipGroupAnalysisGUI_OpeningFcn, ...
    'gui_OutputFcn',  @plusTipGroupAnalysisGUI_OutputFcn, ...
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


% --- Executes just before plusTipGroupAnalysisGUI is made visible.
function plusTipGroupAnalysisGUI_OpeningFcn(hObject, eventdata, handles, varargin)

set(handles.text_copyright, 'String', getLCCBCopyright())

userData = get(handles.figure1, 'UserData');
if(isempty(userData)), userData = struct(); end; % 2014b Hackaton-fix PR
% Get main figure handle and process id
ip = inputParser;
ip.addParamValue('mainFig', [], @isscalar)
ip.parse(varargin{:});
userData.mainFig = ip.Results.mainFig;
userData.handles_main = guidata(userData.mainFig);

% Get current package and process
userData_main = get(userData.mainFig, 'UserData');
userData.ML = userData_main.ML;

% Get icon infomation
userData.questIconData = userData_main.questIconData;
userData.colormap = userData_main.colormap;

% Set test values
testList = {'t-test of the means';'Wilcoxon ranksum test';'Kolmogorov-Smirnov test (K-S test)';...
    'Mean substracted K-S test';'Median substracted K-S test';...
    'Permutation t-test of the means';'Calibrated mean subtracted K-S test'};
testValues=[1 2 10 11 12 20 21];
test_handles = [handles.popupmenu_poolData_testID1...
    handles.popupmenu_poolData_testID2...
    handles.popupmenu_perCell_testID1...
    handles.popupmenu_perCell_testID2];
set(test_handles, 'String', testList, 'UserData', testValues);
set(handles.popupmenu_poolData_testID1, 'Value', 6);
set(handles.popupmenu_poolData_testID2, 'Value', 7);
set(handles.popupmenu_perCell_testID1, 'Value', 1);
set(handles.popupmenu_perCell_testID2, 'Value', 6);

handles.output = hObject;
set(hObject, 'UserData', userData);
guidata(hObject, handles);

% --- Outputs from this function are returned to the command line.
function varargout = plusTipGroupAnalysisGUI_OutputFcn(~, ~, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

userData = get(handles.figure1, 'UserData');
if(isempty(userData)), userData = struct(); end; % 2014b Hackaton-fix PR
if isempty(userData.ML)
    warndlg('At least one movie list is required to perform group analysis.',...
        'Input error', 'modal');
    close(handles.figure1);
    return
end

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

set(handles.figure1, 'UserData', userData);
guidata(hObject,handles);

% --- Executes on button press in pushbutton_run.
function pushbutton_run_Callback(hObject, eventdata, handles)


userData = get(handles.figure1,'UserData');
if(isempty(userData)), userData = struct(); end; % 2014b Hackaton-fix PR

% Load group data
remBegEnd = get(handles.checkbox_remBegEnd,'Value');
disp('Extracting group data');
userData.groupData = plusTipExtractGroupData(userData.ML, remBegEnd);

% Read common value for statistical tests
alpha =str2double(get(handles.edit_alpha, 'String'));
testValues = get(handles.popupmenu_poolData_testID1,'UserData');

% Run within group comparison
if get(handles.checkbox_doWtn,'Value')
    disp('Running within group comparison');
    for i = 1 : numel(userData.ML)
        fprintf(1, 'Movie %g/%g\n', i, numel(userData.ML));
        outputDir = fullfile(userData.ML(i).outputDirectory_,...
            'withinGroupComparison');
        plusTipWithinGroupComparison(userData.groupData, i, outputDir, 1);
    end
end

% Run per-cell group analysis
if numel(userData.ML) > 1 && get(handles.checkbox_poolData,'Value')
    disp('Running pooled data analysis');
    outputDir = fullfile(userData.ML(1).outputDirectory_, 'pooledData');
    plusTipPoolGroupData(userData.groupData, outputDir, 0, 1);
    
    % Perform statistical tests if more than one list is passed
    testID1 = testValues(get(handles.popupmenu_poolData_testID1,'Value'));
    testID2 = testValues(get(handles.popupmenu_poolData_testID2,'Value'));
    plusTipTestDistrib(userData.groupData, outputDir,...
        alpha, testID1, testID2);
end

% Run pooled group analysis
if get(handles.checkbox_perCell,'Value')
    disp('Running per cell analysis');
    outputDir = fullfile(userData.ML(1).outputDirectory_, 'perCell');
    testID1 = testValues(get(handles.popupmenu_perCell_testID1,'Value'));
    testID2 = testValues(get(handles.popupmenu_perCell_testID2,'Value'));
    plusTipGetHits(userData.groupData, outputDir, alpha, testID1, testID2);
end

arrayfun(@save, userData.ML)

% --- Executes on button press in checkbox_poolData.
function checkbox_poolData_Callback(hObject, eventdata, handles)

if get(hObject, 'Value')
    set(get(handles.uipanel_poolData, 'Children'), 'Enable', 'on');
else
    set(get(handles.uipanel_poolData, 'Children'), 'Enable', 'off');
end


% --- Executes on button press in checkbox_perCell.
function checkbox_perCell_Callback(hObject, eventdata, handles)

if get(hObject, 'Value')
    set(get(handles.uipanel_perCell, 'Children'), 'Enable', 'on');
else
    set(get(handles.uipanel_perCell, 'Children'), 'Enable', 'off');
end
