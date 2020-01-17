function costMat_OpeningFcn(hObject, eventdata, handles, parent, ID)
% This function fills various fields of the userData 
%       userData.mainFig - handle to the main figure
%       userData.handles_main - 'handles' of main figure
%       userData.ID - The ID of the cost matrix
%       userData.crtProc - current process
%       userData.parameters - current package
%       userData.questIconData - help icon image information
%       userData.colormap - color map information
%
% Sebastien Besson May 2011
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

% Check input
% The mainFig and procID should always be present
% procCOnstr and procName should only be present if the concrete process
% initation is delegated from an abstract class. Else the constructor will
% be directly read from the package constructor list.
ip = inputParser;
ip.addRequired('hObject',@ishandle);
ip.addRequired('eventdata',@(x) isstruct(x) || isempty(x));
ip.addRequired('handles',@isstruct);
ip.addOptional('parent',[],@ishandle);
ip.addOptional('ID',[],@isscalar);
ip.parse(hObject,eventdata,handles, parent, ID);


set(handles.text_copyright, 'String', getLCCBCopyright());

handles.output = hObject;
userData = get(handles.figure1, 'UserData');
if(isempty(userData)), userData = struct(); end;

% Get main figure handle and process id
userData.procID = ip.Results.ID;
userData.handles_main = guidata(parent);
userData.mainFig = userData.handles_main.figure1;
userData.userData_main = get(userData.handles_main.figure1, 'UserData');
userData.crtProc = userData.userData_main.crtProc;

props = get(parent, {'UserData', 'String', 'Value'});
userData.parameters = props{1}{userData.procID};
set(handles.text_name,'String',props{2}{props{3}});

% Get icon infomation
userData.questIconData = userData.userData_main.questIconData;
userData.colormap = userData.userData_main.colormap;

% ----------------------Set up help icon------------------------

% Set up help icon
set(hObject,'colormap',userData.colormap);
% Set up package help. Package icon is tagged as '0'
set(handles.figure1,'CurrentAxes',handles.axes_help);
Img = image(userData.questIconData); 
set(gca, 'XLim',get(Img,'XData'),'YLim',get(Img,'YData'),...
    'visible','off','YDir','reverse');

[~, filename] = fileparts(get(handles.figure1,'Filename'));
set(Img,'ButtonDownFcn',@icon_ButtonDownFcn,...
    'UserData', struct('class', filename));

% Update user data and GUI data
set(hObject, 'UserData', userData);