function icon_ButtonDownFcn(hObject, eventdata)
% This function call up a help dialog box when user click any of the icons
% in all GUIs.
%
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
handles = guidata(hObject);

ud = get(hObject, 'UserData');

if ~isempty(ud) && isfield(ud, 'class')
%% Open original txt help file. 
% Help files are specified by variable "helpFile"
    
    helpFile=[ud.class '.pdf'];
    if exist(helpFile, 'file')
        if ispc || ismac
            open(helpFile)
        elseif isunix
            helpFilePath = which(helpFile);
            system(['evince ' helpFilePath ' &']);
        end
    else
        warndlg(['Cannot find Help file:', helpFile], 'modal');
        return
    end

else


%% Open GUI based help window
    
    userData = get(handles.figure1, 'UserData');
    % Help dialog from MovieData panel
    splitTag = regexp(get(get(hObject,'parent'), 'tag'), '_','split');

    % Pass handle to userData
    % If from package GUI, call pre-defined help dialog
    % if called from setting GUI, call user-defined help dialog 'msgboxGUI'

    if isfield(userData, 'crtProc')
        copyright = getLCCBCopyright();
        % Help dialog from setting panel
        if ~isempty(userData.crtProc)
            userData.helpFig = msgboxGUI('Text', sprintf([get(hObject,'UserData'), ...
                '\n', copyright ]),'Title',['Help - ' userData.crtProc.getName] );
        else
            userData.helpFig = msgboxGUI('Text', sprintf([get(hObject,'UserData'), ...
                '\n', copyright ]),'Title','Help');
        end


    elseif strcmp(splitTag{1}, 'axes') && length(splitTag) >1

            if strcmpi(splitTag{2}, 'help') % Help icon
                if length(splitTag) < 3
                    % Package help
                    userData.packageHelpFig = msgbox(sprintf(get(hObject,'UserData')), ...
                        ['Help - ' userData.crtPackage.getName], 'custom', get(hObject,'CData'), userData.colormap, 'replace');
                else
                    % Process help
                    procID = str2double(splitTag{3});
                    if ~isnan(procID)

                        procName = regexp(userData.crtPackage.getProcessClassNames{procID}, 'Process','split');
                        userData.processHelpFig(procID) = msgbox(sprintf(get(hObject,'UserData')), ...
                         ['Help - ' procName{1}], 'custom', get(hObject,'CData'), userData.colormap, 'replace');
                    end
                end

            else % Process status icon
                procID = str2double(splitTag{3});
                procName =userData.crtPackage.processes_{procID}.getName;
                userData.statusFig = msgbox(get(hObject,'UserData'), ...
                    'Status', 'custom', get(hObject,'CData'), ...
                    userData.colormap, 'replace');            
            end

    else
        userData.iconHelpFig = msgbox(get(hObject,'UserData'), ...
            'Help', 'custom', get(hObject,'CData'), userData.colormap, 'replace'); 
    end

    set(handles.figure1, 'UserData', userData);
end
