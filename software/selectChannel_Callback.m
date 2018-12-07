function selectChannel_Callback(hObject, eventdata, handles)

% Retrieve  channels properties
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
availableProps = get(handles.listbox_availableChannels, {'String','UserData','Value'});
selectedProps = get(handles.listbox_selectedChannels, {'String','UserData'});

% Find new elements and set them to the selected listbox
newID = availableProps{3}(~ismember(availableProps{1}(availableProps{3}),selectedProps{1}));
selectedChannels = horzcat(selectedProps{1}',availableProps{1}(newID)');
selectedData = horzcat(selectedProps{2}, availableProps{2}(newID));
set(handles.listbox_selectedChannels, 'String', selectedChannels, 'UserData', selectedData);
