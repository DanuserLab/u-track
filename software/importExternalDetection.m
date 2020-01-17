function importExternalDetection(movieData,varargin)
% importExternalDetection imports external detection results into the movie infrastructure
%
% This function copies all the folders defined by the InputData field of
% the input parameters under a folder and registers this output in the
% external process.
%
%     importExternalDetection(movieData) runs the external detection
%     process on the input movie
%
%     importExternalDetection(movieData, paramsIn) additionally takes
%     the input parameters
%
%     paramsIn should be a structure with inputs for optional parameters.
%     The parameters should be stored as fields in the structure, with the
%     field names and possible values as described below
%
%   Possible Parameter Structure Field Names:
%       ('FieldName' -> possible values)
%
%       ('OutputDirectory' -> character string)
%       Optional. A character string specifying the directory to save the
%       detection results to. External detection results for different channels will be copied
%       under this directory.
%
%       ('InputData' -> Positive integer scalar or vector)
%       Optional. A nChanx1 cell array containing the paths of the folders
%       of the input detection results.
%
%
% This function is a wrapper function for ExternalDetectionProcess.m.
% This function is modified from importExternalSegmentation.m.
%
% Qiongjing (Jenny) Zou, Feb 2019
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

importExternalDetectData(movieData, 'ExternalDetectionProcess', varargin{:});
