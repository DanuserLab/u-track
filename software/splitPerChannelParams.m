function pThis = splitPerChannelParams(pAll,iChan)
%SPLITPERCHANNELPARAMS splits per-channel parameters to create single-channel parameter structure
%
% pThisChan = splitPerChannelParams(pAllChan,iChan)
%
% Set up parameter structure for detection on single specified channel,
% given input structure containing variable-format parameters for all
% channels. (designed to accepts output of e.g. prepPerChannelParams.m)
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

%Hunter Elliott
%6/2014


pThis = pAll;
for l = 1:numel(pAll.PerChannelParams)
    if isfield(pAll,pAll.PerChannelParams{l})%Check, some params may use defaults
        if iscell(pAll.(pAll.PerChannelParams{l}))
            pThis.(pAll.PerChannelParams{l}) = pAll.(pAll.PerChannelParams{l}){iChan};
        else
            pThis.(pAll.PerChannelParams{l}) = pAll.(pAll.PerChannelParams{l})(:,iChan);
            if all(isnan(pAll.(pAll.PerChannelParams{l})(:,iChan)))
                pThis.(pAll.PerChannelParams{l}) = [];
            end
        end        
    end
end