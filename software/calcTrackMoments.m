function trackMoments = calcTrackMoments(coordinates,standardDevs,...
    momentOrders,maxLag)
%CALCTRACKMOMENTS calculates the order^th moments of a track (MSD when order = 2)
%
%SYNPOSIS: trackMoments = calcTrackMoments(coordinates,standardDevs,...
%    momentOrders,maxLag)
%
%INPUT  coordinates : nTP-by-nD array of coordinates.
%                     (nTP = number of time points, nD = dimensionality).
%                     Missing coordinates should be indicated by NaN.
%       standardDevs: nTPD-by-nD array of coordinate standard deviations
%                     Optional. Default: 0 (perfect coordinates).
%       momentOrders: Vector of moment orders to be calculated.
%                     Optional. Default: 1 to 5.
%       maxLag      : Maximum lag for moment calculation.
%                     Optional. Default: nTP/4.
%
%OUTPUT trackMoments: Structure array with # entries = # moment orders.
%                     Contains the fields:
%           .momentOrder: Order of moment calculated.
%           .momentValues: column 1 - Moment values from lag 1 to maxLag.
%                          column 2 - Moment standard deviations.
%
%REMARKS The pth moment at lag dt is calculated as <|x(t+dt)-x(t)|^p> where
%        <...> denotes the average over all available |x(t+dt)-x(t)| from
%        a track (i.e. time average).
%
%Khuloud Jaqaman, March 2008
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

%% input
if nargin < 1 || isempty(coordinates)
    disp('--calcTrackMoments: Please input at least track coordinates!');
    return
end

%get number of time points and dimensionality
[numTimePoints,probDim] = size(coordinates);

%assign standard deviations if not input
if nargin < 2 || isempty(standardDevs)
    standardDevs = zeros(numTimePoints,probDim) + eps;
end

%assign momentOrder if not input
if nargin < 3 || isempty(momentOrders)
    momentOrders = 1 : 5;
end
numOrders = length(momentOrders);

%assign maxLag if not input
if nargin < 4 || isempty(maxLag) || maxLag > floor(numTimePoints/4)
    maxLag = floor(numTimePoints/4);
end

%% moments calculation

%reserve memory for moment values
momentValues = NaN(maxLag,numOrders);
momentStds = momentValues;

%go over all lags
for iLag = 1 : maxLag
   
    %calculate displacements over the current lag
    displacements = coordinates(iLag+1:end,:) - coordinates(1:end-iLag,:);
    
    %calculate displacement variances
    dispVar = standardDevs(iLag+1:end,:).^2 + standardDevs(1:end-iLag,:).^2;
    
    %get displacement magnitudes
    dispMag = sqrt(sum(displacements.^2,2));
    
    %get standard deviation of displacement magnitude
    dispMagStd = sqrt( sum(displacements.^2 .* dispVar,2) ) ./ dispMag;
    dispMagStd(dispMagStd==0) = eps;
    
    %keep only non-NaNs
    goodEntries = find(~isnan(dispMag));
    dispMag = dispMag(goodEntries);
    dispMagStd = dispMagStd(goodEntries);
    
    %     %remove outlier values
    %     [outlierIdx,inlierIdx] = detectOutliers(dispMag,20);
    %     dispMag = dispMag(inlierIdx);
    %     dispMagStd = dispMagStd(inlierIdx);
    
    %proceed if there are at least 5 values to calculate the average
    if length(dispMag) >= 5
        %         if length(dispMag) >= 10
    
        %evaluate the moments at this lag
        for iOrder = 1 : numOrders

            %get the current order
            orderVal = momentOrders(iOrder);

            %calculate the moment
            [momentValues(iLag,iOrder),momentStds(iLag,iOrder)] = ...
                weightedStats(dispMag.^orderVal,dispMagStd);
            
        end %(for iOrder = 1 : numOrders)
        
    end %(if length(dispMag) >= 10)
    
end

%% output

%store output for all moments
trackMoments = repmat(struct('momentOrder',[],'momentValues',[]),1,numOrders);
for iOrder = 1 : numOrders
    trackMoments(iOrder).momentOrder = momentOrders(iOrder);
    trackMoments(iOrder).momentValues = [momentValues(:,iOrder) momentStds(:,iOrder)];
end

%% ~~~ the end ~~~
