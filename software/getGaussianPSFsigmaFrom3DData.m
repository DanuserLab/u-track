%[sigma] = getGaussianPSFsigmaFrom3DData(volList) returns the s.d. of the Gaussian PSF estimated from the input data.
% The estimation is performed by running pointSourceDetection3D() with 'sigma' as a free parameter.
%
% Inputs:
%    volList : single vol, cell array of vols, or cell array of vol path strings
%
% Options:
%    'Display' : {true}|false displays the distribution of 'sigma' values
%
% Output:
%        sigma : standard deviation of the Gaussian PSF estimated from point sources in input data
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

% Philippe Roudot after Francois Aguet, October 2014

function sigma = getGaussianPSFsigmaFrom3DData(volList, varargin)
ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('volList');
ip.addOptional('frameRange', []);
ip.addParamValue('Display', true, @islogical);
ip.parse(volList, varargin{:});

if ~iscell(volList)
    volList = {volList};
end
frameRange = ip.Results.frameRange;

if isempty(frameRange)
    nd = numel(volList);
else
    nd = numel(frameRange);
end
svect = cell(1,nd);
parfor i = 1:nd
    if isempty(frameRange)
        if ischar(volList{i})
            vol = double(imread(volList{i}));
        else
            vol = volList{i};
        end
    else
        vol = double(readtiff(volList{1}, frameRange(i)));
    end
    % First pass with fixed sigma
    pstruct = pointSourceDetection3D(vol, [2 1], 'Mode', 'xyac');
    if ~isempty(pstruct)
        pstruct = fitGaussians3D(vol, [pstruct.x' pstruct.y' pstruct.z'],pstruct.A', [2 1],pstruct.c','xyzAsrc');
        isPSF = ~[pstruct.hval_AD] & [pstruct.pval_Ar] < 0.05;
        svect{i} = pstruct.s(:,~isnan(pstruct.s(1,:)) & isPSF);
    end
end
svect = [svect{:}];

opts = statset('maxIter', 200);
try
    w = warning('off', 'stats:gmdistribution:FailedToConverge');
    nTestedComponent=5;
    obj = cell(1,length(nTestedComponent));
    for nIdx = 1:length(nTestedComponent)
        obj{nIdx} = gmdistribution.fit(svect', nTestedComponent(nIdx), 'Options', opts);
    end
    [~,bestNComponent] = min(cellfun(@(i) i.BIC, obj));
    bestObj = obj{bestNComponent};
    %svec = sqrt(squeeze(bestObj.Sigma));
    mixing = bestObj.PComponents;
    [~,maxMixIdx] = max(mixing);
    %[~,idx] = max(amp./(sqrt(2*pi)*svec));
    sigma = bestObj.mu(maxMixIdx,:);

    warning(w);

    if ip.Results.Display
        figure()
        plot(cellfun(@(i) i.BIC, obj));

        figure()
        hold on
        scatter(svect(1,:),svect(2,:),10,'.');
        h = ezcontour(@(x,y)pdf(bestObj,[x y]),[0 5],[0 5]);
        hold off
        xlabel('Lateral standard deviation of 3D Gaussian model (pixels)');
        ylabel('Axial standard deviation of 3D Gaussian model (pixels)');
    end
catch
    fprintf('Could not determine distribution, potentially due to insufficient samples.');
    sigma = mean(svect,2)';
end
