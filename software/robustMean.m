function [finalMean, stdSample, inlierIdx, outlierIdx] = robustMean(data,dim,k,fit,logical_idx)
%ROBUSTMEAN calculates mean and standard deviation discarding outliers
%
% SYNOPSIS [finalMean, stdSample, inlierIdx, outlierIdx] = robustMean(data,dim,k,fit)
%
% INPUT    data : input data
%          dim  : (opt) dimension along which the mean is taken {1}
%          k    : (opt) #of sigmas at which to place cut-off {3}
%          fit  : (opt) whether or not to use fitting to robustly estimate
%                  the mean from the data that includes outliers.
%                  0 (default): mean is approximated by median(data)
%                  1 : mean is approximated by
%                      fminsearch(@(x)(median(abs(data-x))),median(data))
%                      This option is only available for scalar data
%          logical_idx : (opt) logical indicating whether to return binary
%                        indicies for inlierIdx and outlierIdx
%                        false (default): actual idx, but slow due to find
%                        true           : logical idx, faster
%
%
% OUTPUT   finalMean : robust mean
%          stdSample : std of the data (divide by sqrt(n) to get std of the
%                      mean)
%          inlierIdx : index into data with the inliers
%          outlierIdx: index into data with the outliers
%
% REMARKS  NaN or Inf will be counted as neither in- nor outlier
%          The code is based on (linear)LeastMedianSquares. It could be changed to
%          include weights
%
% REFERENCES
%
% 1. PJ Rousseeuw and AM Leroy. Robust Regression and Outlier Detection.
%    New York: John Wiley & Sons, 1987. ISBN: 978-0-471-48855-2.
% 2. PJ Rousseeuw. Least median of squares regression. J. Am. Stat. Ass.
%    79, 871-880 (1984). DOI: 10.1080/01621459.1984.10477105.
% 3. G Danuser and M Stricker. Parameteric Model Fitting: From Inlier
%    Characterization to Outlier Detection. IEEE Trans Pattern Anal. Mach.
%    Intell. Vol 20, No. 2, March 1998. DOI: 10.1109/34.667884.
% 4. K Jaqaman and G Danuser. Linking Data to models: data regression.
%    Nature Rev Mol Cell Bio. 7, 813-819 (Nov 2006). DOI: 10.1038/nrm2030.
%
% c: jonas, 04/04
% Mark Kittisopikul, November 2015
%
% See also mad
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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


% test input

if isempty(data)
    error('Please supply non-empty data to robustMean')
end
if nargin<2 || isempty(dim)
    % make sure that the dimensinon is correct if there's a vector
    if any(size(data)==1) && ismatrix(data)
        dim = find(size(data)>1);
    else
        dim = 1;
    end
end
if nargin < 3 || isempty(k)
    % Cut-off is roughly at 3 sigma (References 1-3) by default
    k = 3;
end
% mkitti: was bug? only four parameters possible
% if nargin < 5 || isempty(fit)
if nargin < 4 || isempty(fit)
    fit = 0;
end
if nargin < 5 || isempty(logical_idx)
    logical_idx = false;
end

if fit == 1
    % check for vector
    if sum(size(data)>1)>1
        error('fitting is currently only supported for 1D data')
    end
end

insufficientData = true;
if(numel(data) >= 4)
    if(all(isfinite(data(1:4))) ...
    || all(isfinite(data(end-3:end))) ...
    || all(isfinite(data(((1:4)+floor(end/2)-2))))   )
        % Quickly check if first, last, or middle four are finite
        insufficientData = false;
    else
        finiteMap = isfinite(data);
        % Only need to find four
        finiteCount = numel(find(finiteMap,4));
        insufficientData = finiteCount < 4;
    end
end

% if sum(isfinite(data(:))) < 4
if insufficientData
    warning('ROBUSTMEAN:INSUFFICIENTDATA',...
        'Less than 4 data points!')
    finalMean = nanmean(data,dim);
    stdSample = NaN(size(finalMean));
    inlierIdx = find(isfinite(data));
    outlierIdx = [];
    return
end


%========================
% LEAST MEDIAN SQUARES
%========================

% Scale factor that relates Median absolute deviation (MAD) to standard deviation
% See mad(X,1)
% mad2stdSq=1.4826^2; %see same publications

% mad2stdSq = 1/norminv(3/4).^2
% mad2stdSq = 2.198109338317732

% backwards compatible constant
% sprintf('%0.9g',1.4826^2)
mad2stdSq = 2.19810276;

% calc median - reduce dimension dim to length 1
if fit
    % minimize the median deviation from the mean
    medianData = fminsearch(@(x)(median(abs(data-x))),median(data));
else
    medianData = nanmedian(data,dim);
end

% calculate statistics
% res2 = (data-repmat(medianData,blowUpDataSize)).^2;
res2 = bsxfun(@minus,data,medianData).^2;

medRes2 = max(nanmedian(res2,dim),eps);

%testvalue to calculate weights
% testValue=res2./repmat(mad2stdSq*medRes2,blowUpDataSize);
testValue = bsxfun(@rdivide,res2,mad2stdSq*medRes2);

% outlierIdx = testValue > k^2;
% Note: NaNs will always be false in comparison
inlierIdx = testValue <= k^2;
% Old outliers: 
% outlierIdx = find(testValue>k^2); % Does not include NaNs
% New outliers:
outlierIdx = ~inlierIdx; % Also includes NaNs

% Prior to Nov 2015, there used to be an if/else statement here depending if
% vector or higher dimensional input was given for data.
nInliers = sum(inlierIdx,dim);

% calculate std of the sample;
if nargout > 1

 %% Obsolete code block with bug   
    % put NaN wherever there are not enough data points to calculate a
    % standard deviation
%         goodIdx = sum(isfinite(res2),dim) > 4;
%mkitti, Oct 29 2015
% I believe the following commented out lines constitute a bug.
% goodIdx does not correctly index res2 in the expected manner.
% Therefore the second output of robustMean.m when supplied with a
% multidimensional input is invalid.
%         stdSample = NaN(size(goodIdx));
%         stdSample(goodIdx)=sqrt(nansum(res2(goodIdx),dim)./(nInliers(goodIdx)-4));

    %% outlierIdx should send NaN to zeros also so nansum not needed
    res2(outlierIdx) = 0;
    stdSample = sqrt(sum(res2,dim)./(nInliers-4));
    stdSample(nInliers <= 4) = NaN;
end

%====END LMS=========

%======
% MEAN
%======

data(outlierIdx) = 0;
finalMean = sum(data,dim)./nInliers;

if(nargout > 2 && ~logical_idx)
    % For backwards compatability only
    inlierIdx = find(inlierIdx);
end
if(nargout > 3)
    % For backwards compatability only, exclude NaNs
    % Above, NaNs are included as outliers
    % Next two line should be equivalent
%      outlierIdx = testValue > k^2;
    outlierIdx(outlierIdx) = ~isnan(testValue(outlierIdx));
    if(~logical_idx)
        outlierIdx = find(outlierIdx);
    end
end
