function compMatrices = discriminationMatrix(dataStructure,testStructure)
%DISCRIMINATIONMATRIX generates discriminationMatrices comparing means and distributions of lists of values in the input data structures
%
% SYNOPSIS compMatrices = discriminationMatrix(dataStructure)
%
% INPUT    dataStructure : n-by-1 structure containing the lists of values
%                          that are to be compared
%                          Fieldnames can be arbitrary
%          testStructure : (opt) 1-by-1 structure containing information on
%                          test to be used. Fieldnames have to be exactly
%                          the same as in dataStructure.
%                          The fields have to contain a vector of length
%                          two, specifying the tests to be used. The first
%                          entry specifies the test shown below the
%                          diagonal, the second entry specifies the test
%                          shown above the diagonal.
%                          1 - t-test for means (default below diag.)
%                          2 - Wilcoxon ranksum test for medians
%                          10 - K-S test for distributions
%                          11 - K-S test for distributions with subtracted
%                               mean (default for above diag.)
%                          12 - K-S test for distributions with subtracted
%                               median
%                          20 - permutation test for means
%                          21 - distribution test - calibrated K-S test
%                               with mean subtraction
%
% OUTPUT   compMatrices  : Structure with fieldnames equal to the
%                     fieldnames of the data structure, containing
%                     n-by-n matrices with p-values for
%                     - below the diagonal: T-test for equality of means
%                     - above the diagonal: Kolmogorov-Smirnov test for
%                         equality of the distributions, shifted by their
%                         means.
%                     Values in the diagonal are 1.01
%
% c: jonas, 12/04
% kathryn, 10/09 - added permutation and distribution tests
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

%=====================
% TEST INPUT
%=====================

defaultTest = [1 11];
% tests:
% If you add tests, please don't forget to reference them here AND in the
% help AND add the entry to the testList. Don't change the default.
% 1 - t-test for means (default below diag.)
% 2 - Wilcoxon ranksum test for medians
% 10 - K-S test for distributions
% 11 - K-S test for distributions with subtracted
%       mean (default for above diag.)
% 12 - K-S test for distributions with subtracted
%       median
% 20 - permutation test for means
% 21 - distribution test - calibrated KS with mean subtraction


testList = [1 2 10 11 12 20 21];

% nargin
if nargin == 0 || isempty(dataStructure) || ~isstruct(dataStructure)
    error('Please input a non-empty data structure in DISCRIMINATIONMATRIX')
end
if nargin < 2 || isempty(testStructure)
    useDefaultTest = 1;
else
    useDefaultTest = 0;
end


% length structure
nGroups = length(dataStructure);
if nGroups < 2
    error('Please try to compare at least two sets of data!')
end

% number of data sets to compare
dataNames = fieldnames(dataStructure);
nData = length(dataNames);

% loop through all the data. Get means, medians and fill testStruct
% (check whether empty, too)
groupMeans = zeros(nGroups,nData,2);
for iData = 1:nData

    for iGroup = 1:nGroups
        if ~isempty(dataStructure(iGroup).(dataNames{iData}))
% subtract nanMean/nanMedian to make sure that we don't lose all data for
% one NaN
            groupMeans(iGroup,iData,1) = ...
                nanmean(dataStructure(iGroup).(dataNames{iData}));
            groupMeans(iGroup,iData,2) = ...
                nanmedian(dataStructure(iGroup).(dataNames{iData}));


        else
            error('no empty entries are allowed in the dataStructure')
        end
    end

    % select tests
    if useDefaultTest
        testStructure.(dataNames{iData}) = defaultTest;
    else
        % first: make sure there is a field of the right name.
        % Then make sure two tests are given.
        % Finally, make sure they are on the testList
        if isfield(testStructure,dataNames{iData})
            if length(testStructure.(dataNames{iData}))==2 && ...
                    all(ismember(testStructure.(dataNames{iData}),testList))
                % all is good.
            else
                error('You either didn''t specify two tests, or you specified an invalid test')
            end
        else
            testStructure.(dataNames{iData}) = defaultTest;
        end
    end
end


%======================



%=====================
% TEST DATA
%=====================

% preassign output
rawMat = repmat(1.01,[nGroups,nGroups]);

% fore each data entry: Compare all the groups
for iData = 1:nData

    % preassign output
    compMatrices.(dataNames{iData}) = rawMat;

    for gi = 1:nGroups
        for gj = gi-1:-1:1

            % test 1: below diagonal
            try
            pValue = test(dataStructure(gi).(dataNames{iData}),...
                dataStructure(gj).(dataNames{iData}),...
                testStructure.(dataNames{iData})(1),...
                groupMeans([gi,gj],iData,:));
            compMatrices.(dataNames{iData})(gi,gj) = pValue;
            catch
                warning('data sets %i and %i could not be compared %s',...
                    gi,gj,lasterr)
            end

            % test 2: above diagonal
            try
            pValue = test(dataStructure(gi).(dataNames{iData}),...
                dataStructure(gj).(dataNames{iData}),...
                testStructure.(dataNames{iData})(2),...
                groupMeans([gi,gj],iData,:));
            compMatrices.(dataNames{iData})(gj,gi) = pValue;
            catch
                warning('data sets %i and %i could not be compared %s',...
                    gi,gj,lasterr)
            end

        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  SUBFUNCTION TEST    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pValue = test(data1, data2, whichTest, groupMeans)
% groupMeans is 2-by-1-by-2; (:,:,1) is the mean, (:,:,2) is the median

switch whichTest
    case 1
        % compare mean:  ttest. Below diagonal
        [dummy,pValue] = ...
            ttest2(data1, data2, 0.05,'both','unequal');
    case 2
        % compare median: Wilcoxon Mann Whitney ranksum test. Remove NaN
        % first
        pValue = ranksum(data1(~isnan(data1)), data2(~isnan(data2)));
    case 10
        % compare distributions: KS-test. No correction
        [dummy,pValue] = ...
            kstest2(data1, data2);
    case 11
        % compare distributions: KS-test. Subtract means
        [dummy,pValue] = ...
            kstest2(data1-groupMeans(1,1,1), data2-groupMeans(2,1,1));
    case 12
        % compare distributions: KS-test. Subtract medians
        [dummy,pValue] = ...
            kstest2(data1-groupMeans(1,1,2), data2-groupMeans(2,1,2));
    case 20
        % compare means: permutation test
        [~,pValue] = permTest(data1, data2, 0.05, 'both');
    case 21
        % compare distributions: calibrated KS-test
        pValue=distribTest(data1,data2); % here pValue is really the confidence value
end
