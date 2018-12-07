function out = catStruct(dim,structName,separator,string)
%CATSTRUCT catenates values from multilevel structures
%
% SYNOPSIS out = catStruct(dim,structName)
%
% INPUT   dim        : dimension along which you want the struct to be catenated
%         structName : (string) full tree of the structure down to the field to be
%                      catenated, e.g. 'struct.fieldname.subfieldname'. 
%                      USE RANGES ONLY FOR THE LAST LEVEL
%         separator  : (opt) number (NaN/Inf possible) or string (if
%                      string==1) that should separate the values in the
%                      catenated list. LENGTH HAS TO BE 1!
%         string     : (opt) 1 if the contents of the fields are strings (def: 0)
%
% OUTPUT  out        : array containtin the catenated values
%
% SAMPLE CALLS: speed = catStruct(1,wt30,'wt30.individualStatistics.summary.antipolewardSpeed');
%               nums  = catStruct(1,wt30,'wt30.convergenceClustering.apNum(:,1:2)', NaN); 
%               for NaN-separated list
%
% c: 2/04 jonas
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
% test input
%=====================

% nargin
if nargin < 2 || isempty(dim) || isempty(structName)
    error('CATSTRUCT needs two non-empty input arguments!')
end


% test string-opt
if nargin > 3 && ~isempty(string)
    % we have a string argument
    if ~(ismember(string,[1 0]))
        error('string has to be either 0 or 1')
    end
    if ~ismember(dim,[1 2])
        error('strings can only be catenated along dims 1 and 2 right now')
    end
else
    % set default
    string = 0;
end

% test separator
if nargin > 2 && ~isempty(separator)
    if string
        if ~ischar(separator)
            error('separator has to be string if ''string'' is 1')
        elseif length(separator)>1
            error('maxLength of separator is 1')
        else
            separator = ['''' separator ''''];
            isseparator = 1;
        end
    else
        if ~isnumeric(separator)
            error('separator has to be numeric if ''string'' is 1')
        elseif length(separator)>1
            error('maxLength of separator is 1')
        else
            separator = num2str(separator);
            isseparator = 1;
        end
    end
else
    separator = [];
    isseparator = 0;
end



% structName
if ~isstr(structName)
    error('please input the structure as a string')
end

if ~strfind(structName,'.')
    error('no ''.'' found in structName. Make sure you input a structure')
end
% we allow range for lowest level only - therefore check whether we have to
% chop off a bit of the variable name
openParenthesis= strfind(structName,'(');
if ~isempty(openParenthesis)
    % if there is only one opening parenthesis and the closing one is at
    % the very end of the input, we are happy, else we let the user run
    % into the next check
    if length(openParenthesis)==1 & strfind(structName(end),')')
        rangeString = structName(openParenthesis:end);
        structName(openParenthesis:end)=[];
    else
        rangeString = '';
    end
else
    rangeString = '';
end

if strfind(structName,'(') | strfind(structName,';') | strfind(structName,':') | strfind(structName,',') | strfind(structName,')')
    error('CATSTRUCT only supports ranges for the very last level')
end

% dim is tested while catenating

%====================
% end test input
%====================



%=================================
% read structure info and catenate
%=================================

% parse structure name to find number of levels & number of entries. if
% only one level, we can use good old cat, otherwise, we will build a
% nested loop (horribile dictu!) to fill in an array

% dimension string
dimString = num2str(dim);

% '.' marks boundaries between levels
levelBreaks = strfind(structName,'.');

% a level is everything above the last field name
numberOfLevels = length(levelBreaks);

% load structure into mFile; make sure it gets the right name!!
topLevelName = structName(1:levelBreaks(1)-1);
eval([topLevelName '= evalin(''caller'',''' topLevelName ''');'])

% make sure the top level is not called 'out'!
if strmatch(topLevelName,'out')
    topLevelName = 'input';
    input = out;
    clear out
end

% check whether out is empty
emptyInput = eval(sprintf('isempty(%s);',topLevelName));
if emptyInput
    out = [];
    return
end


% decide whether we can go easy or not
if numberOfLevels == 1 & isempty(separator) & isempty(rangeString)
    % string check
    if string
        if dim == 2
            evalString = ['out = strcat(' structName ');'];
        else
            evalString = ['out = strvcat(' structName ');'];
        end
    else
        % normal cat
        evalString = ['out = cat(' dimString ',' structName ');'];
    end
else
    
    
    
    % adjust levelBreaks so that we can include the last filename
    levelBreaks = [levelBreaks length(structName)+1];
    
    % init build nested loop
    startString = [];
    middleString = [];
    endString   = []; 
    currentLevelName = topLevelName;
    levelSize = 1; 
    
    breakString = [];
    
    % build nested loop
    for nLevel = 1:numberOfLevels
        
        % get next level name
        iLevel = ['i' num2str(nLevel)];
        nextLevelName = [currentLevelName '(' iLevel ').' structName(levelBreaks(nLevel)+1:levelBreaks(nLevel+1)-1)];
        
        % get number of elements of current level
        eval([iLevel '=1;']);
        levelElms = eval(['numel(' currentLevelName ');']);
        levelSize = levelSize * levelElms;
        
        % write loop
        startString = [startString 'for i' num2str(nLevel) '= 1:numel(' currentLevelName '),if ~isempty(' nextLevelName '),'];
        endString = ['end,end,' endString];
        
        % update currentLevelName
        currentLevelName = nextLevelName;
        
        breakString = [breakString,'break,'];
        
    end % for nLevel = 1:numberOfLevels
    
    % build out-range (to which the new data will be assigned)
    outRangeLoopString = ['outRange = [''(''];',...
            'reassignRange = [''(''];',...
            'dataSize = size(' currentLevelName rangeString '); ',...
            'for nDims = 1:max(length(dataSize),dim),',... % we might be catenating along an other dim!
            'if nDims == 1,',...
            'sepStr = '''';',...
            'else,',...
            'sepStr = '','';',...
            'end,',...
            'if nDims == dim,',...
            'ranStr = ''dataCounter+1:newDataCounter'';',...
            'reRanStr = ''1:dataCounter'';',...
            'else,',...
            'ranStr = '':'';',...
            'reRanStr = '':'';',...
            'end,',...
            'outRange = [outRange sepStr ranStr];',...
            'reassignRange = [reassignRange sepStr reRanStr];',...
            'end,',...
            'outRange = [outRange '')''];',...
            'reassignRange = [reassignRange '')''];'];
    eval([startString outRangeLoopString breakString endString]);    
    
    % introduce separator (we need to know currentLevelName)
    if isseparator
        % the first time the loop runs, we want to calculate the
        % perpendicular dimensions of the data that we catenate, so we can
        % fill the separator array. From then on, we just put it into the
        % new array above the new entry
        separatorString = ['if isFirst4Separator, ',...
                'catDim = size(' currentLevelName rangeString '); ',...
                'catDim(' dimString ') = 1; ',...
                'separatorArray = repmat(' separator ',catDim);',...
                'isFirst4Separator = 0;',...
                'else, ',...
                'newDataCounter = dataCounter + 1;',...
                'out' outRange '=separatorArray; ',...
                'dataCounter = newDataCounter;',...
                'end,'];
    else
        separatorString = '';
    end
    
    % build the middle string with pre-assignment of out
    
    
    
    middleString = ['if isFirst4Init,',...
        'dataSizeInput = size(' currentLevelName rangeString '); ',... % get data size of input
        'dataSize = ones(1,max(length(dataSizeInput),dim));',... % add more dimensions, if necessary
        'dataSize(1:length(dataSizeInput)) = dataSizeInput;',... % fill in input dataSize
        'dataSize(' dimString ') = (dataSize(' dimString ')+isseparator)*levelSize;',... % estimate how big the array has to be
        'dataSizeCrit = dataSize(' dimString ');',... % current critical data size
        'dataSizeIncr = zeros(size(dataSize)); dataSizeIncr(' dimString ') = dataSizeCrit;',...
        'if string,',... % check how to assign
        'out = repmat('' '',dataSize);',... % assign spaces
        'else,',...
        'out = zeros(dataSize);',...
        'end,',... % if string
        'dataCounter = 0;',...
        'isFirst4Init = 0;',... % we do not want to come back
        'end,',... % now check whether we have to reassign out and fill in values
        'dataLength = size(' currentLevelName rangeString ',' dimString '); ',...
        'newDataCounter = dataCounter + dataLength;',...
        'if (newDataCounter + isseparator)>dataSizeCrit,',...
        'outmp = out;',... % store out in tmp array
        'dataSize = dataSize + dataSizeIncr;',...
        'dataSizeCrit = dataSize(' dimString ');',... % current critical data size
        'if string,',... % check how to assign
        'out = repmat('' '',dataSize);',... % assign spaces
        'out' reassignRange ' = outmp' reassignRange '; clear outmp;',...
        'else,',...
        'out = zeros(dataSize);',...
        'out' reassignRange ' = outmp' reassignRange '; clear outmp;',...
        'end,',... % if string
        'end,',... %end if first4init
        'out' outRange ' = ' currentLevelName rangeString ';',... % assign out
        'dataCounter = newDataCounter;',... % update counter and go on
        ' ']; 

    cleanupString = ['out = out' reassignRange,';'];


    % finish loop
    evalString = ['out = []; isFirst4Separator = 1; isFirst4Init = 1;' startString separatorString middleString endString cleanupString];
end

% catenate structure

try
    eval(evalString);
catch
    rethrow(lasterror)
end

