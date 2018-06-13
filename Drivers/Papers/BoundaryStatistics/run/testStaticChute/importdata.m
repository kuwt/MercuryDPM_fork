function [out, delimiter, headerlines] = importdata(varargin)
%IMPORTDATA Load data from a file into MATLAB.
%
% IMPORTDATA(FILENAME) loads data from FILENAME into the workspace.
%
% A = IMPORTDATA(FILENAME) loads data from FILENAME into A.
%
% IMPORTDATA(FILENAME, DELIM) loads data from FILENAME using DELIM as 
% the column separator. DELIM must be a string. Use '\t' for tab.
%
% IMPORTDATA(FILENAME, DELIM, HLINE) where HLINE is a number that
% indicates on which line of the file the header text is located, 
% loads data from line HLINE+1 to the end of the file.
%
% [A D] = IMPORTDATA(...) returns the output structure in A, and 
% the delimiter character in D.
%
% [A D H] = IMPORTDATA(...) returns the output structure in A, the 
% delimiter character in D, and the line number of the header in H.
%
% [...] = IMPORTDATA('-pastespecial', ...) loads data from your 
% computer's paste buffer rather than from a file.
%
% IMPORTDATA looks at the file extension to determine which helper function
% to use.  If the extension is mentioned in FILEFORMATS, IMPORTDATA will
% call the appropriate helper function with the maximum number of output
% arguments.  If the extension is not mentioned in FILEFORMATS, IMPORTDATA
% calls FINFO to determine which helper function to use.  If no helper
% function is defined for this extension, IMPORTDATA treats the file as
% delimited text.  Empty outputs from the helper function are removed from
% the result.
%
% NOTE: When the file to be imported is an ASCII text file and IMPORTDATA has
%     trouble importing the file, try TEXTSCAN with a more elaborate argument
%     set than the rather simple application of TEXTSCAN in IMPORTDATA.
% NOTE: When reading of a rather old Excel format fails, try updating the
%     the Excel file to Excel 2000 or 95 format by opening and saving
%     with Excel 2000 or Excel 95.
%
% Examples:
%
%    s = importdata('ding.wav')
% s =
%
%    data: [11554x1 double]
%      fs: 22050 
%
%   s = importdata('flowers.tif');
%   size(s)
% ans =
%
%    362   500     3
%
% See also LOAD, FILEFORMATS, TEXTSCAN, OPEN, LOAD, UIIMPORT.

% Copyright 1984-2006 The MathWorks, Inc.
% $Revision: 1.17.4.21 $  $Date: 2007/05/23 18:54:34 $

error(nargchk(1,3,nargin,'struct'));
error(nargoutchk(0,3,nargout,'struct'));

FileName = varargin{1};

if nargin > 1
    delim.requested = varargin{2};
else
    delim.requested = NaN;
end

if nargin > 2
    requestedHeaderLines = varargin{3};
else
    requestedHeaderLines = NaN;
end

out = [];

if nargout > 1
    delim.printed = [];
    delimiter = [];
end

if nargout > 2
    headerlines = 0;
end

fromClipboard = false;

if strcmpi(FileName,'-pastespecial')
    % fetch data from clipboard
    cb = clipboard('paste');
    if isnan(delim.requested)
        delim.printed = guessdelim(cb);
        delim.requested = delim.printed;
    else
        delim.printed = sprintf(delim.requested);
    end
    %When importing data from clipboard, do not warn of format mismatches.
    bWarn = false;
    [out.data, out.textdata, headerlines] = parse(cb, ...
        delim, requestedHeaderLines, bWarn);
    out = LocalRowColShuffle(out);
    delimiter = delim.printed;
    fromClipboard = true;
else

    % attempt extracting descriptive information about file
    warnState = warning('off', 'MATLAB:xlsfinfo:ActiveX');
    try
        [FileType,openCmd,loadCmd,descr] = finfo(FileName);
        warning(warnState);
    catch
        warning(warnState);
        error('MATLAB:importdata:FileNotFound', 'Unable to open file.');
    end
    
    % Test success of FINFO call
    if strcmp(descr,'FileInterpretError')

        % Generate a warning that FINFO could not interpret data file
        Message = 'File contains uninterpretable data.';
        warning('MATLAB:IMPORTDATA:InvalidDataSection',Message)
        out.data=[]; % return an empty matrix object
        out.textdata={}; % return an empty cell object
        return
    end

    delim.printed = NaN;
    delimiter = delim.printed;
    %Just in case we found incorrect command, i.e. the name was a
    %coincidence, we'll try to load, but if we fail, we'll try to load
    %using other means.
    loaded = 0;
    if ~isempty(loadCmd) && ~strcmp(loadCmd,'importdata')
        try
            out.data = feval(loadCmd, FileName);
            loaded = 1;
        catch
            err = lasterror;
            err.message = [err.message sprintf('\t Using %s command', loadCmd)];
        end
    end
    if (loaded == 0)
        switch FileType
            case 'xls'
                warnState = warning('off', 'MATLAB:xlsread:Mode');
                warnState(2) = warning('off', 'MATLAB:xlsread:ActiveX');
                % special case for single sheet
                if length(descr) == 1
                    [n,s] = xlsread(FileName,descr{1});
                    if ~isempty(n)
                        out.data = n;
                    end
                    if ~isempty(s)
                        out.textdata = s;
                    end
                    out = LocalRowColShuffle(out);
                else
                    % top level fields so assignments below work right
                    for i = 1:length(descr)
                        [n,s] = xlsread(FileName,descr{i});
                        if ~isempty(n)
                            out = AssignStructureUsingGenvarname('data', out, n, descr{i});
                        end
                        if ~isempty(s)
                            out = AssignStructureUsingGenvarname('textdata', out, s, descr{i});
                        end
                        
                        if ~isempty(s) && ~isempty(n)
                            [dm, dn] = size(n);
                            [tm, tn] = size(s);
                            if tn == 1 && tm == dm
                                out = AssignStructureUsingGenvarname('rowheaders', out, s(:,end), descr{i});
                            elseif tn == dn
                                out = AssignStructureUsingGenvarname('colheaders', out, s(end,:), descr{i});
                            end
                        end
                    end
                end
                warning(warnState); %#ok
            case 'wk1'
                [out.data, out.textdata] = wk1read(FileName);
                out = LocalRowColShuffle(out);
            case 'avi'
                out = aviread(FileName);
            case 'im'
                [out.cdata, out.colormap, out.alpha] = imread(FileName);
            case {'au','snd'}
                [out.data, out.fs] = auread(FileName);
            case 'wav'
                [out.data, out.fs] = wavread(FileName);
            case 'mat'
                wasError = false;
                try
                    if ~isempty(whos('-file',FileName))
                        % call load with -mat option
                        out = load('-mat',FileName);
                    else
                        wasError = true;
                    end
                catch
                    wasError = true;
                end
                if wasError
                    % call load with -ascii option
                    out = load('-ascii',FileName);
                end
                if isempty(out)
                    error('MATLAB:importdata:InvalidFile', 'Unable to load file.  Not a MAT file.');
                end
            otherwise
                % try to treat as hidden mat file
                try
                    out = load('-mat',FileName);
                catch
                    out = [];
                end
                if isempty(out)
                    % file is an unknown format, treat it as text
                    [out, delimiter, headerlines] = LocalTextRead(FileName, delim, requestedHeaderLines);
                end
        end
    end
end

if (length(out) == 1 && isstruct(out))
    % remove empty fields from output struct
    names = fieldnames(out);
    for i = 1:length(names)
        if isempty(out.(names{i}))
            out = rmfield(out, names{i});
        end
    end
    if ~isempty(out)
        % flatten output struct if single variable
        names = fieldnames(out);
        if length(names) == 1
            out = out.(names{1});
        end
    end
end

if isempty(out) || (isstruct(out) && isempty(fieldnames(out)))
    if fromClipboard
        error('MATLAB:importdata:InvalidClipboardData', 'Clipboard does not contain any recognizable data.');
    else
        err = lasterror;
        default_message = 'Unable to load file.\nUse TEXTSCAN or FREAD for more complex formats.';
        if ~isempty(err)
            err.message = [err.message sprintf(['\n' default_message])];
            rethrow(err);
        else
            error('MATLAB:importdata:UnableToLoad', default_message);
        end
    end
end


function [out, delimiter, headerlines] = LocalTextRead(filename, delim, hlines)
% get the delimiter for the file
if isnan(delim.requested)
    fid = fopen(filename);
    str = fread(fid, 4096,'*char')';
    fclose(fid);
    delim.printed = guessdelim(str);
    delim.requested = delim.printed;
    
else
    delim.printed = sprintf(delim.requested);
end
delimiter = delim.printed;
% try load first (it works with tabs, spaces, and commas)
out = [];
if isnan(hlines) && ~isempty(findstr(delim.printed, sprintf('\t ,')))
    try
        out = load('-ascii', filename);
        headerlines = 0;
    catch
        out = '';
    end
end
if isempty(out)
    try
        fileString = fileread(filename);
        bWarn = true;
        [out.data, out.textdata, headerlines] = parse(fileString, delim, hlines, bWarn);
        out = LocalRowColShuffle(out);
    catch
        out = [];
        headerlines = 0;
    end
end

function out = LocalRowColShuffle(in)

out = in;

if isempty(in) || ~isfield(in, 'data') || ~isfield(in,'textdata') || isempty(in.data) || isempty(in.textdata)
    return;
end

[dm, dn] = size(in.data);
[tm, tn] = size(in.textdata);

if tn == 1 && tm == dm
    % use as row headers
    out.rowheaders = in.textdata(:,end);
elseif tn == dn
    % use last row as col headers
    out.colheaders = in.textdata(end,:);
end


function [numericData, textData, numHeaderRows] = parse(fileString, ...
    delim, headerLines, bWarnOnMismatch)
%This is the function which takes in a string and parses it into
%"spreadsheet" type data.
error(nargchk(1,4,nargin, 'struct'));

numericData = [];
textData = {};
numHeaderRows = 0;

% gracefully handle empty
if isempty(fileString) && isempty(regexp(fileString,'\S','once')); 
    %regexp is faster than all(isspace(fileString));
        return;
end

% validate delimiter
if length(delim.printed) > 1
    error('MATLAB:importdata:InvalidDelimeter', 'Multi character delimiters not supported.')
end

if nargin < 3
    headerLines = NaN;
end

%Arbitrarily set the maximum size for a line, used to calculate this but it
%was slow, so we went with all in one line or the old maximum it would
%check.
bufsize = min(1000000, max(numel(fileString),100)) + 5;

% use what user asked for header lines if specified
[numDataCols, numHeaderRows, numHeaderCols, numHeaderChars] = ...
    analyze(fileString, delim, headerLines, bufsize);

% fetch header lines and look for a line of column headers
headerLine = {};
headerData = {};
origHeaderData = headerData;
useAsCells = 1;
    
if numHeaderRows
    firstLineOffset = numHeaderRows - 1;  
    if firstLineOffset > 0    
        try
            headerData = textscan(fileString,'%[^\n\r]',firstLineOffset,...
                'whitespace','','delimiter','\n','bufsize', bufsize);
        catch
            %Just in case.  Should not happen.
            headerData = {''};
        end
        origHeaderData = headerData{1};
        if numDataCols
            headerData = [origHeaderData, cell(length(origHeaderData), numHeaderCols + numDataCols - 1)];
        else
            headerData = [origHeaderData, cell(length(origHeaderData), numHeaderCols)];
        end
    else
       headerData = cell(0, numHeaderCols + numDataCols);
    end
    try
        Data = textscan(fileString,'%[^\n\r]',1,'headerlines',firstLineOffset,...
            'delimiter',delim.requested,'bufsize', bufsize);
    catch
            %Just in case.  Should not happen.
            Data = {''};
    end
    headerLine = Data{1};
    origHeaderLine = headerLine;
    
    useAsCells = 0;
    
    if ~isempty(delim.printed) && ~isempty(headerLine) && ~isempty(strfind(deblank(headerLine{:}), delim.printed))
        cellLine = split(headerLine{:}, delim);
        if (isequal(headerLine{:}(end),delim.printed))
            cellLine(end+1) = {[]};
        end
        if length(cellLine) == numHeaderCols + numDataCols
            headerLine = cellLine;
            useAsCells = 1;
        end
    end
    
    if ~useAsCells
        if numDataCols
            headerLine = [origHeaderLine, cell(1, numHeaderCols + numDataCols - 1)];
        else
            headerLine = [origHeaderLine, cell(1, numHeaderCols)];
        end
    end
end

if isempty(delim.printed)
    formatString = [repmat('%s', 1, numHeaderCols) repmat('%n', 1, numDataCols)];
else
    formatString = [repmat(['%[^' delim.printed ']'], 1, numHeaderCols) repmat('%n', 1, numDataCols)];
end

% now try for the whole shootin' match
try
    if numDataCols
		 %When the delimiter is a space, multiple spaces do NOT mean
		 %multiple delimiters.  Thus, call textsscan such that it will
		 %treat them as one.

		if isequal(delim.printed,' ')
	        Data = textscan(fileString,formatString,'headerlines',numHeaderRows,...
                'CollectOutput', true, 'bufsize', bufsize);
		else
	        Data = textscan(fileString,formatString,'delimiter',delim.requested,...
                'bufsize', bufsize, 'headerlines',numHeaderRows, 'CollectOutput', true);
		end
		if (numHeaderCols)
		    numericData = Data{2};
		else
	        numericData = Data{1};
		end
    end
    wasError = false;
catch
    wasError = true;
end


if nargout > 1
    if numHeaderCols > 0
    	textData = cell(size(Data{1}, 1), numDataCols+numHeaderCols);
        textData(:, 1:numHeaderCols) =  Data{1};
    end

    if ~isempty(headerLine)
        textData = [headerLine; textData];
    end

    if ~isempty(headerData)
        textData = [headerData; textData];
    end    
end
clear('Data');

if (numDataCols && numHeaderCols && (size(textData, 1) ~= numHeaderRows + size(numericData, 1)))
    wasError = true;
end
    
% if the first pass failed to read the whole shootin' match, try again using the character offset
if wasError && numHeaderChars
    wasError = false;
    
    % rebuild format string
    formatString = ['%' num2str(numHeaderChars) 'c' repmat('%n', 1, numDataCols)];

    %When the delimiter is a space, multiple spaces do NOT mean
    %multiple delimiters.  Thus, call textscan such that it will
    %treat them as one.
	if isequal(delim.printed,' ')
        Data = textscan(fileString,formatString,'headerlines',numHeaderRows,...
            'returnonerror',1, 'CollectOutput', true);
	else
        Data = textscan(fileString,formatString,'delimiter',delim.requested,...
            'headerlines',numHeaderRows,'returnonerror',1,'CollectOutput', true);
	end
	
    textCharData = Data{1};
    numericData = Data{2};
    numHeaderCols = 1;
    if ~isempty(numericData)
        numRows = size(numericData, 1);
    else
        numRows = length(textCharData);
    end
    
    if numDataCols
        headerData = [origHeaderData, cell(length(origHeaderData), numHeaderCols + numDataCols - 1)];
    else
        headerData = [origHeaderData, cell(length(origHeaderData), numHeaderCols)];
    end
    
    if ~useAsCells
        if numDataCols
            headerLine = [origHeaderLine, cell(1, numHeaderCols + numDataCols - 1)];
        else
            headerLine = [origHeaderLine, cell(1, numHeaderCols)];
        end
    end
   
    if nargout > 1 && ~isempty(textCharData)
        textCellData = cellstr(textCharData);
        if ~isempty(headerLine)
            textData = [headerLine;
                textCellData(1:numRows), cell(numRows, numHeaderCols + numDataCols - 1)];
        else
            textData = [textCellData(1:numRows), cell(numRows, numHeaderCols + numDataCols - 1)];
        end

        if ~isempty(headerData)
            textData = [headerData; textData];
        end
    end
end

if bWarnOnMismatch && wasError
	warning('MATLAB:importdata:FormatMismatch', ...
		'An unexpected format mismatch was detected.\nPlease check results against original file.');
end

if nargout > 1 && ~isempty(textData)
	textData = TrimTrailing(@(x)cellfun('isempty', x), textData);
end

if ~isempty(numericData) 
  %THOMAS' CHANGE:
    %numericData = TrimTrailing(@(x)(isnan(x)), numericData);
end


function [numColumns, numHeaderRows, numHeaderCols, numHeaderChars] = ...
    analyze(fileString, delim, header, bufsize)
%ANALYZE count columns, header rows and header columns

numColumns = 0;
numHeaderRows = 0;
numHeaderCols = 0;
numHeaderChars = 0;

if ~isnan(header)
    numHeaderRows = header;
end

try
    Data = textscan(fileString,'%[^\n\r]',1,'headerlines',numHeaderRows,...
       'delimiter',delim.requested,'bufsize', bufsize);
catch
    %Just in case.  Should not happen.
    return;
end
thisLine = Data{1};

if isempty(thisLine)
    return;
end
thisLine = thisLine{:};
[isvalid, numHeaderCols, numHeaderChars] = isvaliddata(thisLine, delim);

if ~isvalid
    numHeaderRows = numHeaderRows + 1;
    try
        Data = textscan(fileString,'%[^\n\r]',1,'headerlines',numHeaderRows,...
            'delimiter',delim.requested, 'bufsize', bufsize);
    catch
        %Just in case.  Should not happen.
        return;
    end
    thisLine = Data{1};
    if isempty(thisLine)
        return;
    end
    thisLine = thisLine{1};
    
    [isvalid, numHeaderCols, numHeaderChars] = isvaliddata(thisLine, delim);
    if ~isvalid
        newLines = strfind(fileString,sprintf('\n'));
        
        if isempty(newLines)% must be a MAC
            newLines = strfind(fileString,sprintf('\r'));
            if isempty(newLines) % do not know what to do with this...
                return
            end
        end
        % deal with the case where the file is not terminated with a \n
        if newLines(end) ~= length(fileString)
            newLines(end + 1) = length(fileString);
        end
        numLines = length(newLines);
        newLines(end+1) = length(fileString) + 1;
        while ~isvalid 
            if(numHeaderRows == numLines)
                break;
            end
            
            if (numHeaderRows >= 1000)
                %Assume no data.
                numHeaderRows = numLines + 1;
                break;
            end            
            % stop now if the user specified a number of header lines
            
            if ~isnan(header) && numHeaderRows == header
                break;
            end
            numHeaderRows = numHeaderRows + 1;
            
            thisLine = fileString((newLines(numHeaderRows)+1):(newLines(numHeaderRows+1)-1));
            if isempty(thisLine)
                break;
            end
            [isvalid, numHeaderCols, numHeaderChars] = isvaliddata(thisLine, delim);
        end     
    end
end

% This check could happen earlier
if ~isnan(header) && numHeaderRows >= header
    numHeaderRows = header;
end

if isvalid
    % determine num columns
    %remove trailing spaces.  Spaces are different from other delimiters.
    thisLine = regexprep(thisLine, ' +$', '');
    delimiterIndexes = strfind(thisLine, delim.printed);
    if all(delim.printed ==' ') && length(delimiterIndexes) > 1
        delimiterIndexes = delimiterIndexes([true diff(delimiterIndexes) ~= 1]);
		delimiterIndexes = delimiterIndexes(delimiterIndexes > 1);
    end
    
    % format string should have 1 more specifier than there are delimiters
    numColumns = length(delimiterIndexes) + 1;
    if numHeaderCols > 0
        % add one to numColumns because the two set of columns share a delimiter
        numColumns = numColumns - numHeaderCols;
    end
end


function [status, numHeaderCols, numHeaderChars] = isvaliddata(fileString, delim)
% ISVALIDDATA delimiters and all numbers or e or + or . or -
% what about single columns???

numHeaderCols  = 0;
numHeaderChars = 0;

if isempty(delim.printed)
    % with no delimiter, the line must be all numbers, +, . or -
    status = isdata(fileString);
    return
end

status = 0;
delims = strfind(fileString, delim.printed);
if isempty(delims)
    % a delimiter must occur on each line of data (or there is no delimiter...)
    return
end

% if there is data at the end of the line, it's legit

checkstring = fileString(delims(end)+1:end);
if isdata(checkstring)
    try
        [cellstring indices] = split(fileString, delim); 
        numNonEmptyCols = find(cellfun('isempty',deblank(cellstring)) == false, 1, 'last');
        numHeaderCols = maxNotData(cellstring); 
		% use contents of 1st data data cell to find num leading chars
        if numHeaderCols > 0
            numHeaderChars = indices(numHeaderCols);
        end
        if (numHeaderCols == numNonEmptyCols)
            numHeaderCols = 0;
            numHeaderChars = 0;
        else
    		status = 1;
        end
    catch
        numHeaderCols = 0;
        numHeaderChars = 0;
    end
end

function index = maxNotData(cellstring)
len = length(cellstring);
index = 0;
for i = len:-1:1
    if ~isdata(cellstring{i})
        index = i;
        break;
    end
end
        
function status = isdata(fileString)
%ISDATA true if string can be shoved into a number or if it's allwhite
if isempty(regexp(fileString, '\S', 'once'))
    status = 1;
else
    [a,b,c] = sscanf(fileString, '%g');
    status = isempty(c);
end

function [cellOut indOut] = split(fileString, delim)
%SPLIT rip string apart
if delim.printed == ' '
    %Multiple spaces are often used as a "fixed-width", thus we treat them
    %differently.
    cellOut = textscan(fileString,'%s','delimiter',delim.requested,...
            'multipleDelimsAsOne', 1, 'whitespace','');
    indOut = regexp(strtrim(fileString),' [^ ]') ; 
else
    cellOut = textscan(fileString,'%s','delimiter',delim.requested,...
        'whitespace','');    
	indOut = regexp(fileString,delim.printed) ;
end
cellOut = (cellOut{1})';

function out = AssignStructureUsingGenvarname(name_of_field, out, data, descr)
% Assign into a structure
if ~isfield(out, name_of_field)
    out.(name_of_field) = struct();
end
name = genvarname(descr,fieldnames(out.(name_of_field)));
out.(name_of_field).(name) = data;


function out = TrimTrailing(operation, out)
% Trim trailing that use a certain operation
cols = size(out,2);
while cols >= 1
    if ~all(operation(out(:,cols)))
        break;
    end
    cols = cols - 1;
end
% trim trailing empty rows from textData
rows = size(out,1);
while rows >= 1
    if ~all(operation(out(rows,1:cols)))
        break;
    end
    rows = rows - 1;
end
if rows < size(out,1) || cols < size(out,2)
    out = out(1:rows,1:cols);
end

