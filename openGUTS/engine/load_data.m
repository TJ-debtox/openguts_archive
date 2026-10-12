function [data_all,error_flag] = load_data(filename,filepath,flag_stage,corr_factor)

% Usage: [data_all,error_flag] = load_data(filename,filepath,flag_stage,corr_factor)
%
% This function loads one or more data files in the standard text format.
% It uses a sub-function <text_parser> to read the text file, line by line,
% and translate it to a Matlab-friendly format. The parser function returns
% a cell array in <data> (a cell for each data set, with structures
% underneath) that is subsequently sent to <prepare_data> that produces the
% internal format (a single structured object <data_all>). The <flag_stage>
% tells it which stage of the analysis we are in. This is needed for
% predictions since the text files are different: they lack survival data
% and also lack the first row with the identifier. Flag 4 is reserved for
% batch processing of the LPx: no need for output to screen.
% 
% Inputs
% <filename>    one or more filenames returned by the Matlab file-open GUI (either cell or string!)
% <filepath>    path to the directory in which the input text files are
% <flag_stage>  flag for the stage of the analysis we are in (1-4, last for batch LPx)
% <corr_factor> can be used to correct exposure profiles for predictions to
%               the correct unit (only used when <flag_stage> = 3 or 4).
% 
% Outputs
% <data_all>    a structured object (cell array) with the following fields:
%
% Once for the object:
% <data_all.hb_hat>   = background hazard rate for all controls together
% <data_all.nr_sets>  = number of data sets
% <data_all.nr_treat> = vector with number of treatments in each data set
% 
% Per data set:
% <data_all.time>  = entire time vector for survival data in this data set
% <data_all.c_twa> = vector with time-weighted averages for all treatments
% <data_all.c_max> = vector with maximum concentration over all treatments
% 
% Per treatment:
% <data_all.surv> = matrix with time, survivors, deaths, prob. surv., Wilson score (min and max)
% <data_all.scen> = events-time matrix with time, concentration and slope
% <data_all.name> = name of this treatment (string)
% 
% <error_flag>  flag to signal that we have errors
% 
% Author     : Tjalling Jager 
% Date       : April 2023
% Web support: <http://www.openguts.info>

% =========================================================================
% Copyright (c) 2018-2024, Tjalling Jager (tjalling@debtox.nl). 
% This file is part of the Matlab version of openGUTS.
% 
% This program is free software: you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation, either version 3 of the License, or (at your
% option) any later version.
%  
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
% Public License for more details.
% 
% You should have received a copy of the GNU General Public License along
% with this program. If not, see <https://www.gnu.org/licenses/>.
% =========================================================================

% TO DO: try to remove the globals. 

global GLO % make structure <GLO> global    

%% BLOCK 1. Load text file, and put in easier format.

if ~iscell(filename) % the GUI makes it a cell array only when multiple files are selected ...
    filename = {filename}; % but we also want a cell array if it's only one
end

data_all   = []; % initialise output as empty as errors prompt a premature return
data       = cell(length(filename),1); % initialise empty cell array
error_str1 = {}; % define error string as empty cell array
error_flag = 0; % unless something happens, error flag is zero

if flag_stage == 1 || flag_stage == 2 % loading survival data for calibration/validation
    % I assume there is a specifically formatted text file available with
    % the input data set(s).
    
    for i_d = 1:length(filename) % run through filenames that are entered
        
        % Note: if the Matlab GUI is used for selecting files, there is no
        % need to check if the file exists, I guess.
        if exist([filepath,filename{i_d}],'file') == 2
            % Note: the <text_parser> (sub-function below) is used to parse
            % the text and convert it to a structured cell array <data>.
            [data{i_d},error_str] = text_parser(filename{i_d},filepath);
            error_str1 = cat(1,error_str1,error_str); % add the error to the previous ones
        else % if the file is not available in the right folder, produce an error
            error_str1 = cat(1,error_str1,['There is no data set with filename ',filename{i_d},' in directory input_data']);
        end
        
    end
    
else % loading exposure profiles for predictions
    % I assume that this is a simple text file, without any headers or
    % other info, just two columns (time and concentrations).
    
    for i_d = 1:length(filename) % run through filenames that are entered
        
        if exist([filepath,filename{i_d}],'file') == 2
            data{i_d}.scen = load([filepath,filename{i_d}]); % load the text file and place in cell array
            if size(data{i_d}.scen,2) ~= 2 % check if it has two columns!
                error_str1 = cat(1,error_str1,['File ',filename{i_d},': exposure profiles must contain two columns: 1 time vector and 1 concentration vector']);
            end
            data{i_d}.name = {filename{i_d}(1:end-4)}; % just use filename without extension as name
        else % if the file is not available in the right folder, produce an error
            error_str1 = cat(1,error_str1,['There is no data set with filename ',filename{i_d},' in directory input_profile.']);
        end

    end
    
end

if ~isempty(error_str1)
    errordlg(error_str1,'Errors in input file'); % display a message box with the errors
    error_flag = 1; % signal main script that we ran into errors
    return
end

% If there are no major errors so far, we can go to <prepare_data>.
% Possible correction factor for different units of exposure profile,
% compared to calibration data in <corr_factor>. However, it may be missing
% from the input, and then we need to set it to 1.
if ~exist('corr_factor','var') || corr_factor == 0 % if it does not exist in the input, or is zero ...
    corr_factor = 1; % just set it to 1
end
[data_all,error_str2] = prepare_data(data,flag_stage,corr_factor); % prepare the data into a structured object

if ~isempty(error_str2)
    errordlg(error_str2,'Errors in input file'); % display a message box with the errors
    error_flag = 1; % signal main script that we ran into errors
    return
end

%% BLOCK 2. Print which data sets are loaded for what purpose.

% First, see if there is an output folder for this specific studyname, and
% otherwise create it.
savestr = ['output_report',filesep,GLO.basenm];
if ~exist(savestr,'dir')
    mkdir(savestr); % if not, create it
end

if flag_stage ~= 4 % <flag_stage> 4 is used for LPx batch processing: no need to say which file is loaded
    
    diary ([savestr,filesep,GLO.basenm,'.out']) % collect output in the diary *.out
    switch flag_stage
        case 1
            disp('Following data sets are loaded for calibration:')
            GLO.filename = filename; % make the loaded filenames global so they are available for printing with calibration output in <calc_optim>
        case 2
            disp('Following data sets are loaded for validation:')
        case 3
            disp('Following exposure profiles are loaded for predictions:')
    end
    for i_d = 1:length(filename) % run through data sets loaded
        disp(['    set: ',num2str(i_d),', file: ',filename{i_d}]) % and display their filename
    end
    diary off % close diary
    
end

function [data,error_str] = text_parser(filename,filepath)

% A sub-function to read the formatted text file for the input data. It
% reads one line at a time and interprets is. It returns the structured
% cell array <data>. This is needed for the Matlab version to be able to
% work with formatted text files as proposed for the software. This is thus
% not needed for the software.
%
% Outputs
% <data>            cell array, with one structured block per data set. 
% fields in <data>, for each data set:
%   <data.name>       names of the treatments as cell array of strings
%   <data.surv>       survival data set without headers (matrix, with time in first column)
%   <data.scen>       exposure scenario data without headers (matrix, with time in first column)
%   <data.concunit>   concentration unit of the exposure data set.
% 
% <error_str>  a array with strings of error messages

fileID    = fopen([filepath,filename]); % open file for reading
nxtline   = fgetl(fileID);   % reads one line of text
error_str = {}; % define error string as empty cell array
data      = []; % initialise as empty as errors prompt a premature return

% NOTE: function <contains> was introduced in Matlab R2016b, so earlier
% versions throw an error here.
while ~contains(nxtline,'Survival time','IgnoreCase',true) % continue until we find the line that starts with "Survival time"
    nxtline = fgetl(fileID); % reads one line of text
     if nxtline == -1
         error_str = sprintf('Format error in data set %1s: unexpected end of file',filename);
         return
     end
end
surv_header = strsplit(nxtline,'\t','CollapseDelimiters',0); % split text string at the tab and makes cell array
% Multiple tabs are treated as multiple tabs (hopefully, this catches empty
% cells as well).

if length(surv_header) ~= length(unique(surv_header))
    error_str = sprintf('Format error in data set %1s: treatment names are not unique',filename);
    return
end
if sum(strcmpi(surv_header,'Control')) > 1 % is this already caught by prvious check?
    % Note: <strcmpi> compares strings without caring for case.
    error_str = sprintf('Format error in data set %1s: there can be only one treatment identified as control',filename);
    return
end

nxtline    = fgetl(fileID); % reads the next line from the file
surv_nmbrs = []; % initialise an empty matrix for survivors

% NOTE: function <contains> was introduced in Matlab R2016b, so earlier
% versions throw an error here.
while ~contains(nxtline,'Concentration','IgnoreCase',true) % continue until we find the line that starts with "Concentration"
    
    cell_line = strsplit(nxtline,'\t','CollapseDelimiters',0); % return cell array with expressions, split at the tabs
    if length(surv_header) ~= length(cell_line)
        error_str = sprintf('Format error in data set %1s: survival data does not have same number of columns at each time point',filename);
        return
    end
    
    surv_nmbrs = [surv_nmbrs; nan(1,length(cell_line))]; %#ok<AGROW> initialise a line of NaNs
    for i = 1:length(cell_line) % run through all cells of the array
        surv_i = str2double(cell_line{i}); % convert string to number
        surv_nmbrs(end,i) = surv_i; % note that any text fields will automatically become NaN
    end
    nxtline = fgetl(fileID); % reads next line from file
    if nxtline == -1
        error_str = sprintf('Format error in data set %1s: unexpected end of file',filename);
        return
    end
end 

conc_unit = strsplit(nxtline,'\t'); % split text string at the tab and makes cell array
% The initial version of the data format applied a space between the words
% "Concentration unit:" and the unit itself. We modified it for the
% software to have a tab in between, which is more logical. However, I here
% want to be able to use the old version as well. The C++ software can use
% both formats as well.
if numel(conc_unit) == 1 || isempty(conc_unit{2}) % then the unit is part of the first string
    conc_unit = conc_unit{1}(21:end); % removes the text "Concentration unit: " to leave the unit itself
else % there is a tab between the text and the unit    
    conc_unit = conc_unit{2}; % second element is the unit
end

nxtline = fgetl(fileID); % reads next line from the file

% NOTE: function <contains> was introduced in Matlab R2016b, so earlier
% versions throw an error here.
if ~contains(nxtline,'Concentration','IgnoreCase',true)
    error_str = sprintf('Format error in data set %1s: expected exposure scenario header below concentration unit',filename);
    return
end 
scen_header = strsplit(nxtline,'\t','CollapseDelimiters',0); % split text string at the tab and makes cell array
if length(scen_header) ~= length(surv_header)
    error_str = sprintf('Format error in data set %1s: exposure scenario must have same number of columns as survival data',filename);
    return
end

nxtline    = fgetl(fileID); % reads next line from file
scen_concs = []; % initialise an empty matrix for scenario information

while nxtline ~= -1 % continue until we find the end of file
    
    cell_line = strsplit(nxtline,'\t','CollapseDelimiters',0); % return cell array with expressions
    if length(scen_header) ~= length(cell_line)
        error_str = sprintf('Format error in data set %1s: exposure scenario does not have same number of columns at each time point',filename);
        return
    end
    
    scen_concs = [scen_concs; nan(1,length(cell_line))]; %#ok<AGROW> initialise a line of NaNs
    for i = 1:length(cell_line) % run through all cells of the array
        scen_i = str2double(cell_line{i}); % convert string to number
        scen_concs(end,i) = scen_i; % note that any text fields will automatically become NaN
    end
    nxtline = fgetl(fileID); % reads next line from file
end 

fclose(fileID); % close the file

surv_header = surv_header(2:end); % remove the text in the first column
scen_header = scen_header(2:end); % remove the text in the first column

[tst,loc_scen] = ismember(surv_header,scen_header); % are all identifiers of survival also in the scenarios?
if any(tst==0)
    error_str = sprintf('Format error in data set %1s: some scenarios for survival do not have matching scenario identifier',filename);
    return
end

% If that check goes through, all survival identifiers have a matching
% scenario identifier. However, they may still be in the wrong order. Here,
% I decided to change the order of the scenario data set to match the
% survival data. However, it may be better to produce an error.

% scen_header = scen_header(loc_scen); % put scenario identifiers in same order as survival data
scen_concs(:,2:end) = scen_concs(:,loc_scen+1); % put scenario info in same order (first column is times)

% Collect all information into a structure <data> for output.
data.name     = surv_header;
data.surv     = surv_nmbrs;
data.scen     = scen_concs;
data.concunit = conc_unit;
