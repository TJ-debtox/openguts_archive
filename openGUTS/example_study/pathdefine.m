function pathdefine

%  Syntax: pathdefine
%
% A tiny piece of 'intelligent' code to add the <engine> directory to the
% Matlab Path. Should work on UNIX, Mac and Windows, as long as the script
% from which this function is called is in a subdirectory of the <openGUTS>
% directory. The addition to the path is temporary, and forgotten when
% Matlab is restarted.
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

% first, make sure the folders for the input and output files exists
if ~exist('output_sample','dir')
    mkdir output_sample; % if not, create it
end
if ~exist('output_report','dir')
    mkdir output_report; % if not, create it
end
if ~exist('input_data','dir')
    mkdir input_data; % if not, create it
end
if ~exist('input_profile','dir')
    mkdir input_profile; % if not, create it
end
    
%% Search and add engine folder to the path
% This section checks whether openGUTS/engine is in the path. If not,
% it adds it. This works both on PC and Unix (and probably also on Linux),
% as it makes use of the Matlab built-in <filesep> for the file separator
% on the OS that is used.

% NOTE: function <contains> was introduced in Matlab R2016b, so earlier
% versions throw an error here.
if ~contains(path,['openGUTS',filesep,'engine']) % <openGUTS/engine> is not already in the path ... add it!
   
    curdir = pwd; % ask for path to current directory as string
    k = strfind(curdir,'openGUTS'); % find where <openGUTS> is in the path
    
    if isempty(k) % if it is not in there, produce an error
        error('Make sure you place your scripts somewhere within the folder openGUTS!')
    end
    
    % Calculate the path for <openGUTS>, add the <\engine>, and add to path
    addpath([curdir(1:k+7),filesep,'engine'],'-begin'); % add openGUTS\engine to the path
    
end

%% Remove any BYOM elements from the path
% This will generally be superfluous, but there may be situations where you
% run openGUTS-Matlab and BYOM after each other. This section removes all
% elements with BYOM\engine from the path. This also removes any
% BYOM/engine_par folders from the path.

% NOTE: function <contains> was introduced in Matlab R2016b, so earlier
% versions throw an error here.
if contains(path,['BYOM',filesep,'engine']) % then user has previously ran BYOM!
    path_cells = split(path,';'); % turn path into cell array
    for i = 1:length(path_cells)
        if contains(path_cells{i},['BYOM',filesep,'engine'])
            rmpath(path_cells{i}); % remove BYOM engine from the path
        end
    end
end

