function [scen_worst,error_flag] = calc_batch_lpx(filename,filepath,sdit,corr_factor,flag)

% Usage: [scen_worst,error_flag] = calc_batch_lpx(filename,filepath,sdit,corr_factor,flag)
%
% Function to calculate LPx values (with or without CIs) in batch form. A
% list of file names should be entered as cell array in <filename>. The
% function expects text files with two columns: time (days) and exposure
% concentration. The <corr_factor> is for an optional element to correct
% possible difference in units for the exposure profile (factor is
% multiplied with the concentration in the exposure profile). We decided
% not to use this factor for the software, but I keep it in for possible
% future updates.
%
% Inputs
% <filename>    one or more filenames returned by the Matlab GUI (either cell or string!)
% <filepath>    path to the directory in which the profiles are
% <sdit>        switch for SD (1) or IT (2)
% <corr_factor> optional correction when different units would be used for
%               profiles than for calibration
% <flag>        1) default, without CIs (as in standalone version), 2) with CIs
% 
% Outputs
% <scen_worst>  filename of the scenario yielding the lowest LP50 (string)
% <error_flag>  flag to signal main script that there is an error in an input file
% 
% Author     : Tjalling Jager
% Date       : November 2019
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

global GLO SETTINGS_PRED % make globals available

%% BLOCK 1. Initial things.

if ~iscell(filename) % the Matlab file-input GUI makes it a cell array only when multiple files are selected ...
    filename = {filename}; % but we also want a cell array if it's only one (this would be a bit silly for batch processing but it could happen)
end

disp(' ')
disp('Starting batch processing ...')
diary (['output_report',filesep,GLO.basenm,filesep,GLO.basenm,'.out']) % collect output in the diary *.out
disp(' ')

if isempty(corr_factor) || corr_factor == 0
    corr_factor = 1; % possible correction factor for different units of exposure profile, compared to calibration data
end

F_LPx     = SETTINGS_PRED.F_LPx; % read the effect levels for LPx from the global (vector)
nr_files  = length(filename);    % number of files to run through
LPx_coll  = zeros(nr_files,length(F_LPx)+1); % initialise matrix to collect LPx's and an identifier in first columns for the data set
LPlo_coll = zeros(nr_files,length(F_LPx)+1); % initialise matrix to collect lower CIs of LPx's and an identifier in first columns for the data set
LPhi_coll = zeros(nr_files,length(F_LPx)+1); % initialise matrix to collect upper CIs of LPx's and an identifier in first columns for the data set

%% BLOCK 2. Batch processing.
% At this moment, I do not create a dedicated multi-panel plot (to show all
% profiles together) in this function. I simply use the standard plotting
% function in <plot_main>. Note that I could have loaded all scenarios and
% placed them into one data object (the function <calc_lpx> can deal with
% that). However, that would put unnecesary strain on memory use (people
% might run hundreds of profiles).

% Use the fancy progress bar from Matlab
if flag == 1
    f = waitbar(0,{sprintf('Bar shows fraction of %1.0f exposure profiles calculated',nr_files),'Cancel will stop the analysis *after* finishing this profile!'},'Name','Batch calculation for LPx','CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
else
    f = waitbar(0,{sprintf('Bar shows fraction of %1.0f exposure profiles calculated',nr_files),'Confidence intervals will be calculated, so be patient ...','Cancel will stop the analysis *after* finishing this profile!'},'Name','Batch calculation for LPx','CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
end

for i_d = 1:nr_files % run through scenarios
    
    [data_pred,error_flag] = load_data(filename(i_d),filepath,4,corr_factor); % load and prepare the data set for stage 4 (batch processing)
    if error_flag == 1 % see if there are errors spotted in the input file
        % Note: <load_data> will display a message box with the errors
        disp(['Batch calculations for LPx stopped because of error(s) in file: ',filename{i_d}])
        delete(f)
        return % simply return (stop and do nothing)
    end

    if flag == 1 % standard, no CIs
        [LPx,LPlo,LPhi] = calc_lpx(data_pred,3,sdit); % calculates LPx values (with flag for batch work) (vector, as there is one data set at a time, but two levels of x)
    else % calculate with CIs
        [LPx,LPlo,LPhi] = calc_lpx(data_pred,4,sdit); % calculates LPx values (with flag for batch work with CIs) (vector, as there is one data set at a time, but two levels of x)
        % NOTE: LPlo and LPhi need to be collected into LPx_coll. We can
        % also do that when no CIs are calculated, as they are defaulted as
        % NaNs by <calc_lpx>.
    end
        
    LPx_coll(i_d,:)  = [i_d LPx];  % collect the output into the matrix, with an identification in first column
    LPlo_coll(i_d,:) = [i_d LPlo]; % collect the output into the matrix, with an identification in first column
    LPhi_coll(i_d,:) = [i_d LPhi]; % collect the output into the matrix, with an identification in first column
    
    plot_main(data_pred,sdit,[3 3],LPx) % plot of data (stage 3, plot model on, invisible plot), pass on the calculated LPx values
    % Note that the second 3 is a flag for <plot_main> to plot model curves
    % (always without CIs), and call <make_fig> to make an invisible figure
    % window (figure is saved as file in [output_report]. This is all a bit
    % awkward, and does not have to be implemented like this for the
    % software.
    
    % Progress bar from Matlab: check if cancel is pressed
    if ~isempty(f)
        if getappdata(f,'canceling')
            break % break the <for> loop
        end
        waitbar(i_d/nr_files,f) % update waitbar
    end
    
end

delete(f) % delete progress bar again
if i_d < nr_files % if cancel was pressed, display that on screen
    fprintf('Batch calculation of LPx stopped by user after %1.0f of %1.0f sets \n',i_d,nr_files)
    LPx_coll  = LPx_coll(1:i_d,:);  % only keep values that were calculated
    LPlo_coll = LPlo_coll(1:i_d,:); % only keep values that were calculated
    LPhi_coll = LPhi_coll(1:i_d,:); % only keep values that were calculated
end

%% BLOCK 3. Make a table with sorted LPx values on screen.
% Note that I decided to sort based on LP50. Users might rather see sorting
% on LP10, which could be slightly different. However, that could easily be
% changed below. Note that the standalone software does *not* sort the
% scenarios.

sort_on = 50; % x in LPx to sort on (by default 50%); must be a value available in <F_LPx> 
ind_LPx = find(F_LPx==sort_on/100); % index to the LPx for sorting in the LPx vector
a = sortrows(LPx_coll,ind_LPx+1); % sort by selected LPx, keeping rows together (and accounting for the extra column for identification)

disp(['LPx values, scenarios ordered by descending risk (sorted on LP',num2str(sort_on),')'])
if flag == 2 && isfield(GLO,'n_lim') && ~isempty(GLO.n_lim)
    disp(['  For CIs, limited set of ',num2str(GLO.n_lim),' random sets from outer hull used'])
    disp('  (this underestimates the width of the CIs)')
end

if sdit == 1
    disp('  Special case: stochastic death (SD)');
else
    disp('  Special case: individual tolerance (IT)');
end
if flag == 1 % standard, no CIs
    fprintf('file analysed             LP%1.0f       LP%1.0f \n',F_LPx(1)*100,F_LPx(2)*100);
else % otherwise, print CIs as well
    fprintf('file analysed             LP%1.0f  (CI)                      LP%1.0f  (CI) \n',F_LPx(1)*100,F_LPx(2)*100);
end
disp('======================================================================================')
for i_d = 1:size(LPx_coll,1) % run through scenarios that were calculated
    % Note: this was <1:length(filename)> in v0.1-0.6, but the cancel
    % button in the <waitbar> required a different approach as not all
    % files may be calculated.
    ie = a(i_d,1); % scenario number (order of increasing risk)
    if flag == 1 % standard, no CIs
        fprintf('%-20s   %#8.4g   %#8.4g \n',filename{ie}(1:end-4),LPx_coll(ie,2),LPx_coll(ie,3));
    else % otherwise, print CIs as well
        fprintf('%-20s   %#8.4g (%#8.4g - %#8.4g)  %#8.4g (%#8.4g - %#8.4g) \n',filename{ie}(1:end-4),LPx_coll(ie,2),LPlo_coll(ie,2),LPhi_coll(ie,2),LPx_coll(ie,3),LPlo_coll(ie,3),LPhi_coll(ie,3));
    end        
end
disp('======================================================================================')
if corr_factor ~= 1
    disp(['A correction factor of ',num2str(corr_factor),' was used on the exposure profiles,'])
    disp('   to correct for differences in units with the calibration data.')
end

scen_worst = filename{a(1)}; % return the name of the worst scenario (lowest LP50) as string
diary off % close the diary

