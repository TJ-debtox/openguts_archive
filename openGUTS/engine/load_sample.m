function [pmat,sample,pmat_print,mll] = load_sample(sdit,flag)

% Usage: [pmat,sample,pmat_print,mll] = load_sample(sdit,flag)
%
% This function loads the MAT file saved by <calc_parspace> and returns the
% saved elements required. The <flag> tells it whether it is called to
% return <pmat> and the <sample> from parameter space for post-calibration
% analysis, or to reconstruct a complete calibration from the saved set.
% 
% Inputs
% <sdit>        switch for SD (1) or IT (2)
% <flag>        1) only return <pmat> and sample, 2) allow reconstructing fit/plot from saved set
% 
% Outputs
% <pmat>        parameter matrix from saved set
% <sample>      sample from parameter space for error propagation (matrix)
% <pmat_print>  matrix with information on each parameter (best value and CI info)
% <mll>         minus log-likelihood; this will be the lowest value found in the run
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

global GLO % make structure GLO global
global DATAcal % make data set object for calibration global 

%% BLOCK 1. Construct the correct filename for loading the MAT file.

if sdit == 1 % this creates a string that equals the name of the saved MAT file, based on SD or IT (should match name construction in <calc_parspace>)
    savestr = [GLO.basenm,'_parspace_SD.mat'];
else
    savestr = [GLO.basenm,'_parspace_IT.mat'];
end

% Load the saved MAT file for this analysis (saved by <calc_parspace>).
if exist(['output_sample',filesep,savestr],'file') == 2 % check if it exists first.
    load(['output_sample',filesep,savestr],'pmat','coll_all','pmat_print','coll_prof_pruned','concunit','fix_hb_cal') 
    % Load all of the saved information from the last <calc_parspace> run.
    % 
    % NOTE: in v1.1, DATAcal was loaded as well, but that variable was only
    % defined as global at the end, for flag=2. Matlab 2023a wants to have
    % all global definitions on top. Now, at the end, DATAcal will be
    % separately loaded for that flag. I don't think it is bad to overwrite
    % DATAcal for flag=1 as well, it is just unnecessary. There may be an
    % additional reason why I restricted overwriting DATAcal for flag=2,
    % but I failed to add a specific comment in the code. Thus, I keep the
    % functionality exactly the same as in v1.1.
else % otherwise, produce an error
    error(['There is no confidence set with filename ',savestr,' in output_sample, so run calibration (with correct settings) first.'])
end

%% BLOCK 2. Reconstruct <coll_all> as we saved only the columns for the fitted parameters.

if isempty(coll_all) % it is empty if we saved a <pmat> with all parameters fixed (silly but possible)
    
    mll    = inf; % just return MLL as infinite
    % Note: the MLL is saved as part of <coll_all> (first row is best fit).
    % However, if no parameters are fitted, <coll_all> is saved as empty,
    % so the MLL is lost. If this is unwanted, we need to make sure that
    % <calc_parspace> saves the MLL in some way.
    
else % otherwise, we have to reconstruct the sample
    
    coll_tmp        = zeros(size(coll_all,1),size(pmat,1)+1); % initialise a new empty matrix with needed size
    coll_tmp(:,end) = coll_all(:,end); % first copy the column with the max. lik. to <coll_tmp>
    
    ind_fit = pmat(:,2) == 1; % logical indices to fitted parameters (vector)
    ind_log = pmat(:,5) == 0; % logical indices to log-scale parameters (vector)
    
    p_hat                = pmat(:,1); % extract the best and fixed values from <pmat> (vector)
    p_hat(ind_log)       = log10(p_hat(ind_log)); % put log-parameters on log10 scale (sample contains log-transformed values for log-scale parameters)
    coll_tmp(:,ind_fit)  = coll_all(:,1:end-1); % copy the *fitted* parameters to the correct column
    coll_tmp(:,~ind_fit) = ones(size(coll_all,1),1) * (p_hat(~ind_fit))'; % copy the *fixed* parameters to the correct column (this copies the fixed values to all rows)
    coll_all             = coll_tmp; % replace <coll_all> by <coll_tmp>
    clear coll_tmp % we can clear the temporary variable <coll_tmp> to save some memory
    
    mll = coll_all(1,end); % extract maximum likelihood value from saved set (<coll_all> is sorted, so it is in the first row)

end

%% BLOCK 3. Create the sample for output, or to make the plot from the saved set.

switch flag
    case 1 % When calling from <calc_conf> or <calc_lcx>/<calc_lpx> return <pmat> and the relevant selection from the sample

        if isempty(coll_all) % it is empty if we saved a <pmat> with all parameters fixed
            sample = []; % all parameters were fixed, so just return an empty sample.
            return % and immediately go back to the function that we're called from
        end

        % For error propagation, we'll use a subset from the sample between
        % two chi-square criteria (defined in <initial_setup>). We work
        % here with the log-likelihood itself, so the chi2 criterion needs
        % to be divided by 2.
        ind_prop1 = find(coll_all(:,end) < mll + 0.5 * GLO.crit_prop(1),1,'last'); % last entry within first criterion
        ind_prop2 = find(coll_all(:,end) < mll + 0.5 * GLO.crit_prop(2),1,'last'); % last entry within second criterion
        
        if isempty(ind_prop1) % <GLO.crit_prop(1)> could be zero (if <GLO.crit_prop> is set to return the entire inner rim), which would lead to empty <ind_prop1>
            ind_prop1 = 1; % start from the first value
        end
        sample = coll_all(ind_prop1:ind_prop2,1:end-1); % take the sets in a defined band around inner rim only (matrix), exclude the likelihood column
        
    case 2 % When calling from <calc_optim> we need to make the plot and return <pmat>, <mll> and <pmat_print>
        
        load(['output_sample',filesep,savestr],'DATAcal') % additionally load DATAcal to make the data set from the saved set the new global one
        % NOTE: this was changed after v1.1 to avoid a warning from Matlab
        % R2023a (defining a global after the variable has been given a value).

        GLO.concunit   = concunit; % make the concentration unit from saved set global
        GLO.fix_hb_cal = fix_hb_cal; % make the setting for background hazard treatment from saved set global
        sample = []; % make sure <sample> is defined so it can be output of the function
        
        if ~isempty(coll_all) % it is empty if we saved a <pmat> with all parameters fixed; then we don't need a plot
             % Reconstruct the plot of parameter space from the saved sample.
             plot_grid(pmat,[],coll_all,coll_prof_pruned,[]); % create plot but don't save plot again ... or should we?
        end       
        
end

