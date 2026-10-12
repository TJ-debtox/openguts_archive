function [data_all,error_str] = prepare_data(data,flag_stage,corr_factor)

% Usage: [data_all,error_str] = prepare_data(data,flag_stage,corr_factor)
% 
% This function translates the data sets entered into the internal format.
% Internally, we like to have a single structured object, with the
% treatments in separate matrices to allow missing data points, and to be
% able to make some preparatory calculations (such as the derivation of
% deaths-per-interval from the survivors). Also, some checks on the data
% set are performed. Note that there can be multiple data sets entered,
% which means that the input <data> is a cell array that can contain
% multiple data sets.
%
% Note: for all data types, multiple data sets are allowed here. However,
% for the standalone software we only allow multiple data sets for
% calibration and batch-processing of LPx. Also, the correction factor
% <corr_factor> is not part of the standalone software.
%
% Inputs
% <data>        structured array with datasets as read from file in <load_sample> with
%               survivors, exposure scenarios, names and unit
% <flag_stage>  flag to indicate which stage of the analysis we are called from
% <corr_factor> can be used to correct exposure profiles for predictions to
%               the correct unit (only used when <flag_stage> = 3 or 4).
%
% Outputs
% <data_all>    A structured object (cell array) with the following fields:
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
% <error_str>   A cell array with strings of error messages (empty if there
%               are no errors)
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

global GLO % make <GLO> global

%% BLOCK 1. Initial things
% Set and check the concentration unit. Check that it is the same when
% multiple data sets are entered, and that it is the same for validation
% and calibration.

nr_sets = length(data); % number of data sets that are entered
error_str = {}; % define error string as empty cell array

switch flag_stage % different steps for each analysis stage
    
    case 1 % then we are called with calibration data sets
        GLO.concunit = data{1}.concunit; % use the unit from the first data set and use it for the global
        for i_d = 2:nr_sets % run through rest of the data sets
            if ~strcmp(GLO.concunit,data{i_d}.concunit) % this compares two strings; if they are unequal produce an error
                error_str = cat(1,error_str,'Data sets entered for calibration have different concentration units');
            end
        end
        
    case 2 % then we are called with validation data sets
        for i_d = 1:nr_sets % run through all of the data sets
            if ~strcmp(GLO.concunit,data{i_d}.concunit) % this compares the strings from the validation sets with the global
                error_str = cat(1,error_str,'Data sets entered for validation have a different concentration unit than for calibration');
            end
        end
        
    otherwise % then we are called for predictions (regular or batch), that means no survival data are used
        % Below, I allow for a correction factor to be used, to correct if
        % the prediction data use a different concentration unit. We decided
        % not to use that for the software, but I leave the option in here.
        % 
        % Next, this section creates some dummy survivor data. The survival
        % data sets are used to shape the exposure profiles later in this
        % function. The dummy data avoids having to modify the core of this
        % function, and allows the exposure data to run through the code
        % smoothly.
        
        % Extract the exposure-scenario matrix from the input <data>, and
        % create dummy survival data. Note that <load_data> already checks
        % whether there are multiple treatments entered (there can be only
        % one per data set).
        for i_d = 1:nr_sets % run through all data sets
            % Apply the correction factor to the exposure concentrations (not
            % first column, as those are the time points).
            data{i_d}.scen(:,2) = data{i_d}.scen(:,2:end) * corr_factor;
            
            % Create dummy survival data (this helps to be able to use the rest
            % of the code in this function). Just 10 survivors at t=0 and at
            % end of exposure scenario.
            data{i_d}.surv = [0 10 ; data{i_d}.scen(end,1) 10];
        end
end

%% BLOCK 2. Some general things to put in the output structures

data_all.nr_sets = nr_sets; % store number of data sets in the structure
data_all.hb_hat  = -1; % initialise background hazard rate from controls as -1

% Also prepare a structure for the control data, which will be fitted later
data_hb.nr_sets  = nr_sets; % store number of data sets in the structure
data_hb.nr_treat = zeros(1,nr_sets); % this will store number of control treatments per data set in this data set (vector with 0 or 1)
flag_ctrl_exist  = 0; % this flag will be set to 1 when there is at least one control data set
flag_ctrl_mort   = 0; % this flag will be set to 1 when there is at least one control with some mortality

%% BLOCK 3. Prepare object with survival data and exposure scenarios
% Code below runs through the different data sets. It creates a complex
% structure <data_all> that contains all information for all data sets
% entered, including the exposure scenarios. The data matrices that are
% entered are split up in the individual treatments, and several columns
% are added with information derived from the data.

for i_d = 1:nr_sets % run through all data sets
    
    % BLOCK 3.1. Initial things.
    nr_treat = length(data{i_d}.name);  % number of treatments in this data set
    % Extract info from this survivors data set
    t_surv   = data{i_d}.surv(:,1);     % complete time vector from survivor data set (vector)
    s_surv   = data{i_d}.surv(:,2:end); % survivors in the data set (matrix)
    % Extract info from this set of exposure scenarios
    t_scen   = data{i_d}.scen(:,1);     % time vector from exposure profiles (vector)
    s_scen   = data{i_d}.scen(:,2:end); % exposure concentrations in the data set (matrix)
    
    % Store some info for this data set already in the object: the list of
    % identifiers, the time vector, and the number of treatments.
    data_all.name(i_d,1:nr_treat) = data{i_d}.name; % copy names to <data_all> (cell arrays of strings)
    data_all.time{i_d}            = t_surv; % store entire time vector of this data set
    data_all.nr_treat(i_d)        = nr_treat; % store number of treatments in this data set
    
    % BLOCK 3.2. Preliminary checks on the input data. Producing errors
    % when things are wrong. Note that <load_data> has already performed
    % a number of checks itself!
    
    % Start with checking the survivor block.
    if t_surv(1) ~= 0 % first time point must be a zero
        error_str = cat(1,error_str,sprintf('Data set %1.0f: first time point in your survivor data must be zero',i_d));
    end
    if any(t_surv(1) < 0) % times must be positive
        error_str = cat(1,error_str,sprintf('Data set %1.0f: time points in the survivor data must be positive',i_d));
    end
    if any(isnan(s_surv(1,:))) % first row cannot contain NaNs
        error_str = cat(1,error_str,sprintf('Data set %1.0f: no missing data allowed at t=0 in survivor data',i_d));
    end
    if any(diff(t_surv) <= 0) % make sure time vector is okay
        error_str = cat(1,error_str,sprintf('Data set %1.0f: time vector in survivor data must be unique and increasing',i_d));
    end
    s_test                = s_surv(:); % turn the survivor matrix into a long column vector
    s_test(isnan(s_test)) = []; % remove any NaNs in there
    if any((s_test == floor(s_test)) ~= 1) % and look for non-integer survivors
        error_str = cat(1,error_str,sprintf('Data set %1.0f: number of survivors should all be whole numbers',i_d));
    end
    if any(s_surv(:)<0) % look for negative survivors
        error_str = cat(1,error_str,sprintf('Data set %1.0f: number of survivors should all be positive numbers',i_d));
    end
    
    % Next, check the exposure scenario block.
    if any(diff(t_scen) < 0) % make sure time vector is okay
        error_str = cat(1,error_str,sprintf('Data set %1.0f: time vector for exposure scenarios cannot decrease',i_d));
    end
    if any(t_scen(1) < 0) % times must be positive
        error_str = cat(1,error_str,sprintf('Data set %1.0f: time points in exposure scenarios must be positive',i_d));
    end
    if any(isnan(s_scen(1,:))) % first row cannot contain NaNs
        error_str = cat(1,error_str,sprintf('Data set %1.0f: no missing data allowed at t=0 in exposure scenario',i_d));
    end
    if any(diff(find(diff(t_scen)==0))==1) % also need to catch case with more than two of the same time points after another
        error_str = cat(1,error_str,sprintf('Data set %1.0f: include no more than two entries in exposure scenario at same time point',i_d));
    end
    if t_scen(1) ~= 0 % first time point must be a zero
        error_str = cat(1,error_str,sprintf('Data set %1.0f: first time point in your exposure scenario must be zero',i_d));
    end
    if any(s_scen(:)<0) % look for negative concentrations
        error_str = cat(1,error_str,sprintf('Data set %1.0f: exposure concentrations should be positive numbers',i_d));
    end
        
    if ~isempty(error_str)
        if i_d < nr_sets % add a message to clarify what has been checked (we return before ALL data are checked)
            error_str = cat(1,error_str,sprintf('Checked until data set %1.0f (there may be more errors)',i_d));
        end
        return % don't continue but return to the function that is calling
    end
    
    % BLOCK 3.3. Split up the survivor data per treatment.
    for i_c = 1:nr_treat % run through treatments in data set
        
        % BLOCK 3.3.1. Remove NaNs, add columns with number of deaths and
        % survival probabilities.
        s_ic = s_surv(:,i_c); % vector with survivors for this treatment only
        t_ic = t_surv;        % copy time vector as well (as we'll remove NaNs)
        index_nan = isnan(s_ic); % indices to NaNs in survivor vector
        s_ic(index_nan) = []; % remove time points with missing data for survivors
        t_ic(index_nan) = []; % remove time points with missing data in time vector
        d_ic  = [-diff(s_ic) ; s_ic(end)]; % vector with deaths per interval (incl. the interval after the test)
        ps_ic = s_ic / s_ic(1); % vector with survival probabilities
        
        if any(d_ic<0) % make sure there are no zombies in there ...
            error_str = cat(1,error_str,sprintf('Data set %1.0f, %1s: number of survivors should never increase in time',i_d,data{i_d}.name{i_c}));
            if i_c < nr_treat % add a message to clarify what has been checked (we return before ALL data are checked)
                error_str = cat(1,error_str,sprintf('Checked until data set %1.0f, %1s (there may be more errors)',i_d,data{i_d}.name{i_c}));
            end
            return % don't continue but return to the function that is calling
        end
        
        % BLOCK 3.3.2. Wilson score interval on data probabilities.
        % https://en.wikipedia.org/wiki/Binomial_proportion_confidence_interval#Wilson_score_interval
        z = 1.96; % critical value of normal distribution for 95% CI
        n = s_ic(1); % number of animals to start with
        a = (ps_ic + z^2/(2*n))/(1+z^2/n); % middle of confidence interval (element-wise calculation on the vector <ps_ic>)
        b = z/(1+z^2/n) * sqrt(ps_ic.*(1-ps_ic)/n + z^2/(4*n^2)); % range of interval (element-wise calculation on the vector <ps_ic>)
        b(1) = 0; % at t=0, there is no uncertainty, so remove the range
        a(1) = 1; % and make middle of range 1
        % Note: the interval itself (the min and max) are given by a-b and
        % a+b, which is done in the next block.
        
        % BLOCK 3.3.3. Combine all vectors into the overall object, as a
        % nice matrix for data set <i_d> and treatment <i_c>. columns:
        % time, survivors, nr. deaths, surv. prob., CI min, CI max
        data_all.surv{i_d,i_c} = [t_ic s_ic d_ic ps_ic max(0,a-b) min(1,a+b)]; % creates a matrix in a cell of <data_all.surv>
        
    end
    
    % BLOCK 3.4. Split up the exposure scenarios per treatment. Note: I run
    % through the treatments again here. It could be done in the same loop
    % as the survival data, but a separate loop may be more transparent.
    % Note that <load_data> already performed checks whether all survival
    % treatments are unique, and whether they all have matching exposure
    % scenarios (and places these scenarios in the same order, if needed).
    
    c_twa = zeros(1,nr_treat); % initialise vector to collect TWA concentrations
    c_max = zeros(1,nr_treat); % initialise vector to collect maximum concentrations
    
    % BLOCK 3.4.1. Constant exposure can be entered as a single row. We'll
    % add a final time point, so it will run smoothly through the same
    % piece of code as for time-varying concentrations.
    if size(s_scen,1) == 1 % if it is only a single row (constant)
        % Add final time point with same concentrations.
        s_scen = [s_scen ; s_scen]; %#ok<AGROW> copy the scenario row vector to a second row
        t_scen = [t_scen ; t_surv(end)]; %#ok<AGROW> and add a final time point to the time vector
    end
    
    for i_c = 1:nr_treat % run through treatments of this data set again
        
        % BLOCK 3.4.2. For each scenario, create an events matrix <ev_ic>
        % with places where the linear patterns changes, and a slope in the
        % third column. Start from the time-conc matrix for this scenario,
        % but exclude last time point (things stay constant after last time
        % point). Third column has the linear slopes in each interval.

        ev_ic = [t_scen s_scen(:,i_c)]; % collect time and concentrations for this treatment (matrix)
        ev_ic(isnan(ev_ic(:,2)),:) = []; % remove rows with NaNs in the concentration column
        ev_ic(1:end-1,3) = diff(ev_ic(:,2)) ./ diff(ev_ic(:,1)); % add slopes in third column (last slope is automatically zero)
        
        ev_ic(1+find(diff(ev_ic(1:end-1,3))==0),:) = []; % prune it: remove values were slope remains the same! (but keep last interval)
        % Note (TJ 10/9/2020): I cannot remember why I placed the 1:end-1
        % in there. Everything seems to run smoothly when kicking the last
        % interval out, if the slope remains the same. However, there is no
        % difference in <ev_ic>, as far as I can tell, so I will leave it
        % as is.
        %
        % Note: pruning needs to be done before removing any INFs, because
        % we allow multiple entries at the same time point. We thus can
        % have two intervals with the same slope that are not continuously
        % increasing. Pruning is thereby not always 100% effective when
        % there are such multiple timepoints (e.g., in the controls, there
        % may be multiple time points, but the slope and concentration
        % remain at zero; these are not recognised as 'prunable'). I don't
        % think this leads to serious inefficiency in the calculations.
        
        % Take care of cases where two instances of the same time point are
        % entered (which signals that we want to have an instantaneous
        % change in conc.). That gives a slope of infinity.
        ev_ic(isinf(ev_ic(:,3)),:) = []; % simply remove the rows with slope INF
        ev_ic(isnan(ev_ic(:,3)),:) = []; % simply remove the rows with slope NaN (might occur in some silly cases, when the same concentration is entered as double time points)
        
        % Make sure that the exposure scenario nicely matches the survival data.
        if t_surv(end) < ev_ic(end,1) % does the exposure scenario extend beyond the time vector of the survival data?
            ev_ic(ev_ic(:,1)>t_surv(end),:) = []; % remove any extra events
            % This could well leave <ev_ic> short, but that is fixed with the next step
        end
        if t_surv(end) > ev_ic(end,1) % do we need more points for survival than available in the exposure scenario?
            ev_ic = [ev_ic; t_surv(end) ev_ic(end,2)+ev_ic(end,3)*(t_surv(end)-ev_ic(end,1)) 0]; %#ok<AGROW>
            % Add an extra time point in <ev_ic>, extrapolated from last
            % one, and add slope zero. Note: when the scenario was short
            % already at the input, the final slope is zero, so there is no
            % extrapolation. There is only extrapolation with a non-zero
            % slope if we cut off the exposure scenario in the previous if
            % statement.
        end
        
        % BLOCK 3.4.3. Calculate the maximum and TWA concentration for this
        % treatment. This is useful for scaling plots and for calculating
        % the parameter ranges in <startgrid>.
        c_twa(i_c) = 0; % initialise TWA as zero (as we'll continuously add to it)
        c_max(i_c) = ev_ic(1,2); % start at start concentration for this treatment
        for k = 1:size(ev_ic,1)-1 % run through all time intervals
            delta_t    = ev_ic(k+1,1) - ev_ic(k,1); % length of time interval
            c_start    = ev_ic(k,2); % start concentration in this interval
            slope      = ev_ic(k,3); % slope in this interval
            c_end      = c_start + slope * delta_t; % concentration at end of this interval
            c_twa(i_c) = c_twa(i_c) + mean([c_start c_end]) * delta_t / ev_ic(end,1);
            % Add contribution of this interval to the TWA. Contribution is
            % mean of start and end concentration in the interval, weighted
            % with the length of the interval, as fraction of the total
            % length of the scenario.
            c_max(i_c) = max(c_max(i_c),c_end); % this will be the max exposure concentration across all time points
            % This is done here to make sure that we catch the max (which
            % might be higher than one of the entries of ev_ic due to the
            % possibility of two concentrations at the same time points
            % where the latter is removed).
        end
        
        data_all.scen{i_d,i_c}  = ev_ic; % store the events-time matrix in the overall object
        
    end
    
    data_all.c_twa{i_d} = c_twa; % store the time-weighted averages over all treatments
    data_all.c_max{i_d} = c_max; % store the maximum concnetration over all treatments
    
    % BLOCK 3.4. Also collect the survival data set for the control (where
    % the scenario identifier is 'Control'). This fills a new data object <data_hb>
    % for the control data so it can be easily fitted later on.
    index_ctrl = strcmpi(data_all.name(i_d,:),'Control'); % find where the control is in this data set (logical index)
    % Note: <strcmpi> compares strings without caring for case.
    if sum(index_ctrl) == 1 % if there is a control ... (there can be only one)
        data_hb.surv{i_d,1}   = data_all.surv{i_d,index_ctrl}; % this is its extended survivor data matrix for the control
        data_hb.nr_treat(i_d) = 1; % the control only has one treatment
        data_hb.time{i_d}     = t_surv; % store entire time vector of this data set
        data_hb.scen{i_d,1}   = [0 0 0; t_surv(end) 0 0]; % create and store an events-time matrix (all zeros)
        data_hb.name{i_d}     = 'Control'; % add a name
        flag_ctrl_exist       = 1; % flag that there is at least one control data set
        if data_hb.surv{i_d,1}(end,2) < data_hb.surv{i_d,1}(1,2) % check if there is mortality in this control
            flag_ctrl_mort    = 1; % flag that there is control mortality, so we need fitting
        end
    end
    
end

%% BLOCK 4. Fit background mortality <hb> on all controls
% The structured object <data_hb> contains a sub-data set for the controls
% only. When there are no controls at all, the entry in the data object
% <data_all.hb_hat> will remain at -1, and <startgrid> will produce an
% error when the user requests fixing the background hazard to the
% controls.

if flag_ctrl_exist == 1 % if at least one control data set is located ...
    
    hb_min = 1e-6; % minimum value for hb (best have same as in startgrid!)
    
    if flag_ctrl_mort == 0 % if there is no mortality in the controls at all, there is no need to fit
        data_all.hb_hat = hb_min; % store the minimum <hb> in the structure <data_all>
        % This avoids having a weird value for <hb> such as 2.476e-6 where
        % 1e-6 would fit just as well.
    else % we need to fit!
                
        pmat           = zeros(5,5); % make a fake <pmat> with zeros
        pmat(:,end)    = 1; % but all on normal scale
        pmat(GLO.hb,:) = [0.01 1 hb_min 0.1 1]; % prepare a row in <pmat> for <hb> (0.1 should be high enough as max)
        pfit           = pmat(GLO.hb,1); % parameter vector to be fitted (only <hb>)
        
        % Use <setup_simplex> and the standard <transfer> function to fit the control data.
        [hb_hat,~] = setup_simplex(pfit,1,pmat,1,data_hb); % do a normal optimisation (use SD)
        % NOTE: second, optional, output is fval (the MLL) which could be
        % output to screen if needed in the future.
        data_all.hb_hat = hb_hat; % store the best-fitting <hb> in the structure <data_all>
        
    end
end

