function minloglik = transfer(pfit,pmat,sdit,datacal_input)

% Usage: minloglik = transfer(pfit,pmat,sdit,datacal_input)
%
% This is a transfer function which calculates the minus log-likelihood
% function for a given set of parameters. This function is called by
% <calc_parspace> and by <prepare_data> (for background hazard rate, via
% <fminsearch>).
%
% This transfer function allows certain parameters to be fixed, while
% others are fitted. The free parameters are in <pfit>, all parameters in
% <pmat>, which contains in the second column a 1 where the parameter is to
% be fitted. This format is required to work with the Matlab optimisation
% routines (specifically <fminsearch>). Note that parameters in <pfit> and
% in <pmat> are log-transformed where the fifth column of <pmat> is a 0
% (including the parameters in <pmat> that are NOT fitted).
% 
% Inputs
% <pfit>          vector with to-be-fitted parameter values
% <pmat>          the parameter matrix
% <sdit>          switch for SD (1) or IT (2)
% <datacal_input> input data set, ONLY used for calculating background
%                 hazard rate in <prepare_data>
% 
% Outputs
% <minloglik>     minus log-likelihood associated with values in <pfit>
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

global DATAcal % make calibration data set global

if isempty(datacal_input) % this input will be used when fitting background mortality; in that case, <transfer> will receive a dataset ...
    datacal = DATAcal; % use the global with the calibration data set
else
    datacal = datacal_input; % otherwise, use the data set in the input (for background hazard rate)
end

%% BLOCK 1. First some work on the parameter matrix

p_tmp          = pmat(:,1);      % copy the first column of the parameter matrix to a temporary parameter vector
ind_fit        = pmat(:,2) == 1; % logical indices to fitted parameters (vector)
ind_log        = pmat(:,5) == 0; % logical indices to log-scale parameters (vector)
p_tmp(ind_fit) = pfit ;          % replace to-be-fitted values of the temporary vector with the values in <pfit>

% Note: the code below checks whether the parameter is within bounds on the
% *log-transformed* value (rather than transforming the parameter itself
% and checking with the normal bounds). Reason is that 10^(log10(x))-x is
% not equal to x, but a small number (related to machine precision). This
% means that when an element of <p_tmp> would be exactly on the bound, it
% would still be seen as outside. This approach works better.

% Put parameter bounds for parameters fitted on log scale on log scale.
pmat(ind_log,[3 4]) = log10(pmat(ind_log,[3 4]));
% If any parameter is outside the bounds, return +inf as criterion.
if any(p_tmp<pmat(:,3) | p_tmp>pmat(:,4))
    minloglik = +inf; % than make the <minloglik> infinite (very bad fit)
    return; % and return to the function calling this one (no need to run model)
end
% Put parameters that need to be fitted on log scale back on normal scale.
p_tmp(ind_log) = 10.^(p_tmp(ind_log));
% Note: <pmat> needs to be on normal scale for <calc_model>.

%% BLOCK 2. Get model results
% This block will call <calc_model> for the model predictions, and together
% with the data set, turn that into a log-likelihood value. It has two
% for-loops to loop over the data sets and over the individual treatments
% in each data set. 

nr_treat = datacal.nr_treat; % vector with the number of scenarios in each data set
loglik   = zeros(sum(nr_treat),1); % initialise vector to catch all log-likelihoods
cntr     = 0; % simple counter for number of data sets in total (as index for <loglik>)

for i_d = 1:datacal.nr_sets % run through all data sets
    
    for i_c = 1:nr_treat(i_d) % run through all exposure scenarios in this data set
        
        % Extract the data block for survivors and exposure scenario for this treatment only.
        data_surv = datacal.surv{i_d,i_c}; % survival info from data (matrix)
        ev_ic     = datacal.scen{i_d,i_c}; % events matrix (three columns)
        cntr      = cntr+1; % increase the simple counter (index for <loglik>)
        
        if ~isempty(data_surv) % when the data set is not empty ...
            % I think <data_surv> will be empty when fitting control data
            % on multiple data sets, where a control is missing in some
            % data sets. Since <loglik> is initialised with zeros, the
            % loglik on this position will remain zero, which is fine.
            
            x_out = calc_model(data_surv(:,1),p_tmp,ev_ic,sdit); % use <calc_model> to provide the output for this scenario
            % Note: <calc_model> is called with the first column of
            % <data_surv>, which contains the time vector for this block
            % (NaNs in the data set have been removed). The <x_out> will thus
            % only contain the predictions at the time points where there
            % is data.
            
            deathprob = [-diff(x_out(:,1));x_out(end,1)]; % multinomial death probabilities from model for each interval (incl. interval after the test!)
            deathdata = data_surv(:,3); % extract number of deaths from the data set
            deathprob = max(deathprob,1e-50); % otherwise we get problems taking the logarithm
            % Note: this 1e-50 is unlikely to come into play as we set a
            % non-zero minimum for <hb>. Only for extremely short intervals
            % between observations may this be needed to avoid an error.
            
            loglik(cntr) = deathdata'*log(deathprob);  % and calculate the log likelihood
            % Note: this is vector product. What it does is element-wise
            % product of <deathdata> and ln-<deathprob>, and sum the
            % individual products.
        end
    end
end

%% BLOCK 3. Final things

minloglik = -1 * sum(loglik) ; % sums the log-likelihoods and makes it negative 
% Note: this is needed as <fminsearch> is looking for a MINIMUM, and we
% need to find the MAXIMUM of the likelihood function.

if isnan(minloglik) || ~isreal(minloglik) % when the loglik calculation returns NaN or complex numbers ...
    minloglik = +inf; % give it a really bad likelihood value so the optimisation ignores it
    % I am not sure if this still occurs, but it does not hurt to check and
    % return a +inf.
end

