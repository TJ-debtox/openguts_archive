function x_out = calc_model(t,pmat,ev_ic,sdit)

% Usage: x_out = calc_model(t,pmat,ev_ic,sdit)
%
% This function returns the state variables survival probability and scaled
% damage for a given time vector, parameter set, and exposure scenario
% (time-events matrix). The differential equation for scaled damage is
% always solved analytically (done in a separate function
% <calc_damage_linear>). Mortality due to chemical stress is calculated in
% separate functions, which are also used for the calculation of LCx and
% LPx. For IT and for SD with constant exposure, survival probabilities are
% calculated analytically. For SD with time-varying exposure, hazard rates
% are calculated on a fine time scale and integrated using the trapezium
% rule. 
%
% Inputs
% <t>      the time vector for which output needs to be returned
% <pmat>   the parameter matrix (or vector with parameter values)
% <ev_ic>  the exposure scenario for a single treatment as time-events matrix (three columns)
% <sdit>   switch for SD (1) or IT (2)
%
% Outputs
% <x_out>  matrix with time in rows, and states in columns (survival probability and scaled damage)
%
% Author     : Tjalling Jager
% Date       : May 2020
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

global GLO % allow for global parameters in structure <GLO>

%% BLOCK 1. Initial settings
% We need to consider the time vector for the calculations. For constant
% exposure, we use the one provided in the input <t>. For time-varying
% exposure, it depends on the death mechanism. For IT, we'll add the
% scenario time points in <ev_ic> to the time vector. For SD, we'll make a
% very detailed time vector as we need to numerically integrate hazard
% later.

t     = t(:); % force <t> to be a row vector (it should be already)
t_rem = t;    % remember the original time vector (so we can extract the ones we need as output later)

if all(diff(ev_ic(:,2))==0) && all(ev_ic(:,3)==0) % if all concentrations are equal and all slopes are zero
    % Then we have a constant exposure concentration, and don't need to do
    % anything special to the time vector.
    constant = 1; % set flag that we have constant exposure
else
    constant = 0; % set flag for constant exposure to off

    % =====================================================================
    % NEW after v.1.0: it is safer to ALWAYS use a detailed time vector for
    % time-varying exposure. In v1 and earlier, I assumed that for IT,
    % maximum damage will always occur at a time point that is in the
    % events definition <ev_ic>, but this is not correct.
    
    t = unique([t;ev_ic(:,1)]); % (vector)
    % this adds the time vector of the exposure scenario to the existing
    % vector <t>. This is generally a good idea.
    %     
    % Note: this is a bit awkward when we are called from <calc_lpx>, as
    % the profile <ev_ic> will already be in the <t> that we are called
    % with. So, we'll be adding points and removing them in the same line
    % of code. We could try to speed things up by calling <calc_model> with
    % a flag to signal we're calling from <calc_lpx>. 
    
    % For time-varying exposure, we need to make a longer time vector. This
    % could perhaps be done in <prepare_data>, but since we also use this
    % code for providing model output for plotting, this is the most
    % transparent way to implement it (though probably not the fastest).

    min_t = t(end) * GLO.t_per_day; % minimum length of time vector for calculation (default each 15 minutes)
    % The same global default is used in <calc_lpx>.
    if length(t) < min_t % make sure there are at least <min_t> points
        t = unique([t;(linspace(min(t),max(t),min_t))']); % (vector)
        % this adds a regular-spaced <min_t> time points to the
        % existing vector <t> (the originals still need to be in there,
        % as we need them for the calculation of the likelihood in
        % <transfer>). For LPx calculations, this part will automatically
        % be skipped as the time vector will made be on this detailed
        % basis already in <calc_lpx> (and thus, <t> will not be
        % smaller than <min_t>).
    end
    % =====================================================================
    
end

%% BLOCK 2. Calculate damage
% This part calls the explicit model in <calc_damage_linear> to calculate
% the damage over time. However, if we have constant exposure, we can
% calculate it directly here, which saves a little bit of time.

kd = pmat(GLO.kd,1); % extract the dominant rate constant
if constant == 1 % for constant concentration, no need to call the function
    Dw = ev_ic(1,2) * (1-exp(-kd*t)); % first-order accumulation (vector)
else    
    Dw = calc_damage_linear(t,kd,ev_ic); % call the function to calculate damage <Dw> (vector)
end

%% BLOCK 3. Calculate survival
% We always use the analytical solution for the damage level in
% <simple_linear>. Therefore, only the damage level is available at this
% point, and the survival probabilities are calculated from it for each
% death mechanism separately. Both survival and scaled damage are collected
% in a single matrix that will be the output of this function.
%
% Note: the calls of these functions below are different from each other.
% This is a consequence of the fact that they have to serve various (and
% different) purposes. For SD, the functions need to be used in conjunction
% with <fzero> (root finding) to find LCx or LPx. For IT, these metrics are
% calculated directly.

hb = pmat(GLO.hb,1); % extract background hazard rate from <pmat>

if sdit == 1 % stochastic death
        
        if constant == 1 % for constant exposure, we can use the analytical solution for SD
            [~,Sc] = calc_surv_sd_c(ev_ic(1,2),pmat,t,[]); % only keep survival (<Sc>, vector)
        else % approximate solution with trapezium rule, calculated from <Dw> in a different function
            [~,Sc] = calc_surv_sd_t(1,pmat,t,Dw,[]); % only keep survival (<Sc>, vector)
        end
        
else % individual tolerance, log-logistic distribution is used for GUTS-IT
        
        % Survival under IT at time <t> is determined by the maximum of
        % <Dw> up to that time point. Therefore, we modify <Dw> to a
        % <maxDw> vector that does not decrease over time (dead animals
        % don't become alive). This algorithm is a bit awkward, but pretty
        % fast in Matlab, especially when there are few episodes with
        % decrease in <Dw>.
        
        maxDw = Dw; % copy the vector <Dw> to <maxDw>
        ind   = find([0;diff(maxDw)]<0,1,'first'); % first index to places where <Dw> has decreased
        while ~isempty(ind) % as long as there is a decrease somewhere ...
            maxDw(ind:end) = max(maxDw(ind:end),maxDw(ind-1));  % replace every entry at later time with max of that and previous point
            ind            = find([0;diff(maxDw)]<0,1,'first'); % any decrease left? Then we do this again
        end
        % Note: <maxDw> is a vector with the maximum damage till that
        % point, at each time point.
        
        Sc = calc_surv_it(maxDw,pmat,[],[],[]); % call function to calculate IT survival from <maxDw>
end

S     = Sc .* exp(-hb*t); % multiply survival due to the chemical with survival resulting from background hazard (element-wise: vector)
x_out = [S Dw];           % put both state variables in the output matrix as two columns
% Extract the time points that we need from the total time vector.
if length(t) > length(t_rem)      % only needed if we calculated more points than we were asked to return (not strictly needed, but might be faster)
    [~,loct] = ismember(t_rem,t); % find where the requested time points are in the longer <t> (vector with indices)
    x_out    = x_out(loct,:);     % only keep the rows we asked for
end

