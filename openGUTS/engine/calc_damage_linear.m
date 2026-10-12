function Dw = calc_damage_linear(t,kd,ev_ic)

% Usage: Dw = calc_damage_linear(t,kd,ev_ic)
%
% This function calculates scaled damage for the reduced GUTS model system,
% given a series of linearly-changing exposure concentrations. Note that
% the survival probability due to chemical stress is calculated separately.
% This function sequentially applies the analytical solution for scaled
% damage, for each episode between events. This function is only called
% from <calc_model>. 
% 
% Inputs
% <t>      the time vector
% <kd>     the dominant rate constant
% <ev_ic>  the exposure-events matrix for a single treatment (three columns)
%
% Output
% <Dw>     vector with scaled damage at each value of <t>
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

%% Calculate the model output
% This is the actual damage-dynamics model, specified as explicit function
% that is called sequentially for each episode in <ev_ic>.

Dw    = nan(length(t),1); % initialise output vector for damage with NaNs
Dw0   = 0;   % make initial damage level zero
Dw(1) = Dw0; % and put it in the first spot of the output vector

for i = 1:size(ev_ic,1)-1 % run through all events (not last one, which marks the end of the scenario, to which we extrapolate with the slope)

    ind_t   = (t>ev_ic(i,1) & t<=ev_ic(i+1,1)); % find the logical indices for the period between two events
    te      = t(ind_t) - ev_ic(i,1); % find the new part of the time vector, and make it start at zero again
    c_start = ev_ic(i,2); % take as initial concentration the new one in this period
    slope   = ev_ic(i,3); % take as slope the new one in this period
    
    % Thanks to Bob Kooi for using Maple to obtain the solution for me.   
    Dw(ind_t) = slope * te + c_start - slope / kd + exp(-kd * te) * (Dw0 - c_start + slope / kd);
    
    % Find new starting values at exact moment of start of new period.
    t_end = ev_ic(i+1,1) - ev_ic(i,1); % last time point in this interval is start time for new period
    Dw0   = slope * t_end + c_start - slope / kd + exp(-kd * t_end) * (Dw0 - c_start + slope / kd);
    % No need to do this in the last iteration of the loop, but not sure if
    % adding an <if> statement would really be useful.
end


