function [crit,Sc] = calc_surv_sd_t(mf_try,pmat,t,Dw,Feff)

% Usage: [crit,Sc] = calc_surv_sd_t(mf_try,pmat,t,Dw,Feff)
% 
% This function calculates survival probability <Sc> under stochastic death
% when exposure is time-varying, using numerical integration with the
% trapezium rule. Additionally, it provides a criterion <crit> to be used
% by <fzero> to calculate the LPx. LPx is defined as the factor by which to
% multiply an exposure profile to obtain x% effect at the end of the
% profile. I here make use of the fact that multiplying a profile by a
% factor is equivalent to multiplying damage levels by that factor. This
% works as long as damage is a simple one-compartment thing! The <mf_try>
% needs to be the first input parameter as that is what <fzero> will be
% able to change. Also for <fzero>, <crit> needs to be the first output.
%
% Inputs
% <mf_try>    the multiplication factor to try
% <pmat>      the parameter matrix (or only the first column)
% <t>         the time vector at which to calculate the output (should match <Dw>)
% <Dw>        damage (as vector over time <t>)
% <Feff>      the fraction effect (x/100 in LPx); when not empty, only
%             return <Sc> at the end of the time vector (LPx calculations)
% 
% Outputs
% <crit>      criterion that is zero when survival at end is 1-x/100
% <Sc>        survival probability over time <t>, due to the chemical
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

global GLO % make global the structure <GLO>

% Extract parameters from <pmat>.
mw = pmat(GLO.mw,1); % median of threshold distribution
bw = pmat(GLO.bw,1); % killing rate

Dw     = mf_try * Dw;       % use multiplication on damage
hc     = bw * max(0,Dw-mw); % calculate hazard for each time point
cum_hc = cumtrapz(t,hc);    % integrate the hazard rate numerically over the time vector
Sc     = min(1,exp(-1*cum_hc)); % calculate survival probability due to chemical stress

if isempty(Feff)
    crit = -1; % just return a -1; no need to calculate anything else but <Sc>
else % when a <Feff> is entered ...
    Sc   = Sc(end);       % only return survival at the end of the exposure profile
    crit = Sc - (1-Feff); % calculate criterion: zero when end value is x*100% effect
end

