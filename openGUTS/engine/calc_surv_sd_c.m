function [crit,Sc] = calc_surv_sd_c(c,pmat,t,Feff)

% Usage: [crit,Sc] = calc_surv_sd_c(c,pmat,t,Feff)
% 
% This function calculates survival probability <Sc> under stochastic death
% when exposure is constant. Additionally, it provides a criterion <crit>
% to be used by fzero to calculate the LCx for SD. As LCx is defines as x%
% effect under constant exposure, we can use the analytical solution here
% as well (and ignore background hazard). The <c> needs to be the first
% input parameter as that is what <fzero> will be able to change. Also for
% <fzero>, <crit> needs to be the first output.
%
% Inputs:
% <c>       the concentration at which to calculate the output
% <pmat>    the parameter matrix (or first column only)
% <t>       the time at which to calculate the output (can be a vector)
% <Feff>    the fraction effect (x/100 in LCx)
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
kd = pmat(GLO.kd,1); % dominant rate constant

% For constant exposure, we can use analytical solutions for SD.
minH = zeros(length(t),1); % this will be the minus integrated hazard rate if nothing happens (vector)
if c > mw % if concentration to test is above the threshold ...
    t0   = -log(1-mw/c)/kd; % no-effect-time: time point where effects start
    minH = (bw/kd)*max(0,exp(-kd*t0)-exp(-kd*t))*c - bw*(max(0,c-mw))*max(0,t-t0); % replace <minH> by the calculated effects over time (vector)
end
Sc = min(1,exp(minH)); % turn <minH> into survival probability (vector)

if isempty(Feff)
    crit = -1; % no need to calculate anything else than <Sc>
else
    crit = Sc - (1-Feff); % calculate criterion: zero when end value is x*100% effect (for <fzero> to work on)
end

