function [Sc,LCx,LPx] = calc_surv_it(maxDw,pmat,t,Feff,mf_try)

% Usage: [Sc,LCx,LPx] = calc_surv_it(maxDw,pmat,t,Feff,mf_try)
% 
% This function calculates survival probability <Sc> for the individual
% tolerance model (can be a vector). Additionally, it provides a direct
% calculation for LCx and LPx.
%
% Inputs
% <maxDw>   the maximum damage over time (can be a vector with moving max over time or single value)
% <pmat>    the total parameter matrix (or vector with first column only)
% <t>       the time at which to calculate the output (only used for LCx where it is a single value)
% <Feff>    the fraction effect (x/100 in LCx and LPx)
% <mf_try>  the multiplication factor to use to calculate survival (for use with <calc_lpx>)
% 
% Outputs
% <Sc>      survival probability due to chemical stress for each element of <t>
% <LCx>     the LCx, calculated directly from the parameters in <pmat>
% <LPx>     the LPx, calculated directly from the parameters in <pmat> and <maxDw>
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
mw   = pmat(GLO.mw,1);  % median of threshold distribution
mw   = max(mw,1e-100);  % make sure that the threshold is not exactly zero ... 
% Note: at this moment that is not needed as we set a non-zero lower bound for <mw>.
Fs   = pmat(GLO.Fs,1);  % fraction spread of the threshold distribution
Fs   = max(1+1e-6,Fs);  % make sure it is always slightly above 1
beta = log(39)/log(Fs); % shape parameter for log-logistic from <Fs>
kd   = pmat(GLO.kd,1);  % dominant rate constant

if isempty(mf_try) % if it is empty ...
    mf_try = 1; % use 1 (which implies no multiplication of <maxDw>)
end

Sc = 1 ./ (1+(mf_try*maxDw/mw).^beta); % survival probability for log-logistic distribution
% Note: in this calculation, <maxDw> can be a vector (then, <Sc> will also
% be a vector). Also note that <mf_try> is in here as <calc_lpx> also wants
% to calculate survival at a range of multiplication factors.

if isempty(Feff) % no need to calculate anything more when <Feff> is empty (that may lead to error)
    LCx = -1;
    LPx = -1;    
else % I am pretty sure that in the next calculations, all variables are always single values
    LCx = (mw/(1-exp(-kd*t))) * (Feff/(1-Feff))^(1/beta); % calculate LCx
    LPx = (mw/maxDw) * (Feff/(1-Feff))^(1/beta);          % calculate LPx
end
