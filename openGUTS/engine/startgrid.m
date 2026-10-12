function pmat = startgrid(sdit)

% Usage: pmat = startgrid(sdit)
%
% This function prepares a parameter structure for GUTS analyses, with
% min-max ranges, as required for the parameter-space explorer
% (<calc_parspace>). As input, it needs to know the setting for SD or IT in
% <sdit>.
%
% Note: at this moment, this function does not allow recalculation of
% bounds when the user has fixed some parameters (or set other bounds for
% some parameters). Since some bounds depend on the value of other bounds,
% this could be reconsidered. In any case, the code requires only little
% modification to allow a <pmat> to be entered by the user.
% 
% Inputs
% <sdit>  switch for SD or IT calculation (they have different ranges)
% 
% Outputs
% <pmat>  matrix with parameter info (best, fit/fix, min, max, log/norm)
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

global GLO DATAcal % allow for global parameters in structure glo and data set

%% BLOCK 1. Preliminary things 
% First some preliminary checks.

if ~isfield(GLO,'concunit') % if the concentration unit is not defined (this should not happen anymore with the standard file format)
    GLO.concunit = '[C]'; % set a default
end

% Defining the basic structure of the <pmat> matrix.
pmat(GLO.kd,:) = [0 1 0 0 0]; % dominant rate constant, (1/d)
pmat(GLO.mw,:) = [0 1 0 0 1]; % median threshold for survival ([C])
pmat(GLO.hb,:) = [0 1 0 0 1]; % background hazard rate (1/d)
pmat(GLO.bw,:) = [0 1 0 0 0]; % killing rate (SD only) (1/([C] d))
pmat(GLO.Fs,:) = [0 1 0 0 0]; % fraction spread of threshold distribution (IT only) (-)

% Make sure that right parameters are fitted for each special case.
if sdit == 1 % for SD ...
    pmat(GLO.Fs,[1 2 5]) = [1 0 1]; % never fit the threshold spread, set to 1, normal scale
else % for IT ...
    pmat(GLO.bw,[1 2 5]) = [+inf 0 1]; % never fit the killing rate, set to infinity, normal scale
end

if GLO.fix_hb_cal == 1 % if we want to fix hb to the control value (and it is available)
    if DATAcal.hb_hat == -1
        error('Requested fixing background hazard to control for calibration, but no controls were found.')
    else
        pmat(GLO.hb,[1 2 5]) = [DATAcal.hb_hat 0 1]; % never fit the background hazard, set to fitted value, normal scale
    end
end

% For parameters that are fixed, copy the entered value in the first column
% to both ranges. This may help to obtain smaller search ranges for some of
% the other parameters (when the range is calculated from another range),
% and it looks better when printing results on screen.
pmat(pmat(:,2)==0,3) = pmat(pmat(:,2)==0,1);
pmat(pmat(:,2)==0,4) = pmat(pmat(:,2)==0,1);

%% BLOCK 2. Find min-max time and exposure concentration across all data sets. 
% Make use of the info already stored for each data set. This could also be
% done in <prepare_data> ... however, since the rules for the starting
% ranges may be changed in the future, it is probably safer to keep them
% here.

t_max = 0;   % maximum time
t_min = inf; % minimum *non zero* time
c_max = 0;   % look for maximum peak concentration across all data sets
c_min = inf; % look for the minimum *non-zero* TWA across all data sets

for i_d = 1:DATAcal.nr_sets % run through all data sets for calibration
    c_twa = DATAcal.c_twa{i_d};                   % extract vector with TWAs for this data set
    c_min = min([c_min,min(c_twa(c_twa>0))]);     % look for the minimum *non-zero* TWA across all data sets
    c_max = max([c_max,max(DATAcal.c_max{i_d})]); % look for maximum peak concentration across all data sets
    
    t_max = max([t_max,DATAcal.time{i_d}(end)]); % look for maximum time point across all data sets
    t_min = min([t_min,DATAcal.time{i_d}(2)]);   % look for minimum *non zero* time point across all data sets
    % Since the time vector is unique, and the first entry zero, the second
    % element is the time for the first observation on survival.
end

%% BLOCK 3. Find relevant min-max ranges based on the data set.
% Rules of thumb are used to restrict parameter space, based on the data
% set.

if pmat(GLO.hb,2) == 1 % only if <hb> is fitted (otherwise, a value is provided in the input)
    pmat(GLO.hb,3:4) = [1e-6 , 0.07]; % simply use a fixed upper boundary for <hb> for now
    % Note: the minimum value is also set in <prepare_data> for fitting
    % <hb> on the controls. If we change the value here, it is best also
    % changed there!
end

if pmat(GLO.kd,2) == 1 % only if <kd> is fitted
    % In principle, fixed bounds for <kd> will be used. They are set at 95%
    % of steady state within 0.5 hour (resolution of FOCUS is 1 hour), and
    % 95% of steady state achieved in 5 years. However, if people do tests
    % with a first observation within hours, or run extremely long tests,
    % we need to extend them. The extension is based on the time vector:
    % only 5% of steady state achieved at the end of the test, and 99% of
    % steady state achieved within 1/10 of the first observation. In
    % practice, these seem to be the identifiability limits.
    
    pmat(GLO.kd,3) = min([log(20)/(5*365)  , -log(1-0.05)/t_max]);
    pmat(GLO.kd,4) = max([log(20)/(0.5/24) , -log(1-0.99)/(0.1*t_min)]);
    
end

if pmat(GLO.mw,2) == 1 % only if <mw> is fitted
    % Minimum is body residue that lowest conc. will achieve in 4 hours at lowest value of <kd>.
    pmat(GLO.mw,3) = c_min*(1-exp(-pmat(GLO.kd,3)*(4/24)));
    if sdit == 1 % for SD ...
        pmat(GLO.mw,4) = 0.99*c_max; % max as 99% of maximum peak concentration
    else % for IT ...
        pmat(GLO.mw,4) = 2*c_max; % max as two times of maximum peak concentration
    end
        
end

if pmat(GLO.bw,2) == 1 % only if <bw> is fitted
    % Min is 10% effect at end of test when fast kinetics and threshold
    % zero. Max is that what you gain in damage in one hour, at lowest
    % exposure, at minimum <kd>, is enough to kill you (95% sure) in 1 hour.
    pmat(GLO.bw,3) = -log(0.9)/(c_max*t_max);
    pmat(GLO.bw,4) = (24^2*0.95)/(pmat(GLO.kd,3)*c_max*exp(-pmat(GLO.kd,3)*t_max*0.5));
end

if pmat(GLO.Fs,2) == 1 % only if <Fs> is fitted
    pmat(GLO.Fs,3:4) = [1.05 20]; % this range seems good enough
end

%% BLOCK 4. Finish the parameter matrix
% Value is not important: it is overwritten by <calc_parspace> anyway.
% However, it needs to be within the bounds to avoid errors. 

pmat(:,1) = pmat(:,4); % nothing fancy, just use maximum for the best-value column for now

pmat(:,[1 3 4]) = round(pmat(:,[1 3 4]),4,'significant'); % round initial <pmat> to 4 significant digits


