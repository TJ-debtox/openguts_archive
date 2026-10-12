function [pmat_print,mll] = calc_optim(pmat,sdit,flag)

% Usage: [pmat_print,mll] = calc_optim(pmat,sdit,flag)
%
% Optimisation and finding a likelihood-based joint confidence region
% without the need for starting values. However, the parameter space to be
% searched must limited for the algorithm to be effective. Parameter space
% to be searched is defined as the min-max bounds of the parameters in the
% matrix <pmat>.
% 
% In principle, no output would be needed. However, for some forms of
% testing, output is handy. This function is a shell around the calibration
% algorithm in <calc_parspace>. It prepares for the optimisation, takes
% care of a restart (when slow kinetics is indicated during the
% optimisation), and it plots parameter space for saved sets if needed.
% 
% Note that the function <disp_pmat> is called to do the formatted printing
% (on screen) of <pmat> and <pmat_print>. The <pmat> is displayed before
% optimisation (to show the search ranges), the <pmat_print> after (with
% the optimised values and CIs).
%
% Inputs
% <pmat>        parameter matrix (5x5 with search ranges, fit/fix, etc.)
% <sdit>        switch for SD (1) or IT (2)
% <flag>        1) regular optimisation, 2) use saved set
% 
% Outputs
% <pmat_print>  matrix with parameter info (best fit, CIs, etc.)
% <mll>         minus log-likelihood for the fit
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

global GLO SETTINGS_OPTIM % make structure GLO and settings global
diary (['output_report',filesep,GLO.basenm,filesep,GLO.basenm,'.out']) % collect any output in the diary *.out

%% BLOCK 1. Use the saved set, and don't make a new calculation.

if flag == 2 
    
    [pmat,~,pmat_print,mll] = load_sample(sdit,2); % load info from the saved optimisation
    
    % And display it to screen and diary.
    n_fit = sum(pmat(:,2)); % number of fitted parameters
    disp('=================================================================================')
    disp('Displaying and plotting optimisation results from a saved set')
    disp('=================================================================================')
    fprintf('   Base name    : %s \n',GLO.basenm) % base name of the analysis
    disp(['   Minus log-likelihood: ',num2str(mll),'  AIC: ',num2str(2*n_fit+2*mll)])
    disp('Best fit parameter set and 95% bounds from the saved set:')
    
    disp_pmat(pmat,sdit,pmat_print); % call dedicated function for printing the <pmat> and <pmat_print>
    disp(' ')
    diary off
    return % nothing further to do, so return to the main script

end

%% BLOCK 2. Call <calc_parspace> to do an optimisation.

if any(pmat(:,3) == 0 & pmat(:,5) == 0)
    error('For log parameters, lower bound cannot be zero. Check the parameter structure.')
    % This check is best performed when people are modifying <pmat> in the
    % 'expert setting'.
end
if any(pmat(:,4)-pmat(:,3)<=0 & pmat(:,2)==1)
    error('For at least one parameter, fitting is requested but the search range is of zero or negative width.')
    % This check may seem trivial; it is only needed when someone modifies
    % the <pmat> in the 'expert setting'.
end

disp(' ')
disp('Settings for parameter search ranges:')
disp_pmat(pmat,sdit,[]) % call <disp_pmat> for a formatted printing of the parameter search ranges

tic % turn on the stopwatch to time the optimisations
% Call the optimisation function <calc_parspace> to do the optimisation.
[minmax,mll,stats,pmat_print] = calc_parspace(pmat,sdit);

% The <calc_parspace> will break the analysis if it runs into slow kinetics
% and when <mw> is estimated on normal scale. In that case, it needs to
% restart with <mw> on log-scale, and additionally, the parameter ranges
% can probably be decreased a bit.
if minmax(1) ~= -1 % then we stopped before the end ... this should only happen once (so no 'while' needed)
    
    disp(' ')
    % We suspect slow kinetics (currently the only reason for a restart).
    disp('Slow kinetics indicated: parameter-space explorer is restarting with these settings:')
    % <calc_parspace> returned the min-max of (relevant part of) parameter
    % cloud on normal scale.
    pmat(GLO.mw,5) = 0; % for <mw>, put fit flag on log scale
    % Note: if someone *really* wants to keep <mw> on normal scale, there
    % is currently no simple way to avoid triggering a restart and a
    % log-scale for <mw>.
    
    % Note: for <kd>, we don't need to put it on log scale, as it already
    % is (see <calc_parspace>: if <kd> is on normal scale, there will be no
    % check on slow kinetics anyway, so no restart).
    
    % Adapt parameter bounds in <pmat> where possibel, using the values in <minmax>.
    pmat(GLO.mw,4) = min(pmat(GLO.mw,4) , minmax(GLO.mw,2)*SETTINGS_OPTIM.slowkin_f_mw); % if possible, adapt upper boundary of <mw>
    pmat(GLO.kd,4) = min(pmat(GLO.kd,4) , minmax(GLO.kd,2)*SETTINGS_OPTIM.slowkin_f_kd); % also modify upper bound <kd> a bit if possible ...
    if sdit == 1 % only for SD
        pmat(GLO.bw,3) = max(pmat(GLO.bw,3) , minmax(GLO.bw,1)/3); % also modify lower boundary <bw> a bit if possible ...
    else
        % also modify bounds on <Fs>? Probably not needed
    end
    % The modifications of the upper boundaries for <mw> and <kd>, and
    % lower boundary for <bw>, are bit tricky but seem to work ... if they
    % do not work, that will be visible in the parameter-space plot!
    
    % Display new to-be-fitted parameter ranges on screen.
    disp_pmat(pmat,sdit,[]); % call dedicated function for printing the <pmat>
    
    % And restart the optimisation with <calc_parspace>.
    [~,mll,stats,pmat_print] = calc_parspace(pmat,sdit);
end

%% BLOCK 3. Final printing on screen.

% The stats vector contains some interesting info from the analysis.
ind_final = stats(1); % index to the last element in the sample in the joint 95% CI
ind_inner = stats(2); % index to the last element in the sample in the inner rim
ind_prop  = stats(3); % number of elements in the sample that are used for propagating errors

disp(' ')
disp('=================================================================================')
disp('Results of the parameter estimation')
disp('=================================================================================')
  disp(['   openGUTS Matlab: ',GLO.version])
fprintf('   Base name      : %s \n',GLO.basenm) % base name of the analysis
fprintf('   Analysis date  : %s \n',datetime('now'));
% NOTE: function <datetime> was introduced in this functionality in Matlab
% R2016b, so earlier versions throw an error here.
disp('   Following data sets are loaded for calibration:')
for i_d = 1:length(GLO.filename) % run through data sets loaded
    disp(['     set: ',num2str(i_d),', file: ',GLO.filename{i_d}]) % and display their filename
end
n_fit = sum(pmat(:,2)); % number of fitted parameters
disp(['   Sample: ',num2str(ind_final),' sets in joint CI and ',num2str(ind_inner),' in inner CI.'])
disp(['   Propagation set: ',num2str(ind_prop),' sets will be used for error propagation.'])
% Display minloglik and Akaike Information Criterion.
fprintf('   Minus log-likelihood has reached %#1.2f (AIC=%1.2f). \n',mll,2*n_fit+2*mll)

disp('Best estimates and 95% CIs on single parameters')
disp_pmat(pmat,sdit,pmat_print); % call dedicated function for printing the <pmat>

disp(['Time required for optimisation: ' secs2hms(toc)]) % for testing: print time needed
disp(' ')
diary off  % close the diary that collected the output to screen
