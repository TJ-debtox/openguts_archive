function disp_pmat(pmat,sdit,pmat_print)

% Usage: disp_pmat(pmat,sdit,pmat_print)
%
% This function prints <pmat> and <pmat_print> on screen in a formatted
% manner. This is the only place where the 5 model parameters and their CIs
% are printed to screen.
% 
% Inputs
% <pmat>        parameter matrix
% <sdit>        switch for SD (1) or IT (2)
% <pmat_print>  matrix with information on each parameter (best value and CI info)
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

global GLO

%% BLOCK 1. Initial things.
par_names = GLO.names; % names of the model parameters

% Define strings for the printing units of the parameters.
par_unit{GLO.kd} = '1/d';
par_unit{GLO.mw} = GLO.concunit;
par_unit{GLO.bw} = ['1/(',GLO.concunit,' d)'];
par_unit{GLO.hb} = '1/d';
par_unit{GLO.Fs} = '[-]';

n_par = size(pmat,1); % number of parameters in <pmat> in total

%% Block 2. Print parameter values/ranges/intervals.

if isempty(pmat_print) % then we only print bounds (before an optimisation run)
    
    disp('=====================================================================')
    for i_p = 1:n_par % run through parameters
        if pmat(i_p,5) == 0 % create a string for log or normal scale fitting
            logsc = 'log';
        else
            logsc = 'norm';
        end
        fprintf('%-5s bounds: %10.4g - %10.4g %-8s fit: %1.0f (%s)\n',par_names{i_p},pmat(i_p,3),pmat(i_p,4),par_unit{i_p},pmat(i_p,2),logsc)
    end
    disp('=====================================================================')
    
else % then we print the results of a calibration run and more info needs to be printed
    
    % Replace flagged values in <pmat_print> by an ASCII code for display.
    prt_tst = 32*ones(n_par,2); % make matrix with spaces (two columns)
    prt_tst(pmat_print(:,6:7)==1) = 42; % replace 'hitting bounds' flag by asterisk

    disp('==========================================================================')
    for i_p = 1:n_par % run through parameters
        if pmat(i_p,5) == 0 % create a string for log or normal scale fitting
            logsc = 'log';
        else
            logsc = 'norm';
        end
        % Print the line to screen (broken up in three command lines).
        if pmat(i_p,2) == 1 % for fitted parameters
            fprintf('%-5s best: %#10.4g (%#10.4g%c -',par_names{i_p},pmat_print(i_p,1),pmat_print(i_p,2),char(prt_tst(i_p,1)))
            fprintf('%#10.4g%c) %-10s',pmat_print(i_p,3),char(prt_tst(i_p,2)),par_unit{i_p})
        else % not fitted, so no need to print trailing zeros
            fprintf('%-5s best: %10.4g (%10.4g%c -',par_names{i_p},pmat_print(i_p,1),pmat_print(i_p,2),char(prt_tst(i_p,1)))
            fprintf('%10.4g%c) %-10s',pmat_print(i_p,3),char(prt_tst(i_p,2)),par_unit{i_p})
        end
        fprintf('fit: %1.0f (%s)\n',pmat(i_p,2),logsc)
    end
    disp('==========================================================================')
    if any(prt_tst(:) == 42) % we are printing at least one asterisk
        disp('* edge of 95% parameter CI has run into a boundary')
        disp('  (this may also have affected CIs of other parameters)')
    end
    
end

if sdit == 1 % also print the special case we have just printed
    fprintf('Special case: stochastic death (SD) \n');
else
    fprintf('Special case: individual tolerance (IT) \n');
end
if GLO.fix_hb_cal % did we fix <hb> to the controls?
    disp('Background hazard rate fixed to a value fitted on controls')
end

%% Block 3. After calibration, we need to print some additional information

if ~isempty(pmat_print) % then we print the results of a calibration run and more info needs to be printed
    disp(' ')
    disp('Extra information from the fit:')
    disp('==========================================================================')
    
    for i_p = 1:n_par % run through parameters
        if pmat_print(i_p,8) == 1 % is there a flag for a broken CI?
            disp(['Confidence interval for parameter ',par_names{i_p},' is a broken set.'])
        end
    end
    
    if sdit == 2 % for IT, also show beta of the log-logistic distribution
        beta    = log(39) / log(pmat_print(GLO.Fs,1)); % translate <Fs> to beta
        beta_lo = log(39) / log(pmat_print(GLO.Fs,3)); % also translate CI (note that low <Fs> means high beta)
        beta_hi = log(39) / log(pmat_print(GLO.Fs,2)); 
        fprintf('Beta of the threshold distribution: %#1.3g (%#1.3g - %#1.3g) (-) \n',beta,beta_lo,beta_hi)
    end
    
    % Calculate and print the DRT95 from <kd>.
    DRT95   = log(20)/pmat_print(GLO.kd,1); % calculate DRT95 from <kd>
    DRT95lo = log(20)/pmat_print(GLO.kd,3); % and its CI
    DRT95hi = log(20)/pmat_print(GLO.kd,2);
    fprintf('Depuration/repair time (DRT95)    : %#1.3g (%#1.3g - %#1.3g) days \n',DRT95,DRT95lo,DRT95hi)
    
    % Due to popular demand: print the parameter search ranges again ...
    disp('Parameter ranges used for the optimisation:')
    for i_p = 1:n_par % run through parameters
        fprintf('%-5s range: %10.4g - %10.4g \n',par_names{i_p},pmat(i_p,3),pmat(i_p,4))
    end
    disp('==========================================================================')
    
end

