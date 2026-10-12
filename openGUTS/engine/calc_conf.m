function [model_p,model_dam] = calc_conf(t,ev_ic,sdit,hb_fix,flag_conf)

% Usage: [model_p,model_dam] = calc_conf(t,ev_ic,sdit,hb_fix,flag_conf)
%
% This function creates confidence intervals on model curves from the
% parameter-space sample that was saved by <calc_parspace>. This script
% calculates the best model curve and the bounds from the saved set. There
% is the option to fix the background hazard to a specific value. This is
% helpful when validation uses a fixed <hb> and calibration did not. Also,
% this is helpful when plotting the results for an LPx analysis: to plot
% survival vs. time for a specific multiplication factor. This implies that
% we have to set background hazard to zero. 
%
% Inputs
% <t>          time vector for calculating model output
% <ev_ic>      an exposure-events matrix for a single treatment
% <sdit>       switch for SD (1) or IT (2)
% <hb_fix>     fixed value for the background hazard rate (-1 is not fixed)
% <flag_conf>  1) calculate CIs, 0) do not calculate CIs
% 
% Outputs; two matrices with three columns each (best, minimum, maximum)
% <model_p>    the best curve for survival probability, and CI
% <model_dam>  the best curve for damage, and CI
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

global GLO % make structure <GLO> global

%% BLOCK 1. Initial things

% Load file as saved by <calc_parspace> and prepare sample.
[pmat,sample] = load_sample(sdit,1); % load random sample and parameter matrix from saved MAT file

% Apply fixed background hazard rate when/where needed.
if hb_fix ~= -1 % then we are requested to include a fixed hazard rate (even when <hb> was fitted in calibration)
    pmat(GLO.hb,[1 5]) = [hb_fix 1]; % set background hazard to <hb_fix> and make sure it is normal scale!
    sample(:,GLO.hb)   = hb_fix;     % also make it <hb_fix> in all sets of the sample!
end

%% BLOCK 2. Calculate model curves for each set in the sample.

drawnow % empty plot buffer if there's something in it 
% Note: this should not be needed, but it makes sure that plotting will not
% wait for the calculations here to finish.

% BLOCK 2.1. Start by running the best fitting set (in <pmat>) to produce a
% best-fitting curve.

% Initialise the matrices that will collect the model output (with CIs).
model_p   = nan(length(t),3); % initialise output matrix for survival output with NaNs
model_dam = nan(length(t),3); % initialise output matrix for damage output with NaNs

x_out          = calc_model(t,pmat,ev_ic,sdit); % use <calc_model> to provide the output for this scenario (2-column matrix)
model_p(:,1)   = x_out(:,1); % collect survival probability into <model_p> (first column in matrix)
model_dam(:,1) = x_out(:,2); % collect damage into <model_dam> (first column in matrix)
    
% BLOCK 2.2. Prepare for running through the sample.

n_sets  = size(sample,1); % number of sets from the joint confidence region in the sample
ind_log = pmat(:,5)==0;   % vector with logical indices to parameters on log scale

if flag_conf == 0 % then there is no need to return CIs, just the best model line
    % Note: this returns the CIs as NaN, so they are ignored when plotting
    return
else
    sample = sample(randperm(size(sample,1)),:); % randomly shuffle sets from outer hull
    % Note: this is for the Matlab version only, and only for the
    % calculation of model CIs for LPx calculations, to get more
    % appropriate CIs when the user presses the <cancel> button in the wait
    % bar.
end

% BLOCK 2.3. Run through the sample to calculate CIs on curves.

% Use the fancy progress bar from Matlab
f = [];
if size(ev_ic,1) > 100 % only for long calculation vectors (otherwise bar takes too much time)
    f = waitbar(0,sprintf('Fraction of %1.0f sets calculated',n_sets),'Name','Calculating CIs model curves','CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
end

for i_n = 1:n_sets % run through all sets in the sample
    
    % Progress bar from Matlab: check if cancel is pressed
    if ~isempty(f)
        if getappdata(f,'canceling')
            break % break the <for> loop
        end
        waitbar(i_n/n_sets,f) % update waitbar
    end
    
    pmat(:,1) = sample(i_n,:); % replace values in first column of <pmat> with the <i_n>th set from the sample
    % Put parameters that need to be on log scale back on normal scale (as
    % <calc_model> requires normal scale).
    if sum(ind_log)>0
        pmat(ind_log,1) = 10.^(pmat(ind_log,1));
    end
    % Note: <pmat> is still on normal scale here, but the sample contains
    % the value on a log scale, if a parameter is fitted on log scale.
    % <calc_model> needs normal scale <pmat> (also for the fixed parameters).

    x_out = calc_model(t,pmat,ev_ic,sdit); % use <calc_model> to provide the model output for this scenario
    
    % I don't save all individual curves, but just update min-max for every
    % set. Note that this takes the minimum across two columns from a matrix,
    % for each row.
    model_p(:,2)   = min([model_p(:,2)  , x_out(:,1)],[],2); % smallest value of current lowest and new model line
    model_p(:,3)   = max([model_p(:,3)  , x_out(:,1)],[],2); % largest value of current lowest and new model line
    model_dam(:,2) = min([model_dam(:,2), x_out(:,2)],[],2); % smallest value of current lowest and new model line
    model_dam(:,3) = max([model_dam(:,3), x_out(:,2)],[],2); % largest value of current lowest and new model line
    
end

if ~isempty(f)
    delete(f) % delete progress bar again
    if i_n < n_sets % if cancel was pressed, display that on screen
        fprintf('Calculation of CIs for plotting stopped by user after %1.0f of %1.0f sets \n',i_n,n_sets)
        disp('(note that the sample was randomised before CI calculation)')
    end
end
