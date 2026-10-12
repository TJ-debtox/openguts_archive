%% The openGUTS Matlab version
%
% <<http://openguts.info/images/openguts_banner_small.png>>
%
% *Table of contents*
% 

%% About
% Main script for the openGUTS Matlab version. This version started as the
% prototype for the standalone openGUTS software (the latter developed by
% WSC). This script is a sub-set of the general script [openguts.m] to
% illustrate how that script can be modified to focus on a specific task.
% This script will only load the results from a previous calibration and
% perform batch calculations for LPx.
% 
% Input files for calibration results are read from the folder
% [output_sample], and exposure profiles are read from the folder
% [input_profile]. Figures and log files for the output can be found in the
% folder [output_report].
%
% * Author: Tjalling Jager (DEBtox Research)
% * Date: December 2019
% * Web support: <http://www.openguts.info>
% * Funding: Cefic-LRI, see <http://debtox.nl/projects/project_guts2.html>

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

%% BLOCK 1. Initial things
% The helper function [pathdefine] makes sure that the [engine] folder
% (which holds the actual functions that perform the calculations) is in
% the Matlab path. The script [initial_setup] takes care of some general
% setting up. In this block, you can set a name for the analysis (used for
% identifying output files), and set a switch to use a saved sample from a
% previous calibration rather than loading a file and calibrating it.

pathdefine    % places the [engine] folder in Matlab search path
initial_setup % calls a script that sets things up

% Since we are working in this script with saved MAT files only, there is
% no need to provide a name for the study. Furthermore, there is no need to
% specify what to do with the control in calibration or what the
% concentration unit is. All this information is taken from the saved MAT
% file

%% BLOCK 2. Calibration
% A saved calibration is shown (regenerated from the project MAT file). The
% calibration function [calc_optim] can be used for both purposes. The
% Matlab GUI [uigetfile] is used to select files.
%
% Note: the second and third input arguments for [calc_optim] are a switch
% for SD/IT (1/2) and a flag for optimisation (1) or show saved
% optimisation (2). The [plot_main] also has flags for SD/IT (1/2), and a
% second argument, which is a vector of two flags: for the stage of the
% analysis we are in, and for whether to plot model curves or not.

% BLOCK 2.1. Use saved calibration results.
% Note: Optimised parameter values and sample for confidence region
% will be read from a saved MAT file (selected by the Matlab GUI).
% Therefore, there is no need for calibration (as long as the saved
% file exists).
%
% Note: Calling [calc_optim] to show saved optimisation makes [DATAcal]
% from the saved set global, as well as [GLO.basenm], [GLO.concunit]
% and [GLO.fix_hb_cal]. These are all associated with the previous
% calibration, so they should not be changed after calibration has been
% performed.

filename = uigetfile('output_sample/*.mat','Select a saved project file (selecting either SD or IT will use both)','MultiSelect','off'); % use Matlab GUI to select MAT file
if ~iscell(filename) && numel(filename) == 1 && filename(1) == 0 % if cancel is pressed ...
    return % simply stop
end
GLO.basenm = filename(1:end-16); % remove the last part (_parspace_SD.mat) of the filename so that the basename remains
% First, see if there is an output folder for this specific studyname, and otherwise create it.
savestr = ['output_report',filesep,GLO.basenm];
if ~exist(savestr,'dir')
    mkdir(savestr); % if not, create it
end

% Use [calc_optim] to recreate the plot of parameters space and [plot_main]
% to plot the fit (standard plot). Calling [calc_optim] is needed as it
% will define the global settings for the background hazard and the
% concentration unit from the saved MAT file. However, the calls to
% [plot_main] can be commented out.
calc_optim([],1,2); % SD just use saved set, without new optimisation
calc_optim([],2,2); % IT just use saved set, without new optimisation
% plot_main(DATAcal,2,[1 1]) % IT plot of calibration data with the fit (stage 1, plot model on)
% plot_main(DATAcal,1,[1 1]) % SD plot of calibration data with the fit (stage 1, plot model on)

%% BLOCK 6. Batch predictions for LPx
% The function [calc_batch_lpx] can run a number of profiles sequentially,
% without output to screen (apart from a progress bar) until it is
% finished. The function returns the filename for the scenario that leads
% to the lowest (=worst) LP50 value (that could then be used for more
% detailed scrutiny). Plots for survical vs. time will appear in the
% directory [output_report], but not on screen.
% 
% Note that the call to [calc_batch_lpx] includes a correction factor
% [corr_factor] that can be used to correct for different units of the
% exposure profile, if needed. We decided *not* to use that for the
% software. Also note that the batch calculation can optionally be
% performed with CIs on the LPx values. This option was implemented in the
% standalone v.1.1. Calculations with CIs can be slow, so consider using a
% sub-sample from parameter space with [GLO.n_lim].

% BLOCK 6.1. Select files with exposure profiles for batch prediction.
[filename,filepath] = uigetfile('input_profile/*.txt','Select file(s) with exposure profiles for batch predictions','MultiSelect','on'); % use Matlab GUI to select file(s)
if ~iscell(filename) && numel(filename) == 1 && filename(1) == 0 % if cancel is pressed ...
    return % simply stop
end

% Note: the EFSA opinion also uses the propiconazole data set as a case
% study. However, there is some confusion regarding the concentration
% units. The calibration data are reported as nmol/ml (or uM) whereas the
% exposure profiles are in ug/L. These profiles are distributed on the
% openGUTS download page unchanged, to allow the users to reproduce the
% same numbers as in the opinion. Below, we can set a correction factor.
% However, these profiles are examples only, and not meant to be used for
% research or regulatory application.

corr_factor = 1; % to reconstruct the results in the EFSA opinion (and for the software, we don't use this correction factor)
% corr_factor = 1 / 342.2; % correction factor for ug/L (profiles) -> uM (calibration data)

% BLOCK 6.2. Run batch-calculation of all profiles and return worst one.
sdit       = 1; % select death mechanism: 1) SD 2) IT
GLO.n_lim  = [100]; % option to use random sub-sample of given size from the outer rim to speed things up (set empty [] to use full outer rim)
% Using a sub-sample increases speed, but underestimates the width of the
% CIs. How much it underestimates is difficult to say: it depends on the
% sampling, the parameter space, and on the exposure profile. In most
% cases, the difference for a sub-sample of 100 sets will be <5%. However,
% in one case of slow kinetics, and a profile with many peaks, I found a
% difference of 20%. Use with care. Note that this option is *only* used in
% the batch LPx calculations at the moment.

% NOTE: the last flag in the call below is used to signal whether to do
% batch-calculation of LPx without (1) or with (2) CIs!
scen_worst = calc_batch_lpx(filename,filepath,sdit,corr_factor,1); % call [calc_batch_lpx] to do the calculations, and return name of worst profile

return % the return stops the analysis, so we can run this script in parts

% BLOCK 6.3. Detailed analysis for the worst one from the batch calculations.
% Note: this is identical to the procedure in Block 5, but using the
% filename in [scen_worst] rather than selecting files with the Matlab GUI.
[data_pred,error_flag] = load_data(scen_worst,filepath,3,corr_factor); % load and prepare the data set for stage 3 (returned as structured cell array, so [data_all] can include multiple data sets)
if error_flag == 1 % see if there are errors spotted in the input file
    return % simply return (stop and do nothing)
end

LPx = calc_lpx(data_pred,1,sdit); % calculates LPx values (1 with CI, 2 without) (and makes a plot)
plot_main(data_pred,sdit,[3 1],LPx) % plot of data (stage 3, plot model on), pass on the calculated LPx values
% Note: by default, predictions are plotted without CIs. However, setting
% the flag in [plot_main] from [3 1] to [3 2] creates model curves *with*
% CIs (terribly slow for FOCUS profiles). 
