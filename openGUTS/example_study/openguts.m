%% The openGUTS Matlab version
%
% <<http://openguts.info/images/openguts_banner_small.png>>
%
% *Table of contents*
% 

%% About
% Main script for the openGUTS Matlab version. This version started as the
% prototype for the standalone openGUTS software (the latter developed by
% WSC). This programme follows the workflow as laid down in the 2018 EFSA
% Scientific Opinion on TKTD modelling:
% 
% * Calibration: load one or more data sets from file using the Matlab
% input GUI. The input data set will be shown in a plot, and calibration of
% GUTS-RED-SD and GUTS-RED-IT will commence automatically. Fitting results
% will be shown on screen.
% * Validation: load one or more data sets from file using the Matlab
% input GUI. The results from the calibration will be used to predict the
% survival probability for the exposure scenarios in the validation data
% sets. Results will be shown on screen.
% * Prediction: load one or more exposure profiles from file using the
% Matlab input GUI. The results from the calibration will be used to
% predict the survival probability for the exposure scenarios in the
% exposure profiles. Results will be shown on screen.
% 
% Input files for survival data and exposure profiles are read from the
% respective folders [input_data] and [input_profile]. Information from the
% calibration is saved into the folder [output_sample]. Figures and log
% files for the output can be found in the folder [output_report].
%
% * Author: Tjalling Jager (DEBtox Research)
% * Date: May 2020
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

% rng(2.2) % seed the random number generator with a specific value (this overrides the shuffle in the [initial_setup])
% % Fixing the seed is helpful for testing: it ensures that every run of this file will be exactly the same.

% Provide a name for the study to be used for the plots and saved files.
GLO.basenm       = 'test'; % name used for identifying the analysis
use_saved_sample = 0; % when set to 1, use saved calibration file rather than make a new calibration
GLO.fix_hb_cal   = 1; % when set to 1, fit [hb] to the control data for all calibration data sets, and keep fixed
% Note: for the validation data set, there is a separate switch for fixing [hb] in BLOCK 4.
GLO.intermed_plots = 0; % when set to 1, makes intermediate plots of parameter space to monitor progress

%% BLOCK 2. Calibration
% Based on the setting of [use_saved_sample], a saved calibration is shown
% (regenerated from the project MAT file) or a new calibarion is performed.
% The calibration function [calc_optim] can be used for both purposes. The
% Matlab GUI [uigetfile] is used to select files.
%
% Note: the second and third input arguments for [calc_optim] are a switch
% for SD/IT (1/2) and a flag for optimisation (1) or show saved
% optimisation (2). The [plot_main] also has flags for SD/IT (1/2), and a
% second argument, which is a vector of two flags: for the stage of the
% analysis we are in, and for whether to plot model curves or not.

if use_saved_sample == 1 
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

    % Use [calc_optim] to recreate the plot of parameters space and [plot_main] to plot the fit (standard plot).
    calc_optim([],1,2); % SD just use saved set, without new optimisation
    calc_optim([],2,2); % IT just use saved set, without new optimisation
    plot_main(DATAcal,2,[1 1]) % IT plot of calibration data with the fit (stage 1, plot model on)
    plot_main(DATAcal,1,[1 1]) % SD plot of calibration data with the fit (stage 1, plot model on)  
    
else
    % BLOCK 2.2. Load and prepare calibration data.
    % The data set will be loaded and translated to the internal
    % representation, which is returned in [DATAcal]. [DATAcal] is made
    % global (in [initial_setup]) as it is called many times during
    % optimisation.
    
    [filename,filepath] = uigetfile('input_data/*.txt','Select file(s) for calibration','MultiSelect','on'); % use Matlab GUI to select file(s)
    if ~iscell(filename) && numel(filename) == 1 && filename(1) == 0 % if cancel is pressed ...
        return % simply stop
    end
    
    [DATAcal,error_flag] = load_data(filename,filepath,1); % load and prepare the data set for stage 1 (returned as structured cell array, so [data_all] can include multiple data sets)
    if error_flag == 1 % see if there are errors spotted in the input file
        return % simply return (stop and do nothing)
    end
    plot_main(DATAcal,[],[1 0]) % inspection plot of data (stage 1, plot model off)
    % Note: the 1 in the calls above signal that we are loading, preparing
    % and plotting data for stage 1: calibrations.
    
    % BLOCK 2.3. Ranges for the model parameters. Note: the parameter
    % structure [pmat] will be automatically generated by [startgrid]
    % below. Experts can make changes to the [pmat] after that. 
    %
    % WARNING: some settings of [pmat] will lead to poor performance or
    % even complete failure of the calibration algorithm. Therefore, when
    % messing with this matrix, expertise and care are essential to get
    % reliable results.
    
    % Make separate parameter matrices for SD and IT.
    pmat_SD = startgrid(1); % create reasonable SD parameter ranges from the data set
    pmat_IT = startgrid(2); % create reasonable IT parameter ranges from the data set
    % Below this line, modifications to the parameter matrices can be inserted.

    % BLOCK 2.4. Perform calibration.
    % Note that the starting values in [pmat] are only used for fixed
    % parameters in the analysis; for fitted parameter, only the ranges and
    % other settings are used. [calc_optim] finds the global optimum,
    % confidence intervals on parameters, and the joint confidence region
    % (to be used for predictions). The sample is saved in a MAT file with
    % an extension indicating the analysis (SD or IT), which is
    % automatically used by the post-analyses.
    
    % Use [calc_optim] to perform the calibration and [plot_main] to plot the fit (standard plot).
    calc_optim(pmat_SD,1,1); % SD, new optimisation
    calc_optim(pmat_IT,2,1); % IT, new optimisation
    plot_main(DATAcal,1,[1 1]) % SD plot of calibration data with the fit (stage 1, plot model on)
    plot_main(DATAcal,2,[1 1]) % IT plot of calibration data with the fit (stage 1, plot model on)
    
end

% return % the return stops the analysis, so we can run this script in parts

%% BLOCK 3. Predictions: LCx versus time
% LCx,t values are calculated for several effect percentages and several
% preset time points (defined in [initial_setup]). Optimised parameter
% values and sample to create CIs will be read from the saved MAT file. The
% final timepoint in the call to [calc_lcx] is for plotting only; all
% timepoints defined in [initial_setup] are given on screen.

calc_lcx(16,1); % SD, calculate LCx values with CI (final timepoint t=16)
calc_lcx(16,2); % IT, calculate LCx values with CI (final timepoint t=16)

return % the return stops the analysis, so we can run this script in parts

%% BLOCK 4. Validation: comparison to new data
% The function [plot_main] can directly be used to make predictions (with
% CI) for new scenarios (not used for calibration) and compare them with
% data. As input, it needs a new data set, in the same way as for the
% calibration data set. This also implies that it is possible to use
% multiple data sets (although we allow only 1 for the standalone
% software). Optimised parameter values and confidence region will be read
% from a saved MAT file.
% 
% The data set is treated in the same manner as the calibration data. It is
% translated to the internal representation, which is returned in
% [data_val]. This definition is NOT made global; it is simply passed on to
% the plotting routine.

GLO.fix_hb_val = 1; % when set to 1, fit [hb] to the control data for all validation data sets, and keep fixed

% BLOCK 4.1. Load and prepare validation data.
[filename,filepath] = uigetfile('input_data/*.txt','Select file(s) for validation','MultiSelect','on'); % use Matlab GUI to select file(s)
if ~iscell(filename) && numel(filename) == 1 && filename(1) == 0 % if cancel is pressed ...
    return % simply stop
end

[data_val,error_flag] = load_data(filename,filepath,2); % load and prepare the data set for stage 2 (returned as structured cell array, so [data_all] can include multiple data sets)
if error_flag == 1 % see if there are errors spotted in the input file
    return % simply return (stop and do nothing)
end
plot_main(data_val,1,[2 0]) % inspection plot of data (stage 2, plot model off)
% Note: the 2 in the calls above signal that we are loading, preparing and
% plotting data for stage 2: validations.

% BLOCK 4.2. Call the plotting routine for validation plots (with model output).
plot_main(data_val,1,[2 1]) % SD validation: compare to model predictions (stage 2, plot model on)
plot_main(data_val,2,[2 1]) % IT validation: compare to model predictions (stage 2, plot model on)

return % the return stops the analysis, so we can run this script in parts

%% BLOCK 5. Predictions: LPx
% LPx is calculated for various effect levels (defined in [initial_setup]).
% This is the factor by which a given exposure scenario needs to be
% multiplied to result in x% mortality at the end of the profile, relative
% to the control. Optimised parameter values and confidence region will be
% read from the saved MAT file.
%
% Exposure scenarios need to be entered as text files. However, only one
% treatment is allowed per data set. Multiple data sets can be used
% (although we allow only 1 for the standalone software). The scenario is
% thus a matrix with two columns: time (in days) and exposure concentration
% (in the same units as the data used for calibration!). This will be 
% linearly interpolated.
%
% Note: the second entry in the call to [calc_lpx] is a flag that is
% used for faster calculations (no CI) and batch processing. Note that flag
% setting 2 will probably not be used in the standalone software, but it is 
% handy for testing (as CIs on LPx can take a very long calculation time).
% flag   1) regular, 2) no CI (very fast), 3) batch processing (no output
%        to screen, no diary, no CIs)

% BLOCK 5.1. Load and prepare prediction data.
[filename,filepath] = uigetfile('input_profile/*.txt','Select file(s) with exposure profiles for standard predictions','MultiSelect','on'); % use Matlab GUI to select file(s)
if ~iscell(filename) && numel(filename) == 1 && filename(1) == 0 % if cancel is pressed ...
    return % simply stop
end

[data_pred,error_flag] = load_data(filename,filepath,3); % load and prepare the data set for stage 3 (returned as structured cell array, so [data_all] can include multiple data sets)
if error_flag == 1 % see if there are errors spotted in the input file
    return % simply return (stop and do nothing)
end
plot_main(data_pred,[],[3 0],[]) % inspection plot of data (stage 3, plot model off)
% Note: the 3 in the calls above signal that we are loading, preparing and
% plotting data for stage 3: predictions.

% BLOCK 5.2. Calculate LPx and plot.
sdit = 1; % select death mechanism: 1) SD 2) IT
LPx  = calc_lpx(data_pred,1,sdit); % calculates LPx values (1 with CI, 2 without) and makes a plot
plot_main(data_pred,sdit,[3 1],[ones(size(LPx,1),1) LPx]) % plot of data (stage 3, plot model on), pass on the calculated LPx values
% Note: by default, predictions are plotted without CIs. However, setting
% the flag in [plot_main] from [3 1] to [3 2] creates model curves *with*
% CIs (terribly slow for FOCUS profiles). Note: the vector
% [ones(size(LPx,1),1) LPx] in [calc_lpx] specifies which multiplication
% factors to use for plotting. I add the [ones] here to also plot the
% survival due to the unmodified exposure profile (as is also done in the
% standalone software).

return % the return stops the analysis, so we can run this script in parts

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

return % the return stops the analysis, so we can run this script in parts

%% BLOCK 7. Demonstration of predictions for test design
% The ability to make predictions can also be used for efficient test
% design. Here, an example is given. The user should prepare one or more
% scenarios of choice, e.g., with various pulses. With MF, the user can
% specify multiplication factors to plot, and [plot_main] will plot the
% associated damage and survival over time (here with CIs).
% 
% For testing, after calibration with propiconazole data for constant
% exposure, try the two test profiles in the directory
% [pulses_test_design]. These also work well with the dieldrin case study.

% BLOCK 7.1. Select file(s) with exposure profiles for prediction.
[filename,filepath] = uigetfile('input_profile/*.txt','Select file(s) with exposure profiles for test-design predictions','MultiSelect','on'); % use Matlab GUI to select file(s)
if ~iscell(filename) && numel(filename) == 1 && filename(1) == 0 % if cancel is pressed ...
    return % simply stop
end

[data_pred,error_flag] = load_data(filename,filepath,3); % load and prepare the data set for stage 3 (returned as structured cell array, so [data_all] can include multiple data sets)
if error_flag == 1 % see if there are errors spotted in the input file
    return % simply return (stop and do nothing)
end

% BLOCK 7.2. Plot predictions.
sdit = 1; % select death mechanism: 1) SD 2) IT
MF   = [1 2 2.5 5]; % these are the multiplication factors to plot
MF   = ones(data_pred.nr_sets,1) * MF; % make sure there is the same MF for all profiles that are entered

plot_main(data_pred,sdit,[3 2],MF) % plot of data (stage 3, plot model on), with the different MF values
% Note: the flag in [plot_main] of [3 2] creates model curves *with* CIs.
