%% Initial setup (script)
% 
% Usage: initial_setup
%
% This script is run once when running the main script. It defines a number
% of global settings. Note that some more settings could be added here
% after the testing phase.
%
% Author     : Tjalling Jager
% Date       : March 2024
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

clear, clear global % clear any remaining workspace and global variables
global GLO SETTINGS_PRED SETTINGS_OPTIM % global structures
global DATAcal %#ok<NUSED> make data set object for calibration global 
% Note: since this is a script, the global definition of <DATAcal> will
% affect the main script, and that is where it is needed. Matlab warns that
% this variable is not used in this script, which is correct, but not
% relevant.

GLO.version = 'v1.2 (20 March 2024)'; % this is the current version information
GLO.save_format = 3; % save plot files as 1) Matlab FIG, 2) JPEG, 3) PDF, 4) PNG, 5) Windows meta file

rng('shuffle') % reseed the random number generator (otherwise we get the exact same result after a restart of Matlab)
diary off % turn of the diary function (if it is accidentaly on)

set(0,'DefaultFigureWindowStyle','normal'); % tell Matlab to make separate figure windows (not tabbed)

% Make sure that versions of Matlab 2017a and newer do not automatically
% update figure legends when new info is plotted in the same plot. This
% cannot be set for older versions (they don't update the legend anyway).
if ~verLessThan('matlab','9.2')
    % -- Code to run in MATLAB R2017a and later here --
    set(groot,'defaultLegendAutoUpdate','off');
end

% Look if there are any waiting bars left open (might happen if a run did
% not terminate properly).
F = findall(0,'type','figure','tag','TMWWaitbar');
delete(F) % close them

% % Note
% % To delete any open waitbar that might stop Matlab from exiting:
% set(0,'ShowHiddenHandles','on')
% delete(get(0,'Children'))

GLO.intermed_plots = 0; % by default, don't show intermediate plots in <calc_parspace>

%% BLOCK 1. Define font properties for plots.
% This is done in this way to make sure that all plots have the same look,
% makes it easy to modify fonts, and allows to adapt the font size for
% small (laptop) screens. This is highly Matlab-specific, so not relevant
% for the software.

% use as: <ft.name,ft.ticks> in <set> statements
% 
% NOTE: with Matlab R2023a, you can only use this format in set statements,
% no longer in xlabel, ylabel, text or title statements! Therefore, with
% version 1.1a of openGUTS, these calls had to be modified to use set
% statements for the text formatting.

ft.name   = {'FontSize','FontName','FontWeight'}; % names of the properties we want to set
ft.label  = {12,'Arial','Normal'}; % font definition for axis labels
ft.ticks  = {10,'Arial','Normal'}; % font definition for ticks (values at axes)
ft.legend = {10,'Arial','Normal'}; % font definition for legends
ft.annot  = {10,'Arial','Normal'}; % font definition for annotations (probably not used anymore)
ft.text   = {10,'Arial','Normal'}; % font definition for text added (as overall title with file name)
ft.title  = {10,'Arial','Bold'};   % font definition for panel titles

% On small screens, Matlab does not allow certain figure sizes, so figure
% will be made smaller in make_fig. This also needs smaller fontsizes to
% work well on small screens.
scrsz = get(groot,'ScreenSize'); % get screensize
if scrsz(4) < 870 % this is probably the limit ... (same limit used in <make_fig>)
    ft.label  = {10,'Arial','Normal'};
    ft.ticks  = {8,'Arial','Normal'};
    ft.legend = {8,'Arial','Normal'};
    ft.annot  = {7,'Arial','Normal'};
    ft.text   = {8,'Arial','Normal'};
    ft.title  = {8,'Arial','Bold'};
end

GLO.ft = ft; % put this structure into the global GLO

%% BLOCK 2. Define some properties of the parameter matrix.
% These are the definitions for the parameter matrix <pmat>. Note that
% <pmat> is a regular matrix. This makes it easier in Matlab to extract or
% manipulate columns, and the optimisation routine fminsearch that I use
% requires a regular vector as input.

% Positions for the parameters in matrix <pmat>, so we can refer to values
% by their names.
GLO.kd = 1;
GLO.mw = 2;
GLO.hb = 3;
GLO.bw = 4;
GLO.Fs = 5;

% Also give them a string for printing/plotting.
GLO.names{GLO.kd} = 'kd';
GLO.names{GLO.mw} = 'mw';
GLO.names{GLO.hb} = 'hb';
GLO.names{GLO.bw} = 'bw';
GLO.names{GLO.Fs} = 'Fs';

% Furthermore, I make the minimum number of time points for the extended
% time vector for SD analysis global. This is used in <calc_model> and in
% <calc_lpx> and it is best if they are always the same. The current
% default will lead to short enough time steps for most cases. Also the
% time vector for plotting is a global. This was a fixed 100 points, which
% works well for calibration but not for very peaky FOCUS profiles (peaks
% in damage were sometimes not shown), so every hour is better.
GLO.t_per_day = 4*24; % number of time points per day for model analysis (default every 15 minutes)
GLO.t_plot    = [100 1*24]; % number of time points for the main plot; minimum and optimised number of points per day (default every hour)

%% BLOCK 3. Define criteria for judging the likelihood ratio.
% Define the chi-square criteria for 1-5 df at 95% confidence (so we don't
% need <chi2inv> from the statistics toolbox). In Matlab this table is
% generated with <(chi2inv(0.95,1:5))'>.

% Critical values for 1-5 d.f. at alpha=0.05 (0.95 CI). Note: df=5 will not
% (yet) be used, as we have 4 parameters at maximum.
GLO.crit_table = [3.8415   
                  5.9915   
                  7.8147   
                  9.4877   
                 11.07    ];

% Next criterion selects the chi-square values to select the sub-sample to
% propagate to model predictions. In principle, this could be only the
% samples that sit on the criterium for 0.95, df=1. However, as we have a
% discrete (and limited) sample, it is best to take a bit more. I suggest
% taking a band around the value for 0.95 df=1. At this moment, it is set
% to the 92-96% area.
GLO.crit_prop = [3.0649 4.2179]; % Note: <[chi2inv(0.92,1) chi2inv(0.96,1)]>
             
%% BLOCK 4. Define settings for LCx and LPx calculations.

SETTINGS_PRED.F_LCx      = [0.50 0.20 0.10]; % which x in LCx to calculate
SETTINGS_PRED.F_LCx_plot = [0.50 0.10]; % which x in LCx to plot (must be sub-set of <F_LCx>)
SETTINGS_PRED.T_LCx      = [1 2 3 4 7 10 14 21 27 28 42 50 100 182]; % time points at which to calculate/plot LCx
% NOTE: added 10, 27 and 182 to accommodate analysis in the framework of
% the 2023 EFSA bee guidance document.

SETTINGS_PRED.F_LPx       = [0.10 0.50]; % which x in LPx to calculate
SETTINGS_PRED.S_range     = [0.05 0.999]; % survival range for making a survival vs. multiplication plot
SETTINGS_PRED.len_S_range = 50; % number of calculated points in survival vs. multiplication plot

%% BLOCK 5. Define settings for optimisation
% These settings are used in <calc_parspace>, and the functions that it
% calls. Note that not all tweakable options in <calc_parspace> are turned
% into a global setting (especially for the main round of sampling). I
% don't think there is much need to change those settings, and it would be
% hard to interpret their value outside of the context of the code.
% However, for the software, we can consider collecting them all as global
% settings here.

% Options for parameter optimisation (used in <calc_parspace>). Note that
% many of these settings are vectors: they change with each subsequent
% round of the analysis to allow the sample to contract.
SETTINGS_OPTIM.crit_add = [ 15  7  3   2   1.5   1  0.5  0.25]; % extra on top of chi2-criterion to select ok values for next round
SETTINGS_OPTIM.n_tr     = [NaN 60 40   30   20   10   10   10]; % number of random parameter tries in round 2 and further, for each ok parameter set
SETTINGS_OPTIM.f_d      = [NaN  1 0.65 0.42 0.27 0.18 0.12 0.08]; % maximum step for random mutations, as factor of initial grid spacing
SETTINGS_OPTIM.n_ok     = 200; % number of values that will at least continue to next round
% Stop criterion: minimum number of values within the total joint CI and
% the inner rim, which is made dependent on number of fitted parameters ...
SETTINGS_OPTIM.n_conf_all = [500  500   % for 1 fitted parameters
                            2500 1500   % for 2 fitted parameters
                            5000 2000   % for 3 fitted parameters
                            7000 3000]; % for 4 fitted parameters

% Number of tries per parameter for the initial grid (used in
% <calc_parspace>). These could be made dependent on the range that is to
% be searched, but there is not too much point. Using many more is
% inefficient, and using less is unnecessary.
SETTINGS_OPTIM.tries_1(GLO.kd) = 14;
SETTINGS_OPTIM.tries_1(GLO.mw) = 12;
SETTINGS_OPTIM.tries_1(GLO.hb) = 8;
SETTINGS_OPTIM.tries_1(GLO.bw) = 12;
SETTINGS_OPTIM.tries_1(GLO.Fs) = 12;

% Settings for deciding on a restart for slow kinetics (used in <calc_parspace>).
SETTINGS_OPTIM.slowkin_corr = 0.70; % minimum correlation coefficient between <kd> and <mw>
SETTINGS_OPTIM.slowkin_pars = 0.05; % closeness of <kd> and <mw> to their lower bound (as fraction of total range)
% Settings for defining new upper bounds after restart (used in <calc_optim>)
SETTINGS_OPTIM.slowkin_f_mw = 3;  % factor by which to multiply current-highest <mw> to get new upper bound
SETTINGS_OPTIM.slowkin_f_kd = 10; % factor by which to multiply current-highest <kd> to get new upper bound

% Additional settings for the extra sampling (when profiling shows that
% sampling was poor in parts of the profile).
SETTINGS_OPTIM.crit_add_extra = 3; % continue with sets that are within this value from the MLL (this focusses on parameters within, or close to, the inner rim)
SETTINGS_OPTIM.f_d_extra      = 1; % initial maximum jump size as fraction of the grid spacing
SETTINGS_OPTIM.d_extra        = 0.1; % grid spacing as fraction of the parameter's total range (from joint 95% CI)
SETTINGS_OPTIM.gap_extra      = 0.25; % the gap distance between profile and sample that triggers resampling
SETTINGS_OPTIM.real_better    = 0.05; % difference between new optimum and old optimum (MLL) that leads to printing on screen
