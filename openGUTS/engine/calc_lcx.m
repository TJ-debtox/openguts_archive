function [LCx,LClo,LChi] = calc_lcx(Tend,sdit)

% Usage: [LCx,LClo,LChi] = calc_lcx(Tend,sdit)
%
% Function to calculate LCx,t values, with confidence intervals. For IT,
% this function uses the complete analytical solution for the LCx,t
% (provided in the function <calc_surv_it>). For SD, it uses the analytical
% solution for survival (in <calc_surv_sd_c>), but still needs <fzero> to
% find the LCx. If we ever need to implement the full GUTS model (with
% separate TK and damage) some shortcuts need to be reconsidered!
% 
% The output of this function is used for testing purposes only: matrices
% for LCx,t and the corresponding lower and higher end of the CI. No
% outputs would be needed as this function prints results to screen and to
% the diary.
%
% Inputs
% <Tend>   the maximum time point for the LCx-time plot
% <sdit>   switch for SD (1) or IT (2)
% 
% Outputs
% <LCx>    matrix with LCx,t values (time in rows, increasing, x in columns, decreasing)
% <LClo>   lower edge of the 95% CI; same format as <LCx>
% <LChi>   upper edge of the 95% CI; same format as <LCx>
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

global GLO SETTINGS_PRED % make several things global

%% BLOCK 1. Initial things.

% Create text string for plotting titles and saving plots.
if sdit == 1
    savestr = [GLO.basenm,'_LCxt_SD'];
else
    savestr = [GLO.basenm,'_LCxt_IT'];
end

% Extract some settings from the globals defined in <initial_setup>.
F_LCx      = SETTINGS_PRED.F_LCx; % effect fractions for LCx (x/100)
T_LCx      = SETTINGS_PRED.T_LCx; % time points at which to calculate LCx for output to screen
F_LCx_plot = SETTINGS_PRED.F_LCx_plot; % effect fractions for LCx for plotting
ft         = GLO.ft; % global for font specs in plots

% These checks are no longer needed as we now have fixed time points and
% fixed effect fractions. However, to be on the safe side ... Note: time
% needs to increase and effect fraction need to decrease to reduce the
% range that <fzero> has to search for SD (as LCx,t has to decrease with
% increasing time and with dereasing x).
T_LCx = sort(T_LCx); % make sure that time vector increases
F_LCx = sort(F_LCx,'descend'); % make sure that effect fraction decreases
if T_LCx(1) == 0 % time zero won't work as there is no effect there
    T_LCx(1) = []; % just remove the first entry
end

diary (['output_report',filesep,GLO.basenm,filesep,GLO.basenm,'.out']) % collect output in the diary *.out

%% BLOCK 2. Load calibration file as saved by <calc_parspace>.

[pmat,sample] = load_sample(sdit,1); % load sample (matrix) and parameter matrix from saved MAT file

%% BLOCK 3. Calculate LCx,t
% Several model functions will be used to (help) find the LCx,t. For IT,
% LCx,t can be calculated analytically, but for SD, use of <fzero> (Matlab
% built-in) or <zero> (as used in the standalone) is needed. In the latter
% case, the analytical solution for survival is used, which is possible as
% the LCx is per definition for constant exposure. This implies that
% <calc_model> is NOT used for these calculations. Background hazard is
% ignored as LCx,t is relative to the control.

% Initialise matrices to collect the LCx and CIs over time for the different effect levels.
LCx  = zeros(length(T_LCx),length(F_LCx)); % initialise matrix to collect LCx,t
LClo = nan(length(T_LCx),length(F_LCx));   % initialise with NaNs is always safe
LChi = nan(length(T_LCx),length(F_LCx));   % initialise with NaNs is always safe

if sdit == 2 % for IT, we can use a direct calculation in <calc_surv_it>
    
    for i_t = 1:length(T_LCx) % run through time points for the LCx
        for i_e = 1:length(F_LCx) % run through effect levels requested
            [~,LCx(i_t,i_e)] = calc_surv_it(1,pmat,T_LCx(i_t),F_LCx(i_e),[]); % calculate LCx
        end
    end
    
else % for SD, some more work is needed as we need to use root finding
    
    % NEW: make sure that LCinit is always a range
    LCinit = [pmat(GLO.mw,1)/10 pmat(GLO.mw,1)]; % as search range for LCx,t start with <mw>/10 and <mw>
    % Using <mw> is fine for lower end of range but not for upper, so find
    % a better value for the upper range by trial and error. I do this
    % by shifting this initial range up until it covers the sign change
    % of <crit>.
    crit = 1; % first one will be positive
    while crit > 0
        LCinit = LCinit * 10; % shift lower and upper range by factor of 10
        crit   = calc_surv_sd_c(LCinit(2),pmat,T_LCx(1),F_LCx(1)); % calculate criterion from upper range
    end % this makes sure that there is a range in which there IS an LCx,t for the first time point and highest effect level
    
    for i_e = 1:length(F_LCx) % run through effect levels requested
        for i_t = 1:length(T_LCx) % run through time points for the LCx
            LCx(i_t,i_e) = setup_rootfind(LCinit,1,pmat,T_LCx(i_t),F_LCx(i_e)); % now use wrapper function to find the LCx,t with <fzero> or <zero>
            LCinit = [pmat(GLO.mw,1) 1.01*LCx(i_t,i_e)]; % create new starting *range* (vector with 2 values)
            % Make use of the fact that LCx always decreases in time, and cannot be lower than <mw>!
        end
        LCinit = [pmat(GLO.mw,1) 1.01*LCx(1,i_e)];
        % Before we go to the next effect level, set the upper range for
        % LCx to the first value in this run: since we arranged <F_LCx>
        % from high to low x, the next first value will be lower than this
        % one!
    end
end

%% BLOCK 4. Make a plot of LCx versus time
% Start making a plot, even though we don't have the CIs yet. This gives
% the user something to look at while the calculations proceed.

plotcol= 'ygcrbm'; % six different colours for different x (should be enough)
figh   = make_fig(1,1,2); % create figure of correct size (one panel)
h_ax   = gca; % handle to current axes system
set(h_ax,'LineWidth',1,ft.name,ft.ticks) % adapt axis formatting to the standard font format
hold on % hold figure even when we add data with new plot commands

h_ln = nan(1,length(F_LCx_plot));  % initialise to suppress warning in Matlab editor
L    = cell(1,length(F_LCx_plot)); % initialise to suppress warning in Matlab editor
for i_e = 1:length(F_LCx_plot) % run through all effect levels
    % Note that we're only plotting a subset of the effect percentages.
    [~,loc_plot] = ismember(F_LCx_plot(i_e),F_LCx); % find location of this effect level in total set (index in vector)
    h_ln(i_e) = plot(T_LCx,LCx(:,loc_plot),'ko-','LineWidth',1,'MarkerFaceColor',plotcol(i_e)); % plot LCx over time, and collect the handle for the plotted points
    L{i_e}    = ['LC',num2str(100*F_LCx_plot(i_e))]; % collect legend text for this effect level
end

h_leg = legend(h_ln,L); % create a legend
set(h_leg,ft.name,ft.legend); % put it in the standard font format
h_tmp = xlabel('time (days)'); set(h_tmp,ft.name,ft.label) % add an x-axis label
h_tmp = ylabel(['LCx (',GLO.concunit,')']); set(h_tmp,ft.name,ft.label) % add a y-axis label
xlim([0 Tend])            % set limits of x-axis
ylim([0 1.1*max(LCx(:))]) % set limits of y-axis

drawnow % empty plot buffer if there's something in it (makes sure plot is not stalled until calculations finish)

%% BLOCK 5. Calculate CIs on the LCx,t.
% Note that I do not collect all LCx,t values for all parameter sets in the
% sample. Instead, I only update the min-max for each LCx,t value.

% BLOCK 5.1. Initial things.
n_sets  = size(sample,1); % number of sets in the sample
ind_log = pmat(:,5)==0;   % index to the parameters that are on log scale (vector with logical indexing)

% Use the fancy progress bar from Matlab
if sdit == 1
    f = waitbar(0,sprintf('Fraction of %1.0f sets calculated',n_sets),'Name','Calculating CIs on LCx for SD');
else
    f = waitbar(0,sprintf('Fraction of %1.0f sets calculated',n_sets),'Name','Calculating CIs on LCx for IT');
end

for i_n = 1:n_sets % run through all sets in the sample
    
    waitbar(i_n/n_sets,f) % update waitbar
    
    % BLOCK 5.2. Replace parameters in <pmat> with the values in the next
    % set of the sample.
    pmat(:,1) = sample(i_n,:); % replace values in <pmat> with the next set in the sample
    % Put parameters that need to be on log scale back on normal scale (as
    % the SD/IT functions require normal scale).
    if sum(ind_log)>0
        pmat(ind_log,1) = 10.^(pmat(ind_log,1));
    end
    % Note: <pmat> is on normal scale here, but the sample contains the
    % value on a log scale, if a parameter is fitted on log scale.
    
    % BLOCK 5.3. For IT ...
    if sdit == 2 % for IT, we can use a direct calculation
        for i_e = 1:length(F_LCx) % run through effect levels requested
            for i_t = 1:length(T_LCx) % run through times requested
                [~,LCtmp] = calc_surv_it(1,pmat,T_LCx(i_t),F_LCx(i_e),[]); % calculate LCx
                LClo(i_t,i_e) = min([LClo(i_t,i_e) LCtmp],[],2); % keep smallest value of current lowest and new LCx (matrix)
                LChi(i_t,i_e) = max([LChi(i_t,i_e) LCtmp],[],2); % keep largest value of current largest and new LCx (matrix)
            end
        end
        
    % BLOCK 5.4. For SD ...
    else % for SD, some more work is needed; need to use root finding
        
        % NEW: make sure that LCinit is always a range.
        LCinit = [pmat(GLO.mw,1)/10 pmat(GLO.mw,1)]; % as search range for LCx,t start with <mw>/10 and <mw>
        % Using <mw> is fine for lower end of range but not for upper, so find
        % a better value for the upper range by trial and error. I do this
        % by shifting this initial range up until it covers the sign change
        % of <crit>.
        crit = 1; % first one will be positive
        while crit > 0
            LCinit = LCinit * 10; % shift lower and upper range by factor of 10
            crit = calc_surv_sd_c(LCinit(2),pmat,T_LCx(1),F_LCx(1)); % calculate criterion from upper range
        end % this makes sure that there is a range in which there IS an LCx,t for the first time point and highest effect level
        
        for i_e = 1:length(F_LCx) % run through effect levels requested
            for i_t = 1:length(T_LCx) % run through time points for the LCx
                LCtmp  = setup_rootfind(LCinit,1,pmat,T_LCx(i_t),F_LCx(i_e)); % now use wrapper function to find the LCx,t with <fzero> or <zero>
                LCinit = [pmat(GLO.mw,1) 1.01*LCtmp]; % update <LCinit>, and now make it a range for <fzero> to search
                % Make use of the fact that LC50s always decrease in time, and cannot be lower than <mw>!
                LClo(i_t,i_e) = min([LClo(i_t,i_e) LCtmp],[],2); % smallest value of current lowest and new LCx
                LChi(i_t,i_e) = max([LChi(i_t,i_e) LCtmp],[],2); % largest value of current largest and new LCx

                if i_t == 1 % for the first time point ...
                    LCinit_new = LCinit; % collect the updated search range as start value for next effect level
                end
            end
            LCinit = LCinit_new;
            % For the next effect level, set the upper range for LCx to the
            % first value in this run: since we arranged <F_LCx> from high
            % to low x, the next first value will be lower than this one!
        end
    end
        
end

delete(f) % delete progress bar again

%% BLOCK 6. Put results on screen.

disp(['Results table for LCx,t (',GLO.concunit,'), with 95% CI'])
if sdit == 1
    disp('  Special case: stochastic death (SD)');
else
    disp('  Special case: individual tolerance (IT)');
end
disp('==========================================================================================')

fprintf('time (d)       ')
for i_e = 1:length(F_LCx) % run through effect levels requested
    fprintf('LC%1.0f                        ',F_LCx(i_e)*100) % print table headers
end
fprintf('\n')
disp('==========================================================================================')

for i_t = 1:length(T_LCx) % run through times requested
    fprintf('%4g  ',T_LCx(i_t))
    for i_e = 1:length(F_LCx) % run through effect levels requested
        fprintf('%#7.4g (%#7.4g - %#7.4g) ',LCx(i_t,i_e),LClo(i_t,i_e),LChi(i_t,i_e))
    end
    fprintf('\n')
end
disp('==========================================================================================')

%% BLOCK 7. Finish the plot off with CIs.

figure(figh) % make the figure we already have current, just to make sure we are plotting in the right place

for i_e = 1:length(F_LCx_plot) % run through all effect levels
    
    % Note that we're only plotting a subset of the effect percentages.
    [~,loc_plot] = ismember(F_LCx_plot(i_e),F_LCx); % find location of this effect level in total set (index in vector)
    % Plot CIs with a dotted line.
    plot(T_LCx,LClo(:,loc_plot),'k:','LineWidth',1)
    plot(T_LCx,LChi(:,loc_plot),'k:','LineWidth',1)
    
    % Little trick to fill the area between the two curves, to obtain a
    % coloured confidence band.
    t2  = [T_LCx,fliplr(T_LCx)]; % make a new time vector that is old one, plus flipped one
    Xin = [LClo(:,loc_plot);flipud(LChi(:,loc_plot))]; % do the same for the plot line, hi and lo
    fill(t2,Xin,plotcol(i_e),'LineStyle','none','FaceAlpha',0.25) % and fill this object
    
end

for i_e = 1:length(F_LCx_plot) % run through effect levels requested
    uistack(h_ln(i_e),'top'); % use the previously collected handles to the plotted points to place them on top of the CI bands
end

set(h_ax, 'Layer', 'top') % make sure the axes are visible (the yellow object wants to obscure it)
if ~isnan(LChi(1,1)) % if we did not fit anything, we have no CIs, and no LChi
    ylim([0 1.1*max(LChi(:))]) % otherwise, reset limits on y-axis using upper CI
end

% Make an overall title for the plot with the filename, and save it.
axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off','Visible','off','Units','normalized', 'clipping','off');
h_txt = text(0.5, 1,['File: ',savestr],'Interpreter','none','HorizontalAlignment','center','VerticalAlignment', 'top');
set(h_txt,ft.name,ft.text);
drawnow
save_plot(figh,['output_report',filesep,GLO.basenm,filesep,savestr],GLO.save_format) % save LCx,t plot in output folder
snapnow % for publishing Matlab scripts: make plot appear in the context of the screen output

disp(' ')
diary off

