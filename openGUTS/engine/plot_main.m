function plot_main(data_all,sdit,flag,mf)

% Usage: plot_main(data_all,sdit,flag,mf)
%
% This is the main plotting function. It plots three rows of panels. First
% row is exposure scenarios, second row damage vs. time, third row survival
% vs. time. 
%
% Inputs
% <data_all>    the data object to be plotted (structured cell array)
% <sdit>        switch for SD (1) or IT (2)
% <flag>        a two-value vector, see below
% <mf>          multiplication factor (matrix, only used for
%               predictions); in general, this will be the LPx factors 
%               calculated by <calc_lpx> (for each data set, for several
%               values of x)
% 
% <flag>: first element of vector indicates the stage of the analysis: 
% 1) called from the input data/calibration stage 
% 2) called from the validation stage
% 3) called from the predictions stage
%
% <flag>: second element of vector indicates whether model curves need to be
% calculated/shown: 
% 0) inspect data (shows only survival data and exposure scenario) 
% 1) show model curves as well (for survival and damage), with CIs for stage 1 and 2 
% 2) for predictions: force showing model curves with CIs (terribly slow for FOCUS profiles, so this is standardly not used at the moment)
% 3) for batch processing of LPx: make an invisible plot without CIs (saved but not shown on screen)
% 
% Author     : Tjalling Jager
% Date       : March 2023
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

flag_stage = flag(1); % flag for which stage of the workflow we are in (calibration/validation/prediction)
flag_model = flag(2); % flag for whether to plot with model curves (1-3) or without (0), see above

ft = GLO.ft; % retrieve global for font specs in plots

%% BLOCK 1. Multiplot of scenarios, damage and survival versus time.
% A multiplot is made with endpoint vs. time for each treatment separately.

nr_treat = data_all.nr_treat; % extract vector with the number of scenarios in each data set

for i_d = 1:data_all.nr_sets % run through all data sets
    
    % BLOCK 1.1. Initial things
    
    % Extract useful info from the object for this data set.
    number_scen = nr_treat(i_d); % extract the number of treatments in this data set
    scen_name   = data_all.name; % extract the scenario names (cell array with all names in this data set)
    
    % Extract maximum time and concentration for scaling plots consistently.
    t_max  = data_all.time{i_d}(end);       % last time point in this data set
    c_max  = 1.02*max(data_all.c_max{i_d}); % maximum concentration in this data set
    % Add a little extra factor for <c_max> to make it fit better in the figure.
    d_max  = 0;                   % this will collect the maximum damage, for scaling plots
    h_pl_d = cell(number_scen,1); % create empty cell array to collect handles to damage plots (to allow scaling after all is plotted)
    
    % For predictions, there is only 1 treatment per data set. We will
    % re-use the option to plot various treatments to plot various
    % multiplication factors!
    if flag_stage == 3 && ~isempty(mf)   % then we plot predictions with multiplication factors (MF)
        c_max = c_max * max(mf(i_d,:));  % rescale <c_max> using the highest MF supplied for this data set
        number_scen = length(mf(i_d,:)); % every value in <mf> will become a new scenario, so make a new <number_scen> as the nr of columns of <mf>
    end
    
    if flag_model == 3 % flag 3 is only used for LPx batch processing, where we don't want to see figure windows popping up all the time
        figh = make_fig(3,number_scen,0); % create an INVISIBLE figure window of correct size
    else
        figh = make_fig(3,number_scen,1); % create a figure window of correct size
    end
    hold on
    
    % Create a time vector for the model curves. This creates a long time
    % vector for smooth plots (this is amended with the time vector for the
    % survival data sets in the next loop).
    t_mod = (linspace(0,t_max,max(GLO.t_plot(1),t_max*GLO.t_plot(2))))'; % this vector is used for plotting model results (to create smooth curves)
    % Note: this creates a plotting point at every hour, or at least a 100
    % ones (defaults set in <initial_setup>). For FOCUS scenarios, more
    % than 100 points are needed to properly show damage for very peaky
    % scenarios and fast kinetics.
    
    for i_c = 1:number_scen % run through all scenarios for data set <i_d>
        
        % BLOCK 1.2. Extract the data for survivors and the exposure
        % scenario for this treatment only.
        if flag_stage == 3 && ~isempty(mf)    % then we plot predictions with multiplication factors
            data_surv = data_all.surv{i_d,1}; % always take first treatment (there is only one) (matrix)
            ev_ic     = data_all.scen{i_d,1}; % always take first exposure scenario (there is only one) (matrix)
            ev_ic(:,[2 3]) = ev_ic(:,[2 3]) * mf(i_d,i_c); % apply next multiplication factor to concentrations and slopes of exposure-time matrix (2nd and 3rd column)
        else
            data_surv = data_all.surv{i_d,i_c}; % take survivors for treatment <i_c> from data set <i_d> (matrix)
            ev_ic     = data_all.scen{i_d,i_c}; % take scenario for treatment <i_c> from data set <i_d> (matrix)
        end
        
        % Also ensure that the actual time points for the survival data are
        % in the long time vector. This is handy so we can extract the
        % modelled survival probabilities for the comparison to observed
        % ones, without needing to interpolate.
        t_mod = unique([t_mod;data_surv(:,1)]); % add long time vector and the one for this treatment, keep unique values and sort

        % BLOCK 1.3. Calculate model curves.
        if flag_model > 0 % we're showing fits, so we need model curves
            
            % For validation, we can also use the <hb> estimated from the
            % validation controls. However, then we need to tell
            % <calc_conf> that we want to fix <hb> to a new value for the
            % predicted output. For predictions, we need to fix <hb=0>. For
            % calibration, this is not needed. When we fix <hb> to the
            % controls that is already taken care of in <startgrid>.
            % Therefore, the sample already only contains this fixed value.
            hb_fix    = -1; % by default, <hb_fix> set to -1
            flag_conf = 1;  % by default, calculate confidence intervals
            if flag_stage == 2 && GLO.fix_hb_val == 1 % if we are in validation stage and want to fix <hb>
                if data_all.hb_hat == -1 % if there is no <hb> in the data object
                    error('Requested fixing background hazard to control for validation, but no controls were found.')
                else
                    hb_fix = data_all.hb_hat; % make <hb_fix> equal to the value in the current data object
                end
            end
            if flag_stage == 3 % for predictions, always fix <hb> to zero
                hb_fix = 0;
                if flag_model ~= 2 % unless we specifically ask for it ...
                    flag_conf = 0; % don't calculate CIs for predictions, that's much too slow
                end
            end
            
            % Use <calc_conf> to calculate the model curves for survival
            % and damage, and bounds on the model curves, from the saved
            % calibration file (two three-column matrices).
            [model_p,model_dam] = calc_conf(t_mod,ev_ic,sdit,hb_fix,flag_conf);
            d_max = max([d_max;model_dam(:)]); % simply take max of previous max and all new damage calculations (best and CIs)
            % Note: this will be used to scale the damage plots afterwards.
            
            if flag_stage ~= 3 % don't do this for predictions
                % Collect information for later observed-predicted plots and
                % goodness-of-fit statistics. Find and store predictions for
                % survival probability at time points in data set (with CI).
                [~,loct] = ismember(data_surv(:,1),t_mod); % location (vector) in the long time vector of the time points for this treatment in the data set
                
                % Collect model predictions with CIs for the time points in
                % the survival data for this treatment; simply add them as
                % three new columns to the data object. Last two columns
                % are the CIs on predicted survival probability.
                data_all.surv{i_d,i_c} = [data_all.surv{i_d,i_c} , model_p(loct,:)];
            end
        end
        
        %% BLOCK 1.4. First row of panels to plot is the exposure scenarios
        % The exposure scenario is always plotted, for each analysis stage.
        
        set(0,'CurrentFigure',figh); % make sure figure is current (don't use <figure> here as will make it visible for batch work as well!)
        h_pl = subplot(3,number_scen,i_c); % create subplot in correct position and return handle
        hold on % keep the plotted data when there is a subsequent plotting command
        set(h_pl,'LineWidth',1,ft.name,ft.ticks) % adapt axis formatting to the standard setting
        
        if number_scen > 1 % shrink white space in plot when there are more than 2 rows (Matlab specific, to make multiplots look better)
            p = get(h_pl,'position'); % get position of the subplot
            p([3 4]) = p([3 4])*1.1;  % add 10 percent to width and height
            set(h_pl, 'position', p); % and set it to that value
        end
        
        if flag_stage == 3 % then we plot predictions, and need the multiplication factor in the title
            if isempty(mf) % it is empty for inspection plots, and then we just plot MF=1
                h_tmp = title(['set: ',num2str(i_d),', ',scen_name{i_d,1},', MF=1'],'Interpreter','none');
            else
                h_tmp = title(['set: ',num2str(i_d),', ',scen_name{i_d,1},', MF=',num2str(round(mf(i_d,i_c),1))],'Interpreter','none');
            end
        else % otherwise, just put the number of the data set and the name of the treatment in the title
            h_tmp = title(['set: ',num2str(i_d),', ',scen_name{i_d,i_c}],'Interpreter','none');
        end
        set(h_tmp,ft.name,ft.title)
        
        if i_c == 1 % only put y labels on first column
            h_tmp = ylabel(['exposure conc. (',GLO.concunit,')']); set(h_tmp,ft.name,ft.label)
        else
            set(h_pl,'YTickLabel',[]); % remove tick labels on y-axis for other columns
        end 
        % No x-labels for first row, unless making inspection plots for predictions (which have no survival data).
        if flag_stage == 3 && flag_model == 0
            h_tmp = xlabel('time (days)'); set(h_tmp,ft.name,ft.label) % make label and ticks for inspecting prediction exposure scenarios
        else            
            set(h_pl,'XTickLabel',[]); % remove tick labels on x-axis otherwise
        end
        % Next, all plots within this data set will be scaled in the same manner.
        xlim([0 t_max]) % set limit of x-axis
        ylim([0 c_max]) % set limit of y-axis
        
        % Plot exposure scenario This is done here in a bit of an awkward
        % manner as I allowed multiple entries for the same time point to
        % indicate a sudden switch in exposure concentration. Since these
        % double entries are removed in <prepare_data>, we need to
        % 'reconstruct' them. This is rather inefficient, especially for
        % FOCUS profiles where it is not needed (there should be no double
        % entries in the profiles). However, I don't think there is a
        % pressing need to create shortcuts (as an alternative, we can
        % consider also saving the un-modified scenario in the data object,
        % saving a flag in <prepare_data> to let the plot routine know that
        % there were no double entries, or even re-load the input text
        % file).
        scen_tmp = []; % create an empty matrix to catch the reconstructed profile (two columns as in input)
        for i_t = 1:size(ev_ic,1)-1   % run through all time steps, not last
            c_start = ev_ic(i_t,2);   % start concentration in this interval
            c_next  = ev_ic(i_t+1,2); % start concentration of next interval
            slope   = ev_ic(i_t,3);   % slope in this interval
            delta_t = ev_ic(i_t+1,1) - ev_ic(i_t,1); % length of this interval
            c_end   = c_start + slope * delta_t; % end concentration of this interval
            
            % Add one or two entries to the scenario matrix.
            scen_tmp = [scen_tmp ; ev_ic(i_t,1) c_start]; %#ok<AGROW> always add start time and concentration
            if c_next ~= c_end % but if there is a step-change at the end of this interval ....
                scen_tmp = [scen_tmp ; ev_ic(i_t+1,1) c_end]; %#ok<AGROW> add extra entry with end concentration
            end
        end
        if ev_ic(end,1) ~= scen_tmp(end,1) % if final time point in <scen_tmp> is not the last one in <ev_ic> ...
            scen_tmp = [scen_tmp ; ev_ic(end,[1 2])]; %#ok<AGROW> add last concentration of the scenario
            % This will happen almost always. Exception is when someone
            % enters an exposure scenario with a double entry at the end of
            % the survival data set (which would be silly).
        end
                
        area(scen_tmp(:,1),scen_tmp(:,2),'FaceColor','b'); % plot concentration profile as a filled area (option 'FaceAlpha',0.25 removed)
        
        %% BLOCK 1.5. Second row of multiplot is for damage predictions.
        % This is not needed for inspection plots (there is no damage in
        % the data set).
        
        if flag_model > 0 % not for inspection plots
            h_pl = subplot(3,number_scen,1*number_scen+i_c); % create subplot in correct position and return handle
            hold on % keep the plotted data when there is a subsequent plotting command
            set(h_pl,'LineWidth',1,ft.name,ft.ticks) % adapt axis formatting to standard setting
            h_pl_d{i_c} = h_pl; % collect handle into a cell array, so that we can scale all y-axes later (after we know the maximum)
            
            if number_scen > 1 % shrink white space in plot when there are more than 2 rows
                p = get(h_pl,'position'); % get position of the subplot
                p([3 4]) = p([3 4])*1.1; % add 10 percent to width and height
                set(h_pl, 'position', p); % and set it to that value
            end
            
            if i_c == 1 % only put y labels on first column
                h_tmp = ylabel(['scaled damage (',GLO.concunit,')']); set(h_tmp,ft.name,ft.label)
            else % otherwise y-axes are removed
                set(h_pl,'YTickLabel',[]); % remove tick labels on x-axis
            end
            set(h_pl,'XTickLabel',[]); % remove tick labels on x-axis; no x-labels for second row
            
            % Next, all plots within this data set will be scaled in the same manner.
            xlim([0 t_max]) % set limit of x-axis
            ylim([0 c_max]) % set temporary limit of y-axis (rescaling later)
            
            % Plot max and min model lines as calculated by <calc_conf>.
            % Note that if no CIs are (or can be) calculated, the 2nd and
            % 3rd column of the matrix <model_dam> will only contain NaNs
            % and CIs are not plotted).
            plot(t_mod,model_dam(:,2),'k:','LineWidth',1) % plot each lower confidence line as broken line
            plot(t_mod,model_dam(:,3),'k:','LineWidth',1) % plot each upper confidence line as broken line
            
            % Little trick to fill the area between the two curves, to
            % obtain a coloured confidence interval as a band.
            t2  = [t_mod;flipud(t_mod)]; % make a new time vector that is old one, plus flipped one
            Xin = [model_dam(:,2);flipud(model_dam(:,3))]; % do the same for the plot line, hi and lo
            fill(t2,Xin,'g','LineStyle','none','FaceAlpha',0.25) % and fill this object
            
            plot(t_mod,model_dam(:,1),'k-','LineWidth',1) % plot the model line for the best-fitting parameters on top
        end
        
        %% BLOCK 1.6. Third row of multiplot is for survival data/predictions.
        
        if ~(flag_stage == 3 && flag_model == 0) % don't plot survival for the inspection plots of the prediction stage
            
            h_pl = subplot(3,number_scen,2*number_scen+i_c); % create subplot in correct position and return handle
            hold on % keep the plotted data when there is a subsequent plotting command
            set(h_pl,'LineWidth',1,ft.name,ft.ticks) % adapt axis formatting to standard
            
            if number_scen > 1 % shrink white space in plot when there are more than 2 rows
                p = get(h_pl,'position'); % get position of the subplot
                p([3 4]) = p([3 4])*1.1;  % add 10 percent to width and height
                set(h_pl, 'position', p); % and set it to that value
            end
            
            if i_c == 1 % only put y labels on first column
                h_tmp = ylabel('survival prob. (-)'); set(h_tmp,ft.name,ft.label)
            else % for other columns, remove them
                set(h_pl,'YTickLabel',[]); % remove tick labels on x-axis
            end
            h_tmp = xlabel('time (days)'); set(h_tmp,ft.name,ft.label) % we are in the last row, so use x-labels for each column
            
            % Next, all plots within this data set will be scaled in the same manner.
            xlim([0 t_max]) % set limit of x-axis
            ylim([0 1.02])  % set limit of y-axis
            
            if flag_model ~= 0 % then we plot model curves with CIs
                % Plot max and min model lines as calculated by
                % <calc_conf>. Note that if no CIs are (or can be)
                % calculated, the 2nd and 3rd column of the matrix
                % <model_p> will only contain NaNs and CIs are not
                % plotted).
                plot(t_mod,model_p(:,2),'k:','LineWidth',1) % plot each lower confidence line as broken line
                plot(t_mod,model_p(:,3),'k:','LineWidth',1) % plot each upper confidence line as broken line
                
                % Little trick to fill the area between the two curves, to
                % obtain a coloured confidence interval as a band.
                t2  = [t_mod;flipud(t_mod)]; % make a new time vector that is old one, plus flipped one
                Xin = [model_p(:,2);flipud(model_p(:,3))]; % do the same for the plot line, hi and lo
                fill(t2,Xin,'g','LineStyle','none','FaceAlpha',0.25) % and fill this object
                
                plot(t_mod,model_p(:,1),'k-','LineWidth',1) % plot the model line for the best-fitting parameters
            end
            
            % Next, plot data points with Wilson score CI.
            if flag_stage ~= 3 % don't do this for predictions, which have no survival data
                if flag_model == 0 % for inspection plot, use a connecting dotted line
                    plot(data_surv(:,1),data_surv(:,4),'ko:','MarkerFaceColor','k') % plot the data points for survival probability, with dotted line
                else % but not for other plots
                    plot(data_surv(:,1),data_surv(:,4),'ko','MarkerFaceColor','k') % plot the data points for survival probability, without dotted line
                end
                % Plot Wilson score CI as simple vertical line through the data points.
                for i_t = 1:length(data_surv(:,1)) % run through time points in this survival data set
                    plot([data_surv(i_t,1) data_surv(i_t,1)],[data_surv(i_t,5) data_surv(i_t,6)],'k-') % and plot Wilson score
                end
            end
        end
    end
    
    % BLOCK 1.7. Make an overall title for the plot with the filename, and
    % save it.
    
    % Scale all damage plots to new-found maximum, but only do that when
    % the new maximum is substantially lower than the max exposure conc.
    if d_max < 0.75 * c_max && flag_model > 0 % not when difference is small and not for inspection plots
        for i_c = 1:number_scen % run through the scenarios
            set(h_pl_d{i_c},'Ylim',[0 1.02*d_max]); % scale the y-axis to new maximum
        end
    end
        
    % First, create a name for the plot (used for saving), using the <basenm>
    savestr = GLO.basenm;
    switch flag_stage % the name should contain the stage of the analysis
        case 1
            savestr = [savestr,'_cal']; %#ok<AGROW>
        case 2
            savestr = [savestr,'_val']; %#ok<AGROW>
        case 3
            savestr = [savestr,'_pred_',scen_name{i_d,1}]; %#ok<AGROW> for predictions, add scenario name to the filename
    end
    if flag_model == 0 % inspect plots need to be indicated (but don't really need to be saved, I guess)
        savestr = [savestr,'_inspect_',num2str(i_d)]; %#ok<AGROW>
    else % when model has been plotted, also include whether it is SD or IT
        if sdit == 1
            savestr = [savestr,'_SD']; %#ok<AGROW>
        else
            savestr = [savestr,'_IT']; %#ok<AGROW>
        end
        if flag_stage ~= 3 % also add data set number unless we're doing predictions (for which we already have the filename in the <savestr>)
            savestr = [savestr,'_',num2str(i_d)]; %#ok<AGROW>
        end
    end
    
    % Create an overall title.
    axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off','Visible','off','Units','normalized', 'clipping','off');
    h_txt = text(0.5, 1,['File: ',savestr],'Interpreter','none','HorizontalAlignment','center','VerticalAlignment', 'top');
    set(h_txt,ft.name,ft.text); % use standard formatting for this header
    drawnow % empty the plotting buffer
    save_plot(figh,['output_report',filesep,GLO.basenm,filesep,savestr],GLO.save_format) % save plot in output folder
    snapnow % for publishing Matlab scripts: make plot appear in the context of the screen output
end

%% BLOCK 2. Predicted versus observed plots.
% Plots are made for survival probability and numbers of deaths.
% Additionally, the goodness-of-fit metrics are calculated. That is done
% here, as the plotting routine already calculated the modelled survival
% probability (with CI) for us. No need to do it again elsewhere.

if flag_stage ~= 3 && flag_model ~= 0 % only do this when calibrating/validating and when model curves have been calculated
    
    markers = {'s','^','o','d','v','p','>','<'}; % a complete set of 8 markers
    plotcol = 'grbckw'; % 6 colors should be enough ... (magenta thrown out as I don't like it)
    
    for i_d = 1:data_all.nr_sets % run through all data sets
        
        % BLOCK 2.1. Make a new observed-predicted survival probability plot for each data set.
        figh = make_fig(1,3,2); % initiate a two-panel figure plus a legend panel
        h_pl = subplot(1,3,1);  % first panel, get handle to current axis
        hold on
        set(h_pl,'LineWidth',1,ft.name,ft.ticks) % adapt axis formatting to standard setting
        h_tmp = xlabel('observed survival prob.'); set(h_tmp,ft.name,ft.label) % create an x-label
        h_tmp = ylabel('predicted survival prob.'); set(h_tmp,ft.name,ft.label) % create a y-label
        h_tmp = title(['set: ',num2str(i_d)]); set(h_tmp,ft.name,ft.title) % create a title with data set number
        
        h_plotpts = []; % clear this set of handles, as it creates problems when a next data set has less treatments
        surv_max  = 1;  % this is for scaling the plot (would only be needed if we have to move to survivor *numbers* at some point)
        for i_c = 1:nr_treat(i_d) % run through all treatments for data set <i_d>
            data_surv = data_all.surv{i_d,i_c}; % extract the extended survivor matrix for this treatment
            obs_surv  = data_surv(:,4); % only use observed survival probability
            pred_surv = data_surv(:,end-2:end); % extract predicted survival probabilities (collected in the plotting part as last three columns of the survival cell)
            
            % Catch cases where a user entered only survival data at t=0
            % (and missing data for other time points). Need to catch that
            % as otherwise an error is produced on plotting.
            if numel(obs_surv) < 2
                continue % go to next iteration of for loop
            end
            
            % Plot predicted vs. observed survival probabilities. We only
            % have 8 markers, but we'll change colour when we run out of
            % them.
            nextmarker = ceil(i_c/length(plotcol)); % note: <ceil> rounds up, so this makes sure that marker 1 is used for the first colors, and then marker 2, etc.
            nextcolor  = i_c - (nextmarker-1) * length(plotcol); % take the correct color (restart at 1 for each new marker)
            h_plotpts(i_c) = plot(obs_surv(2:end),pred_surv(2:end,1),['k',markers{nextmarker}],'MarkerFaceColor',plotcol(nextcolor),'MarkerSize',8); %#ok<AGROW>
            % Note: the 2:end makes sure that t=0 is not plotted. The
            % handles are now collected in h_plotpts (they are only used to
            % make sure that the plotting symbols are on top of all error
            % bars).
            
            % Plot error bars on the predicted survival probability as vertical line.
            for i_r = 2:length(obs_surv) % start at 2, to avoid plotting bars at t=0
                plot([obs_surv(i_r) obs_surv(i_r)],[pred_surv(i_r,2) pred_surv(i_r,3)],'k-','LineWidth',1)
            end
        end

        plot([0 surv_max],[0 surv_max],'k--','LineWidth',1) % plot the 1:1 line
        xlim([0 surv_max]) % scale the x-axis
        ylim([0 surv_max]) % scale the y-axis
        
        % Plot plus and minus 20% in the plot (roughly equivalent to what
        % is plotted in the EFSA opinion).
        plot([0 surv_max],[0 surv_max]-0.2*surv_max,'k:','LineWidth',1) % plot x% beneath the 1:1 line
        plot([0 surv_max],[0 surv_max]+0.2*surv_max,'k:','LineWidth',1) % plot x% above the 1:1 line

        uistack(h_plotpts,'top'); % place symbols on top again (otherwise they are covered by the error bars)
        
        % BLOCK 2.2. Make a new observed-predicted deaths plot for each data set.
        h_pl = subplot(1,3,2); % second panel, get handle to current axis
        hold on
        set(h_pl,'LineWidth',1,ft.name,ft.ticks) % adapt axis formatting to standard settings
        h_tmp = xlabel('observed deaths/interval'); set(h_tmp,ft.name,ft.label) % create an x-label
        h_tmp = ylabel('predicted deaths/interval'); set(h_tmp,ft.name,ft.label) % create a y-label
        h_tmp = title(['set: ',num2str(i_d)]); set(h_tmp,ft.name,ft.title) % create a title with data set number
        
        h_plotpts    = []; % clear this set of handles, as it creates problems when a next data set has less treatments
        legend_label = []; % also clear legend labels (which we will collect now)
        death_max    = 0;  % this will catch the largest number of deaths, for scaling the plot
        for i_c = 1:nr_treat(i_d) % run through all treatments for this data set
            data_surv = data_all.surv{i_d,i_c}; % extract the extended survivor matrix for this treatment
            
            % Turn all predicted survival probabilities into death
            % probabilities, and subsequently in numbers (without CIs).
            pred_deaths = data_surv(1,2) * [-diff(data_surv(:,end-2));data_surv(end,end-2)]; % calculate number of deaths from survivors (last interval is after end of test)
            % Note: we forget about CIs for the absolute numbers of deaths.
            % Also note that the end-2 column is the predicted survival
            % probability as collected above in the plotting section.
            
            % For the max number of deaths, I take the max over the previous
            % max and all deaths plotted.
            death_max = max([death_max;data_surv(:,3);pred_deaths]); 
            
            % Plot predicted vs. observed deaths per interval. We only
            % have 8 markers, but we'll change colour when we run out of
            % them.
            nextmarker = ceil(i_c/length(plotcol)); % note: <ceil> rounds up, so this makes sure that marker 1 is used for the first colors, and then marker 2, etc.
            nextcolor  = i_c - (nextmarker-1) * length(plotcol); % take the correct color (restart at 1 for each new marker)
            h_plotpts(i_c) = plot(data_surv(:,3),pred_deaths(:,1),['k',markers{nextmarker}],'MarkerFaceColor',plotcol(nextcolor),'MarkerSize',8); %#ok<AGROW>
            legend_label{i_c} = scen_name{i_d,i_c}; %#ok<AGROW> collect names to make the legend later
        end
        
        plot([0 death_max],[0 death_max],'k--','LineWidth',1) % plot 1:1 line
        xlim([0 death_max]) % scale x-axis
        ylim([0 death_max]) % scale y-axis
        
        % Plot plus and minus 20% in the plot (roughly equivalent to what's plotted in the EFSA opinion)
        plot([0 death_max],[0 death_max]-0.2*death_max,'k:','LineWidth',1) % plot x% beneath the 1:1 line
        plot([0 death_max],[0 death_max]+0.2*death_max,'k:','LineWidth',1) % plot x% above the 1:1 line
        
        uistack(h_plotpts,'top'); % place symbols on top again (otherwise they are covered by the error bars)
        
        % BLOCK 2.3 Use a third panel for the legend. This is a rather
        % tricked-out Matlab-specific code to make a legend and place it in
        % a separate figure panel (so it does not hide the plotted points).
        h_annot  = subplot(1,3,3); % create third panel and extract its handle
        h1       = gca; % remember the current axis handle
        dim      = h_annot.Position; % the position of the subplot is also the correct place for the textbox
        
        h_leg  = legend(h_plotpts,legend_label); % create a legend from the plot handles and the collected labels
        dimleg = h_leg.Position; % position of last legend made
        h_leg.Position = [dim(1) dim(2)+dim(4)-dimleg(4) dimleg(3:4)]; % move it to the new subplot and set it to top-left of panel
        delete(h1) % throw away the axes of the plot, the legend is enough
        
        % BLOCK 2.4. Make an overall title for the plot with the filename, and save it
    
        % First, create a name for the plot (used for saving), using the
        % <GLO.basenm>.
        savestr = GLO.basenm;
        switch flag_stage
            case 1
                savestr = [savestr,'_obspred_cal']; %#ok<AGROW>
            case 2
                savestr = [savestr,'_obspred_val']; %#ok<AGROW>
            case 3
                % does not occur
        end
        if sdit == 1 % also include whether it is SD or IT and number of data set
            savestr = [savestr,'_SD_',num2str(i_d)]; %#ok<AGROW>
        else
            savestr = [savestr,'_IT_',num2str(i_d)]; %#ok<AGROW>
        end
        
        % Make an overall title for the plot with the filename, and save it
        axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off','Visible','off','Units','normalized', 'clipping','off');
        h_txt = text(0.35, 1,['File: ',savestr],'Interpreter','none','HorizontalAlignment','center','VerticalAlignment', 'top');
        set(h_txt,ft.name,ft.text); % use standard formatting for this header
        drawnow % empty plotting buffer
        save_plot(figh,['output_report',filesep,GLO.basenm,filesep,savestr],GLO.save_format) % save plot in output folder
        snapnow % for publishing Matlab scripts: make plot appear in the context of the screen output
    end
    
    %% BLOCK 2.5. Calculate goodness-of-fit measures
    % Do this once for all data sets together. We skip the PPC% as it is
    % not defined when we don't calculate the sampling error.
    
    coll_p     = []; % used to collect all observed and predicted survival probabilities
    coll_s     = []; % used to collect all observed and predicted survivor numbers
    coll_p_end = []; % used to collect all observed and predicted survival probabilities at the end of the test

    % =====================================================================
    % NEW after v.1.0: ALSO CALCULATE NSE AND NRMSE WITH t=0 IN THERE
    coll_p0    = []; % used to collect all observed and predicted survival probabilities
    coll_s0    = []; % used to collect all observed and predicted survivor numbers
    % =====================================================================
    
    for i_d = 1:data_all.nr_sets  % run through all data sets
        for i_c = 1:nr_treat(i_d) % run through all treatments for this data set <i_d>
            
            data_surv = data_all.surv{i_d,i_c}; % extract extended survival matrix for this treatment
            
            % Make a long matrix with all observed and predicted survival
            % probabilities. Note that the end-2 column is the predicted
            % survival probability added above in the plotting part.
            coll_p = [coll_p ; data_surv(2:end,4) data_surv(2:end,end-2)]; %#ok<AGROW> append observed and predicted surv. prob.
            % Also make a long matrix with all observed and predicted survivor *numbers*.
            coll_s = [coll_s ; data_surv(2:end,2) data_surv(1,2) * data_surv(2:end,end-2)]; %#ok<AGROW> append observed and predicted surv. numbers
            % Note: I exclude the observed-predicted values at t=0 (as they are perfectly fitted anyway).
            
            % =====================================================================
            % NEW after v.1.0: ALSO CALCULATE NSE AND NRMSE WITH t=0 IN THERE
            coll_p0 = [coll_p0 ; data_surv(:,4) data_surv(:,end-2)]; %#ok<AGROW> append observed and predicted surv. prob.
            coll_s0 = [coll_s0 ; data_surv(:,2) data_surv(1,2) * data_surv(:,end-2)]; %#ok<AGROW> append observed and predicted surv. numbers
            % =====================================================================
            
            % And also make a long matrix for the probabilities at the end.
            % Also store data set number and treatment number with the
            % probabilities for this one.
            coll_p_end = [coll_p_end; data_surv(end,4) data_surv(end,end-2) i_d i_c]; %#ok<AGROW> append observed and predicted surv. prob. at the end of the test
            
        end
    end
    
    % Do an r-square calculation on the survival *probabilities*. This is
    % also known as the 'model efficiency'.
    ssq_fit  = sum((coll_p(:,1)-coll_p(:,2)).^2); % sum-of-squares of the residuals
    ssq_tot  = sum((coll_p(:,1)-mean(coll_p(:,1))).^2); % sum-of-squares of the total (compared to mean)
    r_square = 1 - ssq_fit / ssq_tot; % r-square
    
    % Do NRMSE as in EFSA document on survival *numbers*, expressed as percentage
    ssq_fit = sum((coll_s(:,1)-coll_s(:,2)).^2); % sum-of-squares of the residuals
    nrmse   = 100 * sqrt(ssq_fit/length(coll_s(:,1))) / mean(coll_s(:,1));
    
    % =====================================================================
    % NEW after v.1.0: ALSO CALCULATE NSE AND NRMSE WITH t=0 IN THERE
    ssq_fit   = sum((coll_p0(:,1)-coll_p0(:,2)).^2); % sum-of-squares of the residuals
    ssq_tot   = sum((coll_p0(:,1)-mean(coll_p0(:,1))).^2); % sum-of-squares of the total (compared to mean)
    r_square0 = 1 - ssq_fit / ssq_tot; % r-square
    ssq_fit   = sum((coll_s0(:,1)-coll_s0(:,2)).^2); % sum-of-squares of the residuals
    nrmse0    = 100 * sqrt(ssq_fit/length(coll_s0(:,1))) / mean(coll_s0(:,1));
    % =====================================================================
    
    % Do SPPE as in EFSA document. Expressed as percentage. This is the
    % deviation between observed and predicted survival *probability* at
    % the end of the test, for each treatment. The first two columns are
    % for the number of the data set and the number of the treatment.
    sppe = [coll_p_end(:,[3 4]) 100 * (coll_p_end(:,1) - coll_p_end(:,2))];
    
    % Finally, put them on-screen in a formatted manner.
    diary (['output_report',filesep,GLO.basenm,filesep,GLO.basenm,'.out']) % collect output in the diary *.out
    if flag_stage == 1
        disp('Goodness of fit measures for calibration')
    elseif flag_stage == 2
        disp('Goodness of fit measures for validation')
    end
    if sdit == 1
        disp('  Special case: stochastic death (SD)');
    else
        disp('  Special case: individual tolerance (IT)');
    end
    if data_all.nr_sets > 1 % add a note on screen when using multiple data sets
        disp('  Note: NSE and NRMSE are for all entered data sets together')
    end
    disp('=====================================================================')
    % =====================================================================
    % NEW after v.1.0: ALSO CALCULATE NSE AND NRMSE WITH t=0 IN THERE
    disp('NSE and NRMSE EXcluding results at t=0')
    fprintf('  Model efficiency (NSE, r-square)          : %#0.4f \n',r_square)
    fprintf('  Normalised root-means-square error (NRMSE): %#0.3f %% \n',nrmse)
    disp('NSE and NRMSE INcluding results at t=0')
    fprintf('  Model efficiency (NSE, r-square)          : %#0.4f \n',r_square0)
    fprintf('  Normalised root-means-square error (NRMSE): %#0.3f %% \n',nrmse0)
    % =====================================================================
    disp(' ')
    disp('Survival probability prediction error (SPPE) for each treatment')
    disp('     Data set     treatment    value')
    for ie = 1:size(sppe,1)
        fprintf('        %i      %10s     %+#0.3g %% \n',sppe(ie,1),scen_name{sppe(ie,1),sppe(ie,2)},sppe(ie,3))
    end
    disp('=====================================================================')
    disp('Warning: these measures need to be interpreted more qualitatively')
    disp('  as they, strictly speaking, do not apply to quantal data')
    disp(' ')
    diary off
    
end
