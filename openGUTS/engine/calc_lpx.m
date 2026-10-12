function [LPx,LPlo,LPhi] = calc_lpx(data_pred,flag,sdit)

% Usage: [LPx,LPlo,LPhi] = calc_lpx(data_pred,flag,sdit)
%
% Function to calculate LPx values, with confidence intervals. LPx is the
% 'exposure-multiplication factor': the factor by which the exposure
% profile needs to be multiplied to obtain a certain percentage mortality
% after a specified duration. Output of LPx is needed as this is passed to
% the plotting routine (to plot model results for profile times LPx).
% Output of LPlo and LPhi is used for batch calculations in
% <calc_batch_lpx> (which can be done with CIs, at least in the Matlab
% version).
%
% I make use of the fact that the scaled damage level is proportional to
% the exposure concentration. Multiplying an entire exposure profile is
% thus equivalent to the multiplication of the damage levels over time by
% the same factor. We therefore only need to calculate <Dw> once for every
% parameter set. If we ever need to implement the full GUTS model (with
% separate TK and damage) some shortcuts need to be reconsidered!
%
% Inputs
% <data_pred> data object containing the exposure scenario (e.g., FOCUS profile)
% <flag>      1) regular with CI on LPx, 2) no CI (very fast), 
%             3) batch processing (no output to screen, no diary, no CIs),
%             4) batch processing WITH CIs (slow)
% <sdit>      switch for SD (1) or IT (2)
% 
% Outputs
% <LPx>    matrix with LPx values (nr of data sets x nr of effect levels)
% <LPlo>   lower edge of the 95% CI; same format as <LPx>
% <LPhi>   upper edge of the 95% CI; same format as <LPx>
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

global GLO SETTINGS_PRED % make globals available

%% BLOCK 1. Initial things.

% Create text string for plotting titles and saving plots.
if sdit == 1
    savestr = [GLO.basenm,'_LPx_SD'];
else
    savestr = [GLO.basenm,'_LPx_IT'];
end

% Extract some settings from the globals defined in <initial_setup>.
F_LPx       = SETTINGS_PRED.F_LPx; % vector with x-values to use for LPx
S_range     = SETTINGS_PRED.S_range; % range for making a survival vs. multiplication plot
len_S_range = SETTINGS_PRED.len_S_range; % number of points for survival vs. multiplication plot
ft          = GLO.ft; % global for font specs in plots

tic % turn on the stopwatch to time this function

%% BLOCK 2. Load calibration file as saved by <calc_parspace>.

[pmat,sample] = load_sample(sdit,1); % load sample and parameter matrix from saved MAT file
pmat(GLO.hb,[1 5]) = [0 1]; % set background hazard to zero and normal scale!
% This is needed initially as we use the survival for the <pmat> for our
% initial guess of LPx. No need to set <hb> to zero in the sample as LPx is
% calculated without going through <calc_model> (well, it goes through only
% once to calculate damage, but not for the survival). 
if flag == 2 || flag == 3 % if we don't need CIs on LPx ...
    sample = []; % make the sample empty
else
    sample = sample(randperm(size(sample,1)),:); % randomly shuffle sets from outer hull
    % Note: this is for the Matlab version only, to get more appropriate
    % CIs when the user presses the <cancel> button in the wait bar.
end

% For batch calculations, it is possible to calculate the CIs with a
% smaller sample. Otherwise, calculations can be VERY time consuming! A new
% global setting GLO.n_lim is used which has the number of elements of the
% sub-sample. For now, this is ONLY possible for calculation of LPx with
% CIs in *batch mode*. In the future, we could consider using a limited set
% for plotting survival-at-LPx with CI as well. In principle, this could
% also be used more widely (then best moved to <load_sample> for other CIs
% as well (calibration, validation, LC50), though those calculations are
% not very time consuming.
if flag == 4 && isfield(GLO,'n_lim') && ~isempty(GLO.n_lim)
    n_lim = GLO.n_lim;
    if n_lim <= size(sample,1) % cannot take more than there are
        sample = sample(1:n_lim,:); % take <n_lim> sets from outer hull (since sample has been shuffled, this is randomised)
    end
end


%% BLOCK 3. Calculate LPx and MF-survival for the best-fitting parameter set.
% This block does two things: calculate MFs corresponding to a range of
% survival probabilities (to make a sort of dose-response curve), and to
% calculate the MF for specific values of survival (the LPx). Doing this
% separately is practical for SD where we cannot directly calculate MFs
% from a given survival probability.
% 
% Note: there can now be multiple data sets (in the Matlab version), but
% only one scenario per data set! For the software, we decided on one data
% set for predictions, but it is good to be prepared.

% BLOCK 3.1. Initial things. Initialise matrix to catch the LPx for the
% best-fitting parameter set, and the CI.

LPx  = zeros(data_pred.nr_sets,length(F_LPx)); % initialise LPx as zeros matrix (nr of data sets x nr of effect levels)
LPlo = nan(data_pred.nr_sets,length(F_LPx));   % initialise LPlo as matrix of NaNs
LPhi = nan(data_pred.nr_sets,length(F_LPx));   % initialise LPhi as matrix of NaNs

for i_d = 1:data_pred.nr_sets % run through all data sets
    
    % BLOCK 3.2. Initial preparation before each data set.
    
    ev_ic = data_pred.scen{i_d,1}; % only extract first scenario from this data set (there should be only 1)  
    t_end = ev_ic(end,1); % final time point for this scenario
    % Note: the final time is thus fixed to the last time point in the
    % exposure profile. If we ever decide to make this free (e.g., to use a
    % moving time window over the profile), this needs some modifications
    % in the code (especially below, in the calculation of <t_mod>).
    % Alternatively, for future time-window calculations, shortened <ev_ic>
    % matrices can be supplied to this function.
    
    coll_survmf = zeros(len_S_range,4); % initialise the matrix that will collect the survival-MF data, with its CI in column 3-4
    coll_survmf(:,3) = inf; % the third column will collect the minimum (for the CI), so already initialise with INF
    
    % =====================================================================
    % NEW after v.1.0: safer to always use a long time vector.
    % Create a time vector for the calculations as <t_mod>.
    % Make sure that <calc_model> returns an extra long time vector.
    t_mod = linspace(0,t_end,t_end*GLO.t_per_day); % minimum length of time vector for calculation (default each 15 minutes)
    % The same global default is used in <calc_model>.
    t_mod = t_mod(:); % make sure it is a column vector
    t_mod = unique([t_mod;ev_ic(:,1)]); % (vector)
    % this adds the time vector of the exposure scenario to the existing
    % vector <t_mod>. This is probably overkill and may slow down
    % calculations. However, especially for IT, it may give a bit of extra
    % accuracy.
    %     
    % NOTE: there is some redundancy in the code as <calc_model> will also
    % increase the time vector in the same manner. However, we here like to
    % have <Dw> with great precision (rather than survival at specific time
    % points), and <calc_model> only returns the time points that we ask it
    % to.
    % =====================================================================
    
    % BLOCK 3.3. Start by calculating damage (<Dw>) for this scenario, for
    % multiplication factor 1
    Xout = calc_model(t_mod,pmat,ev_ic,sdit); % survival and damage at scenario <ev_ic>
    Dw   = Xout(:,2); % scaled damage levels are in the second column
    
    % BLOCK 3.4. For IT we can use a direct calculation: for each <S>, we
    % can ask for the corresponding MF.
    if sdit == 2 % for IT
        
        maxDw = max(Dw); % find the maximum value for <Dw> over time
        % Create a linear range of survival probabilities between the
        % values that we aim for for the survival-MF plot. Place this range
        % in the second column of <coll_survmf>.
        coll_survmf(:,2) = (linspace(S_range(2),S_range(1),len_S_range))'; 
        for i_S = 1:len_S_range % run through this range
            [~,~,coll_survmf(i_S,1)] = calc_surv_it(maxDw,pmat,0,1-coll_survmf(i_S,2),[]); % calculate LPx in a sub-function
            % Note that x/100 = 1-S (at LP10, S=0.90). Each MF is placed in
            % the first column of <coll_survmf>.
        end
        % Next, calculate LPx for the effect levels that we need to report.
        for i_e = 1:length(F_LPx) % run through all effect levels
            [~,~,LPx(i_d,i_e)] = calc_surv_it(maxDw,pmat,0,F_LPx(i_e),[]); % calculation of LPx in the sub-function
        end
        % Note: both loops could have been avoided in the Matlab version as
        % the sub-function can be modified to handle input of a vector (and
        % returns a vector). However, using loops is probably more
        % transparent.
    
    % BLOCK 3.5. For SD, we need to work the other way around. When we have
    % an MF, we can calculate the corresponding <S>. Since we want to know
    % the MF for the LPx, we need to use <fzero>. The strategy to follow is
    % to first go roughly, in steps of 10x, through the MFs to find out
    % where the desired range (<S_range>) can be found. We will refine later.
    else % For SD, some more work is needed ... 
        
        % BLOCK 3.5.1. Initial things.
        mf_try = 1; % start with multiplication factor set to 1
        
        % No need to do the first one again, we already did that!
        S          = Xout(end,1); % survival at the end of the profile for best parameter set
        coll_rough = [mf_try S]; % start a two-column matrix with the rough factors (elements will be added to it as we go along)
        
        % BLOCK 3.5.2. Next, roughly find the range of multiplication
        % factors where the survival probabilities in <S_range> are (take
        % steps of 10x).
        while S > S_range(1) % go up until survival is very low (beyond first element of range)
            mf_try = mf_try * 10; % multiply last <mf_try> by factor of 10
            [~,S]  = calc_surv_sd_t(mf_try,pmat,t_mod,Dw,1); % calculate survival probability at end of profile for new MF
            coll_rough = [coll_rough ; mf_try S]; %#ok<AGROW> store new result in matrix, at bottom
        end
        mf_try = 1; % set MF back to 1
        S      = coll_rough(1,2); % and return survival to first entry
        while S < S_range(2) % go down until survival is very high (beyond last element of range)
            mf_try = mf_try / 10; % divide last <mf_try> by factor of 10
            [~,S]  = calc_surv_sd_t(mf_try,pmat,t_mod,Dw,1); % calculate survival probability at end of profile for new MF
            coll_rough = [mf_try S ; coll_rough]; %#ok<AGROW> store new result in matrix, on top
        end
                
        % BLOCK 3.5.3. Use root finding to zoom in on exact place of
        % survival probabilities in <S_range>, first small one (5% survival
        % by default).
        
        % First find the range in <coll_rough> in which the small survival point is located.
        index_hi = find(coll_rough(:,2)<S_range(1),1,'first'); % index to upper boundary for small survival
        index_lo = find(coll_rough(:,2)>S_range(1),1,'last');  % index to lower boundary for small survival
        mf_range = coll_rough([index_lo index_hi],1); % the small survival point is in this range of MFs
        mf_small = setup_rootfind(mf_range',2,pmat,t_mod,Dw,1-S_range(1)); % find the MF with root finding, from the reduced range
        
        % Next find the range in which the high survival point is located.
        % This follows the same principle.
        index_hi = find(coll_rough(:,2)<S_range(2),1,'first'); % index to upper boundary for high survival
        index_lo = find(coll_rough(:,2)>S_range(2),1,'last');  % index to lower boundary for high survival
        mf_range = coll_rough([index_lo index_hi],1); % the high survival point is in this range of MFs
        mf_large = setup_rootfind(mf_range',2,pmat,t_mod,Dw,1-S_range(2)); % find the MF with root finding, from the reduced range
        
        % BLOCK 3.5.4. Now, we can create a finer grid of multiplication
        % factors in the right spot.
        coll_survmf(:,1) = (logspace(log10(mf_large),log10(mf_small),len_S_range))'; 
        % Put the multiplication factors in the first column, log-spacing works better.
        
        for i_mf = 1:len_S_range % run through the MF values to try
            [~,S] = calc_surv_sd_t(coll_survmf(i_mf,1),pmat,t_mod,Dw,1); % calculate the associated survival probability
            coll_survmf(i_mf,2) = S; % collect it in the matrix in the second column
        end
       
        % BLOCK 3.5.5. Calculate the LPx values for the best parameter set,
        % for each specific effect level requested, using <fzero>.
        for i_e = 1:length(F_LPx) % run through effect levels requested
            % Use the initial range from the MF-survival matrix already
            % calculated to provide a narrow range for fzero.
            index_hi     = find(coll_survmf(:,2)<(1-F_LPx(i_e)),1,'first'); % index to upper boundary for this effect level
            index_lo     = find(coll_survmf(:,2)>(1-F_LPx(i_e)),1,'last');  % index to lower boundary for this effect level
            mf_range     = coll_survmf([index_lo index_hi],1); % the effect level is in this range of MFs
            LPx(i_d,i_e) = setup_rootfind(mf_range,2,pmat,t_mod,Dw,F_LPx(i_e)); % find the LPx with root finding, from the reduced range
            
        end
        
    end
   
    %% BLOCK 4. Output results for best-fitting parameter set.
    
    if flag < 3 % don't do this when we are doing batch processing
        
        % BLOCK 4.1. Print LPx to screen.
        if flag == 2 % print output to screen and diary
            % This is now only done for cases where we want no CI. If we
            % want CIs, we'll have to wait until the calculations for the
            % CI are finished before we see the LPx. Depending on
            % calculation speed, we may want to reconsider this.

            diary (['output_report',filesep,GLO.basenm,filesep,GLO.basenm,'.out']) % collect output in the diary *.out
            disp(' ')
            disp(['Results table for LPx (-) with 95% CI: ',data_pred.name{i_d}])
            if sdit == 1
                disp('  Special case: stochastic death (SD)');
            else
                disp('  Special case: individual tolerance (IT)');
            end
            disp('==================================================================')
            for i_e = 1:length(F_LPx) % run through effect levels requested
                fprintf('LP%1.0f at end of profile: %#5.4g \n',F_LPx(i_e)*100,LPx(i_d,i_e))
            end
            disp('==================================================================')
            diary off
        end
        
        % BLOCK 4.2. make a plot of survival vs. MF. If CIs are requested,
        % they will be plotted in the same figure later on.
        figh = make_fig(1,1,2); % create an empty figure window of correct size
        hold on % hold it even when more things are plotted
        h_ax = gca; % handle for current axis system
        set(h_ax,'LineWidth',1,ft.name,ft.ticks) % adapt axis formatting to standard settings
        h_tmp = title(['set: ',num2str(i_d),', ',data_pred.name{i_d}],'Interpreter','none'); set(h_tmp,ft.name,ft.title) % add a title with name of file with exposure profile
        
        h_tmp = ylabel('survival probability (-)'); set(h_tmp,ft.name,ft.label) % add label to y-axis
        h_tmp = xlabel('multiplication factor LPx (-)'); set(h_tmp,ft.name,ft.label) % add label to x-axis
        
        h_plot = plot(coll_survmf(:,1),coll_survmf(:,2),'k-','LineWidth',1.5); % plot best line for survival vs. MF
        % set(h_ax, 'XScale', 'log') % set x-axis to log-scale, which might work better for these plots
        
        if flag > 1 % then we don't need CIs, so make an overall title for the plot with the filename, and save it
            axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off','Visible','off','Units','normalized', 'clipping','off');
            h_txt = text(0.5, 1,['File: ',savestr],'Interpreter','none','HorizontalAlignment','center','VerticalAlignment', 'top');
            set(h_txt,ft.name,ft.text);
            drawnow
            save_plot(figh,['output_report',filesep,GLO.basenm,filesep,savestr],GLO.save_format) % save lcx plot in output folder
            snapnow % for publishing Matlab scripts: make plot appear in the context of the screen output
        end
        drawnow % empty plotting buffer
    end
    
    %% BLOCK 5. Calculate confidence intervals.
    % Use the saved sample from the parameters-space explorer, and take min
    % and max of results. Note: if we don't need CIs, we'll still run
    % through this code. However, the sample is made empty, so the <for>
    % loop won't even start.
    
    % BLOCK 5.1. Initial things
    n_sets  = size(sample,1); % number of sets in the sample
    ind_log = pmat(:,5)==0;   % index to parameters on log scale (vector of logical indices)
        
    % Use the fancy progress bar from Matlab
    if n_sets > 0 && flag < 3 % only if there are sets in the sample, and not for batch calculations
        f = waitbar(0,sprintf('Bar shows fraction of %1.0f sets calculated',n_sets),'Name','Calculating CIs on LPx','CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
    else
        f = []; % no waitbar, so make the handle empty!
    end
    
    pmat_tmp = pmat; % make a copy of <pmat> as we will change it for each set in the sample
    % This is needed if we have multiple data sets to run through.
    
    for i_n = 1:n_sets % run through all sets in the sample
        
        % Progress bar from Matlab: check if cancel is pressed
        if ~isempty(f)
            if getappdata(f,'canceling')
                break % break the <for> loop
            end
            waitbar(i_n/n_sets,f) % update waitbar
        end
        
        % BLOCK 5.2. Initial things for each set from the sample.
        pmat_tmp(:,1) = sample(i_n,:); % replace values in first column of <pmat> with the next set in the sample
        
        % Put parameters that need to be on log scale back on normal scale
        % (as the SD/IT calculations require normal scale).
        if sum(ind_log)>0
            pmat_tmp(ind_log,1) = 10.^(pmat_tmp(ind_log,1));
        end
        % Note: <pmat> is on normal scale here, but the sample contains the
        % value on a log scale, if a parameter is fitted on log scale.
        
        Xout = calc_model(t_mod,pmat_tmp,ev_ic,sdit); % survival and damage at scenario <ev_ic>
        Dw   = Xout(:,2); % scaled damage levels for this parameter set
        
        % BLOCK 5.3. For IT we can use a direct calculation of <S>.
        if sdit == 2 % for IT ... 
            
            maxDw = max(Dw); % find the maximum value for <Dw> over time
            
            % First the range of <S>-values for the survival-MF plot.
            for i_mf = 1:len_S_range % run through the MF values to try 
                S = calc_surv_it(maxDw,pmat_tmp,0,[],coll_survmf(i_mf,1)); % calculate survival probability in sub-function
                coll_survmf(i_mf,3) = min([coll_survmf(i_mf,3) S]); % update the matrix in the third column (min CI)
                coll_survmf(i_mf,4) = max([coll_survmf(i_mf,4) S]); % update the matrix in the fourth column (max CI)
            end
            
            % Next, the specific LPx values that we need to report.
            for i_e = 1:length(F_LPx) % run through all effect levels
                [~,~,LPtmp] = calc_surv_it(maxDw,pmat_tmp,0,F_LPx(i_e),[]); % calculation of LPx in sub-function
                LPlo(i_d,i_e) = min([LPlo(i_d,i_e),LPtmp]); % keep the lowest value as minimum
                LPhi(i_d,i_e) = max([LPhi(i_d,i_e),LPtmp]); % keep the highest value as maximum
            end
            
        % BLOCK 5.4. For SD, things are more complex (as it was for the
        % best parameter set). We won't do a rough search again, but only
        % use the range of MF values that we already established.
        else % for SD, some more work is needed ... need to use <fzero>
                     
            % BLOCK 5.4.1. Use the finer grid of multiplication factors
            % from the best fit again. Calculate the corresponding survival
            % for each set from the sample and keep highest/lowest.
            coll_c = zeros(length(mf_try),1); % initialise a vector to collect the <S> at each MF
            for i_mf = 1:len_S_range % run through the MF values to try
                [~,S] = calc_surv_sd_t(coll_survmf(i_mf,1),pmat_tmp,t_mod,Dw,1); % calculate survival probability for this MF
                coll_survmf(i_mf,3) = min([coll_survmf(i_mf,3) S]); % update the min in the matrix in the third column
                coll_survmf(i_mf,4) = max([coll_survmf(i_mf,4) S]); % update the max in the matrix in the fourth column
                coll_c(i_mf) = S; % also collect the survival itself for zooming in on specific LPx values below
            end
            
            % BLOCK 5.4.2. Calculate the LPx values for this set, for each
            % effect level requested. Make use of the MF-S combinations
            % already calculated (<S> in <coll_c>, MF in first column of
            % <coll_survmf>).
            for i_e = 1:length(F_LPx) % run through all effect levels
                % Extract an initial range from the MF-survival matrix already calculated.
                index_hi = find(coll_c<(1-F_LPx(i_e)),1,'first'); % index to upper boundary for this effect level
                index_lo = find(coll_c>(1-F_LPx(i_e)),1,'last');  % index to lower boundary for this effect level
                
                % NEW: There can be cases when <coll_c> does not go low
                % enough to catch an LP10 for a set in the sample (i.e.,
                % that the first entry in <coll_c> already has more than
                % 10% effect). That is no problem for <fzero>, however, it
                % is for the algorithm used in the C++ version. Therefore,
                % I try to find new bounds by changing the last one by a
                % factor of 10.
                if ~isempty(index_lo) && ~isempty(index_hi) % then all is okay
                    mf_range = coll_survmf([index_lo index_hi],1); % the effect level is in this range of MFs in first column
                end
                if isempty(index_lo) % then no lower bound for MF was found
                    mf_range = [coll_survmf(1,1) coll_survmf(index_hi,1)]; % start from first entry for lower bound
                    crit = -1;
                    while crit < 0
                        mf_range(1) = mf_range(1) / 10; % divide first entry by 10 till it crosses the threshold
                        crit = calc_surv_sd_t(mf_range(1),pmat_tmp,t_mod,Dw,F_LPx(i_e)); % calculate new criterion for first entry in <mf_range>
                    end
                end
                if isempty(index_hi) % then no upper bound for MF was found (which would be strange, but better safe than sorry)
                    if isempty(index_lo) % then they were both empty ... which would be really, really silly
                        mf_range(2) = coll_survmf(end,1); % then keep the first entry (which was filled above) and start from last value in <coll_surfmf>
                    else
                        mf_range = [coll_survmf(index_lo,1) coll_survmf(end,1)]; % start from last value in <coll_surfmf> for upper range
                    end
                    crit = 1;
                    while crit > 0
                        mf_range(2) = mf_range(2) * 10; % multiply second entry by 10 till it crosses the threshold
                        crit = calc_surv_sd_t(mf_range(2),pmat_tmp,t_mod,Dw,F_LPx(i_e)); % calculate new criterion for second entry in <mf_range>
                    end
                end    
                                
                LPx_tmp = setup_rootfind(mf_range,2,pmat_tmp,t_mod,Dw,F_LPx(i_e)); % find the LPx with root finding, from the reduced range
                LPlo(i_d,i_e) = min([LPlo(i_d,i_e),LPx_tmp]); % keep lowest value as minimum
                LPhi(i_d,i_e) = max([LPhi(i_d,i_e),LPx_tmp]); % keep highest value as maximum
            end
        end
        
    end
    
    if ~isempty(f)
        delete(f) % delete progress bar again
    end
    
    %% BLOCK 6. Print results to screen and finish the plot with CIs.
    
    if flag == 1 % only for regular analyses with CI
        % for <flag=2>, the results without CIs have already been printed,
        % and for <flag=4> we don't do output to screen here (only in
        % <calc_batch_lpx>).
        
        % BLOCK 6.1. Print to screen.
        diary (['output_report',filesep,GLO.basenm,filesep,GLO.basenm,'.out']) % collect output in the diary *.out
        disp(['Results table for LPx (-) with 95% CI: ',data_pred.name{i_d}])
        if sdit == 1
            disp('  Special case: stochastic death (SD)');
        else
            disp('  Special case: individual tolerance (IT)');
        end
        if i_n < n_sets
            fprintf('  Calculation of CIs stopped by user after %1.0f of %1.0f sets \n',i_n,n_sets)
            disp('  (note that the sample was randomised before LPx calculation)')
        end
        disp('==================================================================')
        fprintf('LPx     best       CI')
        fprintf('\n')
        disp('==================================================================')
        for i_e = 1:length(F_LPx) % run through effect levels
            fprintf('LP%1.0f: %#5.4g (%#5.4g - %#5.4g)   ',F_LPx(i_e)*100,LPx(i_d,i_e),LPlo(i_d,i_e),LPhi(i_d,i_e));
            fprintf('\n');
        end
        disp('==================================================================')
        disp(['Time required for the LPx calculations: ' secs2hms(toc)])
        disp(' ')
        diary off
        
        % BLOCK 6.2. Finish figure by plotting CIs in there.
        figure(figh) % make the figure we already have current, just to make sure we are plotting in the right place
        
        plot(coll_survmf(:,1),coll_survmf(:,3),'k:','LineWidth',1) % plot each confidence line as broken line
        plot(coll_survmf(:,1),coll_survmf(:,4),'k:','LineWidth',1) % plot each confidence line as broken line
        
        % little trick to fill the area between the two curves, to obtain a
        % coloured confidence band
        t2  = [coll_survmf(:,1);flipud(coll_survmf(:,1))]; % make a new time vector that is old one, plus flipped one
        Xin = [coll_survmf(:,3);flipud(coll_survmf(:,4))]; % do the same for the plot line, hi and lo
        fill(t2,Xin,'g','LineStyle','none','FaceAlpha',0.25) % and fill this object
        
        set(h_ax, 'Layer', 'top') % make sure the axes are visible (the yellow object wants to obscure it)
        uistack(h_plot,'top');    % place best model line on top again, using the plot handles
        
        % Make an overall title for the plot with the filename, and save it
        axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off','Visible','off','Units','normalized', 'clipping','off');
        h_txt = text(0.5, 1,['File: ',savestr],'Interpreter','none','HorizontalAlignment','center','VerticalAlignment', 'top');
        set(h_txt,ft.name,ft.text);
        drawnow
        save_plot(figh,['output_report',filesep,GLO.basenm,filesep,savestr],GLO.save_format) % save MF-survival plot in output folder
        snapnow % for publishing Matlab scripts: make plot appear in the context of the screen output
        
    end
end

