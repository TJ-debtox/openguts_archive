function [phat,fval,exitflag] = setup_simplex(pfit,option,varargin)

% Usage: [phat,fval,exitflag] = setup_simplex(pfit,option,varargin)
% 
% A wrapper function around the optimisation routine. This allows to change
% the optimisation routine without needing to change all the calls in the
% rest of the software. There is currently a choice between Matlab's
% <fminsearch> and the <nelmin> that is used for the C++ code.
% 
% Inputs
% <pfit>        vector with parameters to be fitted (initial values)
% <option>      1) standard optimisation, 0) rough optimisation
% <varargin>    additional arguments to be passed on to <transfer>
% 
% Outputs
% <phat>        vector with optimised parameter values
% <fval>        optimised value of minus-log-likelihood
% <exitflag>    1) convergence, 0) no convergence or errors
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

use_optim = 2; % 1) Matlab's <fminsearch>, 2) <nelmin> from web

% Here, we need to catch cases where optimisation is requested but no
% parameters are fitted. This happens for example when we fit only one
% parameter and then calculate a profile likelihood. For <fminsearch> this
% is not an issue, but <nelmin> will return <fval> as empty which causes
% problems.
if isempty(pfit) % then we don't need to fit anything
    phat     = []; % so <phat> is empty as well
    fval     = transfer(pfit,varargin{:}); % calculate MLL
    exitflag = 1; % signal that convergence is achieved
    return % and return to where we came from
end
    
if use_optim == 1 % use <fminsearch>
    
    % Make options for the standard and rough simplex optimisations.
    oldopts = optimset('fminsearch'); % retrieve default options for <fminsearch>
    opts    = optimset(oldopts,'FunValCheck','on','Display','off'); % function check on, show no output at all
    if option == 0 % options for rough simplex
        opts = optimset(opts,'TolX',1e-2,'TolFun',1e-2,'MaxFunEvals',30*length(pfit)); % rough Simplex options, based on number of fitted parameters
        % Note: these rough options will use a smaller total number of function
        % evaluations: quick and dirty.
    end
    
    % Use <fminsearch> and the standard <transfer> function to fit.
    [phat,fval,exitflag] = fminsearch('transfer',pfit,opts,varargin{:}); % do a normal optimisation
    
elseif use_optim == 2 % use <nelmin>
    
    n      = length(pfit); % number of fitted parameters
    start  = pfit; % starting values for the fit parameters (vector)
    step   = 0.05*pfit; % intitial simplex; as in Matlab 5% increase
    % If parameters can be zero, this needs to be adapted. Parameters will
    % never be zero for the GUTS software at this moment, but it is good to
    % add a failsafe just in case.
    step(step==0) = 0.00025; % this value is also used in <fminsearch>
    
    konvge = 1; % check for convergence every <konvge> function evaluations
    if option == 0 % options for rough simplex
        kcount = 30*length(pfit); % maximum number of function evaluations
        reqmin = 1e-8; % minimum variance of the simplex point for stopping
    else % options for standard simplex
        kcount = 200*length(pfit); % maximum number of function evaluations
        reqmin = 1e-15; % minimum variance of the simplex point for stopping
        % The min variance cannot be made much lower as that triggers a
        % large number of superfluous restarts (for some reason). The rough
        % options don't have this issue as restarting is disabled in
        % <nelmin> for small max numbers of evaluations.
    end
    
    % Now we can call <nelmin> with its original function call format.
    [phat,fval,~,~,ifault] = nelmin(@transfer,n,start,reqmin,step,konvge,kcount,varargin{:});
    % NOTE: third and fourth output can be icount and numres that are not
    % used here.
    
    phat = phat(:); % make sure it is a column vector (like output of <fminsearch>)
    
    if ifault == 0 % no faults found?
        exitflag = 1; % signal that convergence is achieved
    else
        exitflag = 0; % no convergence
    end
    
end