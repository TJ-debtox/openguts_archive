function root_value = setup_rootfind(prange,option,varargin)

% Usage: root_value = setup_rootfind(prange,option,varargin)
% 
% A wrapper function around the zero-finding routine. This allows to change
% the zero-finding routine without needing to change all the calls in the
% rest of the software. There is currently a choice between Matlab's
% <fzero> and the <zero> that is used for the C++ code.
% 
% Inputs
% <prange>      vector with parameter range to be searched
% <option>      1) use for LCx, 2) use for LPx
% <varargin>    additional arguments to be passed on to <transfer>
% 
% Outputs
% <root_value>  value of the parameter for which the criterion is zero
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

use_optim = 2; % 1) Matlab's <fzero>, 2) <zero> from web

if use_optim == 1 % use <fzero>
    
    % Use <fzero> function to find the root (standard options).
    % Note: the empty matrix is the place for optional options for <fzero>.
    if option == 1 % need to calculate an LCx, so use function <calc_surv_sd_c>
        root_value = fzero(@calc_surv_sd_c,prange,[],varargin{:}); 
    elseif option == 2 % need to calculate an LPx, so use function <calc_surv_sd_t>
        root_value = fzero(@calc_surv_sd_t,prange,[],varargin{:}); 
    end
    
elseif use_optim == 2 % use <zero> from web
    
    tol = 0; % I don't think any tolerance above machine precision is needed (setting <tol> to zero gives same precision as Matlab)
    
    % Note: <eps> is the machine precision in Matlab.
    if option == 1 % need to calculate an LCx, so use function <calc_surv_sd_c>
        root_value = zero(prange(1),prange(2),eps,tol,@calc_surv_sd_c,varargin{:});
    elseif option == 2 % need to calculate an LPx, so use function <calc_surv_sd_t>
        root_value = zero(prange(1),prange(2),eps,tol,@calc_surv_sd_t,varargin{:});
    end
    
end