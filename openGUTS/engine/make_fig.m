function figh = make_fig(m,n,flag)

% Usage: figh = make_fig(m,n,flag)
%
% A little function that creates an empty figure window of a certain
% default size, with multiple panels (<m> by <n>). This is for the Matlab
% version only.
% 
% This function only works when the setting is that figures are NOT docked!
% Figures of size 1x1, 1x2, 1x3, 2x2 and 2x3 will be plotted such that the
% individual plots are (or at least: should be) equal sized. Larger plots
% receive a large standard size that should be okay in most cases.
%
% For the standard plots in <main_plot>, an additional window size is used
% (<flag> = 1). The height <m> is always the same (as we have three rows),
% and the width is based on the number of treatments (in <n>) up to a
% maximum (otherwise it won't fit on the screen). The figure handle is
% returned by the function in <figh>.
%
% The position setting is [left bottom width height]. You might want to
% adapt these settings to your monitor. Left is the number of pixels from
% the left edge of the screen and bottom from the bottom edge.
%
% <flag>:  0) invisible standard plot, 1) visible standard plot, 2) visible
%          custom plot.
%
% Author     : Tjalling Jager
% Date       : May 2020
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

if flag == 0
    figh = figure('visible','off'); % new INVISIBLE figure with default size and position
else
    figh = figure; % new figure with default size and position
end

if strcmp(get(0,'DefaultFigureWindowStyle'),'docked') % this should not be needed as <initial_setup> makes sure it is undocked
    return % when figures are docked, no need to resize (won't work anyway), so simply return
end

scrsz = get(groot,'ScreenSize'); % get screen size that user is using
% =========================================================================
% % NEW after v.1.0: a little trick to make sure that normal figures are
% made in batch mode, instead of smaller ones. When Matlab is run without a
% display, the screensize is set very small. This still does does not
% produce plots that are exactly the same as on a regular screen, for some
% reason ...
if usejava('desktop') == 0 % this setting is zero when run in batch mode (nodisplay)
    scrsz = [1           1        1680        1050]; % just a regular large screen
end 
% =========================================================================

if flag ~= 2 % then we have the standard plot with three rows
    
    figsz = [min(n*300,1500) 770]; % set a default figure size
    extrsz = [100 200]; % pixels from top-right screen corner, for top-right corner of figure window
    
    % On small screens, Matlab does not allow certain figure sizes, so just make them small enough ...
    if scrsz(4) < 870 % this is probably the limit ...
        figsz  = floor(figsz/1.6); % decrease size
        extrsz = [50 100]; % modify position for small screens
    end
    figh.Position = [scrsz(3)-figsz(1)-extrsz(1) scrsz(4)-figsz(2)-extrsz(2) figsz(1) figsz(2)]; % and give it the new size
       
else % we have another type of plot, size is set based on the panel structure (<m> by <n>)
    
    % Figure size for 1x1, 1x2, 1x3, 2x2, 2x3, larger square, larger rectangular
    figsz  = [420 330;932 330;1455 330;935 749;1455 749;1000 770;1500 770];
    extrsz = [100 200]; % pixels from top-right screen corner, for top-right corner of figure window
    
    % On small screens, Matlab does not allow certain figure sizes, so just make them small enough ...
    if scrsz(4) < 870 % this is probably the limit ...
        figsz(3:end,:) = floor(figsz(3:end,:)/1.6); % only decrease the ones with more than 1 row, and the 1x3!
        extrsz = [50 100]; % modify position for small screens
    end

    if m == 1
        if n == 1
            figh.Position = [scrsz(3)-figsz(1,1)-extrsz(1) scrsz(4)-figsz(1,2)-extrsz(2) figsz(1,1) figsz(1,2)]; % and give them a default size for 1x1
        elseif n == 2
            figh.Position = [scrsz(3)-figsz(2,1)-extrsz(1) scrsz(4)-figsz(2,2)-extrsz(2) figsz(2,1) figsz(2,2)]; % and give them a default size for 1x2
        elseif n == 3
            figh.Position = [scrsz(3)-figsz(3,1)-extrsz(1) scrsz(4)-figsz(3,2)-extrsz(2) figsz(3,1) figsz(3,2)]; % and give them a default size for 1x3
        end
    elseif m == 2
        if n == 2
            figh.Position = [scrsz(3)-figsz(4,1)-extrsz(1) scrsz(4)-figsz(4,2)-extrsz(2) figsz(4,1) figsz(4,2)]; % and give them a default size for 2x2
        elseif n ==3
            figh.Position = [scrsz(3)-figsz(5,1)-extrsz(1) scrsz(4)-figsz(5,2)-extrsz(2) figsz(5,1) figsz(5,2)]; % and give them a default size for 2x3
            % this may need modification as it is determined using the BYOM sensitivity plots which are modified in size
        end
    else % it is higher than two sub-plots, so pretty big
        if m == n % if it is square
            figh.Position = [scrsz(3)-figsz(6,1)-extrsz(1) scrsz(4)-figsz(6,2)-extrsz(2) figsz(6,1) figsz(6,2)]; % and give them a default size for 3x3 and larger
        else % if it is rectangular
            figh.Position = [scrsz(3)-figsz(7,1)-extrsz(1) scrsz(4)-figsz(7,2)-extrsz(2) figsz(7,1) figsz(7,2)]; % and give them a default size for 3x4 and larger
        end
    end
    
end

