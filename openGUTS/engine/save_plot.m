function save_plot(figh,savenm,flag)

% Usage: save_plot(figh,savenm,flag)
%
% For saving figure windows in various styles. This is for the Matlab
% version only.
% 
% As inputs the figure handle <figh>, the name to save the plot with in 
% <savenm>. The <flag> is used to signal the file type: 
% 1) save as Matlab FIG
% 2) JPEG
% 3) PDF
% 4) PNG
% 5) Windows meta file
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

switch flag
    case 1 % save it as a fig file
        saveas(figh,savenm)
    case 2 % save it as a jpeg file
        print(figh,'-djpeg',savenm)
    case 3 % save as PDF; some trickery to remove extraneous white space
        warning off
        figh.PaperPositionMode = 'auto';
        fig_pos = figh.PaperPosition;
        figh.PaperSize = [fig_pos(3) fig_pos(4)];
        print(figh,savenm,'-dpdf')
        warning on
    case 4 % save as PNG
        print(figh,'-dpng',savenm)
    case 5 % save as Windows meta file
        print(figh,'-dmeta',savenm)
end
