# openguts_archive
This repository is the archive for the openGUTS Matlab version. This archives version 1.2 of 20 March 2024.

## ⚠️ Project Status: Archived

**This repository has been archived following the closure of business operations at DEBtox Research and is no longer actively maintained. The code is provided strictly "as-is".**

This repository only concerns version 1.2 of openGUTS. I expect there to be updates of openGUTS in the future, but these will not be offered or archived here (under TJ-debtox). This archive will likely be removed when a new version is released; the openGUTS website ([https://www.openguts.info/](https://www.openguts.info/)) will always link to the most recent version.

### Bug Reports & Security Warnings
Since this repository is locked (read-only), new issues or pull requests cannot be opened. However, if a critical bug or error is discovered in these MATLAB files in the future, a warning notice will be posted here:

* **Future updates/warnings will be shared via:** [https://www.openguts.info/download.html](https://www.openguts.info/download.html)
* You can reach me at: **tjalling [at] debtox [dot] nl** (which I will try to keep active) or [via LinkedIn](https://nl.linkedin.com/pub/tjalling-jager/5/906/147)


## Licensing

The openGUTS project as a whole is licensed under the **GNU General Public License v3.0** (GPLv3) — see the [LICENSE](LICENSE) file in the root directory for the full terms. 

### Author

The Matlab code for openGUTS versions 1-1.2 has been written by Tjalling Jager (DEBtox Research). Exceptions are the third-party functions listed below, where the original authors are identified in the file header.

### Funding

Version 1.0 of the software was funded by Cefic-LRI through the project ECO39.2. Version 1.1 and 1.2 were funded by Syngenta Crop Protection AG. 

### Third-party licenses

Please note that that the third-party engine functions nelmin.m and zero.m are licensed under the GNU LGPL license, and the function sec2hms.m under an MIT license. Each third-party file retains its original author's copyright headers including specific license information. Modifications to the code are noted in the comments of these files.

## What is openGUTS?

OpenGUTS is the implementation of the General Unified Threshold model for Survival (GUTS) into software. GUTS is a toxicokinetic-toxicodynamic (TKTD) model for the endpoint survival. For more information, see the supporting web site (see link below).

OpenGUTS is available as a freely-downloadable, open source, standalone for Windows, and as a set of Matlab scripts and functions. This repository is for the Matlab version. The separate repository for the standalone version contains the Windows installer and C++ source code, and can be found at [https://github.com/TJ-debtox/openguts_standalone/](https://github.com/TJ-debtox/openguts_standalone/).

The Matlab version has the same functionality as the standalone version, but lacks the user interface and the option to make formatted output reports. Both versions are able to use the same input files for the data sets and the exposure profiles. The Matlab version will be handy as a platform for modification and experimentation. Furthermore, since the Matlab version was originally developed as prototype for the C++ source code, both follow the same structure, and the Matlab version can be used for future development of the software. 

### About this archive

As of 1 December 2026, I (Tjalling Jager) have retired my business DEBtox Research. I also maintain the [website for openGUTS](https://www.openguts.info/). The website remains in the air, but the software was moved to GitHub. This archive on GitHub contains the last release that I have contributed to: version 1.2 of 20 March 2024. Future updates will be hosted elsewhere (see supporting website for the latest version and the correct download link).

### Supporting web site

The supporting web site for openGUTS is [https://www.openguts.info/](https://www.openguts.info/). This website no longer contains the downloads for the Matlab files, but still hosts manuals, code 'walk-throughs', description of openGUTS, etc. This page will also collect bugs, errors, and other issues that may arise after archiving the repository on GitHub.

## Installation and folder structure

To install the openGUTS Matlab package, unpack the downloaded ZIP file to a location of your choice. Do NOT rename the directories `/openGUTS/` and `/openGUTS/engine/`. You can move the complete folder `/openGUTS/` (with all its subfolders) out of the main folder `/openguts_archive-main/` if you like. Do NOT set the path manually to any of the openGUTS folders, and please make sure that the directory where your main script is located is the active directory (if you run a script when it is not in the active directory please select <Change folder> and NOT <Add to Path>).

The reason for not changing the folder names and modifying the Matlab path is that this is done automatically in `pathdefine.m`, which should be present in any folder where you run the script file `openguts.m` or `openguts_batch.m`. The folder names are hard-coded in `pathdefine.m`. This function makes sure that the path is set correctly. It will also remove redundant path settings from previous runs. The changes to the path are temporary. Also note that the `pathdefine.m` of openGUTS will remove any path settings made by the BYOM Matlab version. This allows both platforms to be used in the same Matlab session without path conflicts.

The `/engine/` directory contains all the Matlab functions for fitting the GUTS model to data, plotting results, statistical analysis etc. There should be no need to modify these files. 

## Getting started

The file to start with is `openguts.m`, run it, and select an input data set, e.g., `/input_data/propiconazole_constant.txt` and watch the analysis proceed. You can find a detailed manual on the download page: [https://www.openguts.info/download.html](https://www.openguts.info/download.html) that also explains how to create and run your own data sets.
