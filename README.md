<!-- ######################################################################## -->

***
# MTWSPy
***

MTWSPy is a Python implementation of the Morphological Time Window Selection (MTWS) method first published in _Geophysical Journal International_ by Lei Li et al., ([2023](https://doi.org/10.1093/gji/ggad338)), including preceeding data processing workflow.

<!-- ######################################################################## -->

***
## Description
***

This package includes codes to download and process seismic data, create corresponding synthetics and measure travel times using MTWSPy and cross correlation.

* We use [iris-fetch-scripts](https://github.com/EarthScope/fetch-scripts?tab=readme-ov-file) for download of seismic data.
* We use [SPECFEM3D_GLOBE](https://specfem.org/) v8.0.0 available [here](https://github.com/SPECFEM/specfem3d_globe/releases/tag/v8.0.0) for calculation of global synthetic seismograms.
* We use [ObsPy](https://docs.obspy.org/) for observed and synthetic data processing.

The codes are designed to be run on a year-by-year basis to build up a travel time database progressively, thereby avoiding data storage issues where possible.

HPC operations are setup using an example from the PSMN cluster at ENS Lyon, France (96 cores per node) using the SLURM batch scheduler.

Parallelisation is used to speed up post-processing of observed and synthetic data and arrival time picking (python concurrent.futures) using the MTWS algorithm.


<!-- ######################################################################## -->

***
## Installation
***

### Data & Synthetics

1. Ensure SPECFEM3D_GLOBE can be compiled and run using parameters appropriate to the system used to compute synthetics (typically HPC). v8.0.0 installed fresh from Github.
2. Ensure the following files are available on your PATH e.g., in a local bin directory : ~/bin/:
    * FetchMetadata (iris link above)
    * FetchEvent (iris link above)
    * FetchData (iris link above)
    * parallel_proc_specfem_seis.py (MTWSPy/utils/specfem)

### MTWSPy

`conda env create -f environment.yaml`

`pip install -r requirements.txt`

`pip install (-e) . for --editable`

`May need to: export PYTHONPATH=./MTWSPy:$PYTHONPATH`

`conda activate venv_MTWSPy`

`pytest` (should pass all tests)

`MTWSPy_main` (should run code on test data provided)


<!-- ######################################################################## -->

***
## Usage
***

### Data Download

navigate to `utils/EarthScope_FetchData`

* `FetchData.py` &rarr; launches obspy download via batch scheduler for given year (using run_FetchData.bash).

### Generate Synthetics

navigate to `utils/gcmt`

* `python download_cmtsolutions.py` &rarr; Download from gcmt catalog CMT solutions in Specfem format for 0-700km depth, 5.5-7.5Mw for given year.

navigate to `utils/specfem`

* `python mk_SPECFEM_STATIONS_obspy.py` &rarr; Creates station file for Specfem using datacenter quieries (often more frequently used).
* Par file  &rarr; Verify parameters chosen are appropirate for your system (e.g., minimum period). Normally using 2hrs seismogram length.
* Station file &rarr; The STATION file will be updated from the process described above.
* `get_specfem_synthetics.sh` &rarr; Launches all SPECFEM simulations for CMT solutions in gcmt directory using batch scheduler (and run_specfem_slurm_TEMPLATE.bash). May need modification to suit local resources.
* `pyproc_specfem_synthetics.sh` &rarr; Launches post processing of SPECFEM simulation data using python parallel (and run_pyproc_specfem_TEMPLATE.bash) and batch scheduler.

### MTWSPy _Main_

navigate to `home`

#### Contents:

* __params_in.yaml__ &rarr; Parameter file for main code (must be consulted prior to running)
* __MTWSPy/v01_phasenames.py__ &rarr; Dictionary of phases desired
* __MTWSPy/toolkit.py__ &rarr; Global functions
* __MTWSPy/match_catalog.py__ &rarr; Make useable csv file of possible events
* __MTWSPy/find_twin.py__ &rarr; Find twins for observed & synth data
* __MTWSPy/match_twin_files.py__ &rarr; Find matching twin files obs/synth 
* __MTWSPy/phase_association.py__ &rarr; Associate observed & synth twin with a phase
* __MTWSPy/correlate_twin.py__ &rarr; Correlate observed and synth time windows
* __MTWSPy/MTWSPy_main.py__ &rarr; Execute main code
* __MTWSPy/run_MTWSPy.bash__ &rarr; Slurm batch scheduler 
* __MTWSPy/post_processing/process_tdl_files.py__  &rarr; Time delay processing code
* __MTWSPy/post_processing/compare_picks.py__  &rarr; Time delay comparison code
* __MTWSPy/post_processing/plot_phase.gmt__  &rarr; Time delay plotting code


#### Description:

The overall philospy of the code is to have a setup where we can just hit `GO` or `MTWSPY_main from the directory where the code is installed.

To do this we rely heavily on an input parameter file `params_in.yaml` which must be studied and adapted for purpose before using the code. All parameters are described/commented and are grouped with each relevant part of the code. This file is loaded before the execution of any part of the code.

-> open params_in.yaml
-> Check home, data_loc, synth_loc, fmt_data_loc, cmt_infile, year, parallel, cores

Where possible the code also relies on a toolkit of global functions in `toolkit.py` these should not be changed but are crucial to the effective running of the code (e.g., initialising parallelisation of commands).

Although each step of the code can be launched individually (using the executable created by conda at `XXX.py`) recommended usage is carried out by executing the full code: `MTWSPy_main` 

All steps of the code that utilise parallelisation (also written to work in serial) have a common input format:

* __main_function__ &rarr; Usually `process_one_event`, checks/prepares specific inputs for each function(_not to be changed_)
* __input_directory__ &rarr; Location of input data (_set based on parameter file_)
* __evt_id_tab__ &rarr; Loaded events (_pandas dataframe_)
* __functions__ &rarr; List of specifc functions for given code step (_unlikely to be changed_)
* __params_in__ &rarr; Loaded parameter file (_dict_)
* __phases__ &rarr; Loaded phase file (_dict_)

#### Detailed description:

__params_in.yaml__ &rarr; Parameter file containing all variables that can be changed within the main code.

__MTWSPy/v01_phasenames.py__ &rarr; Dictionary of phases desired, could be modified to search for alternative seismic phases

__MTWSPy/toolkit.py__ &rarr; Global function definitions used throughout the code to speed up execution.

__MTWSPy/match_catalog.py__ &rarr; Make useable csv file of possible earthquakes starting with gcmt catalog file at `MTWSPy/data/gcmt_catalog.csv`

__MTWSPy/find_twin_obs.py__ &rarr; For observed & synthetic data: get travel times, (obs: primary quality control: snr_pow, snr_amp), detect time windows, filter time windows, plot time windows, save to twin file

__MTWSPy/match_twin_files.py__ &rarr; Find matching twin files between observed and synth data, save to file

__MTWSPy/phase_association.py__ &rarr; Associate unique observed & synthetic twin with unique phase within a given limit of each time window, save to twin file

__MTWSPy/correlate_twin.py__ &rarr; Correlate observed and synth time windows using cross-correlation of Zaroli et al (2010) formulation, select windows based on given quality criteria, save time delay files.

__MTWSPy/MTWSPy_main.py__ &rarr; Execute main code sequentially using params_in.yaml

__MTWSPy/run_MTWSPy.bash__ &rarr; Execute code on cluster using SLURM batch scheduler

__MTWSPy/post_processing/process_tdl_files.py__ &rarr; Time delay processing code to write station means or phase difference times from .tdl files produced by main code

__MTWSPy/post_processing/compare_picks.py__ &rarr; Time delay processing code to compare datasets using that from different parameters produced by main code

__MTWSPy/post_processing/plot_phase.gmt__ &rarr; Time delay plotting code

<!-- ######################################################################## -->

***
## Credits
***

The original serial code written in __matlab__ by Lei Li is available on request to stephanie.durand@ens-lyon.fr

<!-- ######################################################################## -->

***
## License
***

MIT License

Copyright (c) 2024 ALISTAIR BOYCE

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
