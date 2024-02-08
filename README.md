_________________
# MTWSPy
_________________

MTWSPy is a Python implementation of the Morphological Time Window Selection (MTWS) method first published in _Geophysical Journal International_ by Lei Li et al., ([2023](https://doi.org/10.1093/gji/ggad338)), including preceeding data processing workflow.

***

_________________
## Description
_________________

This package includes codes to download and process seismic data, create corresponding synthetics and measure travel times using MTWSPy.

* We use [obspyDMT](https://kasra-hosseini.github.io/obspyDMT/) for download of seismic data.
* We use [SPECFEM3D_GLOBE](https://specfem.org/) v8.0.0 available [here](https://github.com/SPECFEM/specfem3d_globe/releases/tag/v8.0.0) for calculation of global synthetic seismograms.
* We use [parallel](https://www.gnu.org/software/parallel/), [SAC](http://www.iris.edu/dms/nodes/dmc/software/downloads/sac), [saclst](http://geophysics.eas.gatech.edu/people/zpeng/Software/sac_msc.tar.gz) and [ttimes](https://ds.iris.edu/pub/programs/iaspei-tau/) for observed and synthetic data processing.

The codes are designed to be run on a year-by-year basis to build up a travel time database progressively, thereby avoiding data storage issues where possible.

HPC operations are setup using an example from the PSMN cluster at ENS Lyon (France) using the SLURM batch scheduler.

Parallelisation is used to speed up post-processing of observed and synthetic data (GNU parallel) and arrival time picking (python concurrent.futures) using the MTWS algorithm.

<!-- 
## Table of Contents

__This will also be bold__
_This will also be italic_ -->

***

_________________
## Installation
_________________

1. Create a fresh python3 environment (used 3.12) with the following pacakges available ('obspy','numpy', 'scipy', 'pandas', 'matplotlib', 'sys', 'glob', 'shutil', 'os', 'time', 'warnings', 'datetime', 'inspect', 'yaml'). Obspy 1.4 installed using pip.
2. Ensure obspyDMT is installed and active on the system used for data download. (Version 2.2.11 available using pip)
3. Ensure SPECFEM3D_GLOBE can be compiled and run using parameters appropriate to the system used to compute synthetics (typically HPC). v8.0.0 installed fresh from Github
4. Ensure GNU parallel, SAC, saclst and ttimes are available on system used for data processing. Small programs are available in the utils directory.
5. Ensure the following 10 files are copied into a local bin directory and are available on your PATH:
    * get_julday.py
    * get_CMTSOLUTION_tshift.py
    * obspy_rotate_comps.py
    * obspyDMT_rem_inst_resp.py
    * obspyDMT_data_par_proc_all.sh
    * obspyDMT_data_proc_event.sh
    * parallel_proc_specfem_seis.sh
    * proc_specfem_seis_file.sh
    * check_obspyDMT_SPECFEM_install.sh
    * check_obspyDMT_SPECFEM_install.py

***

_________________
## Usage
_________________

### Data Download

navigate to `dmt`
get_obspy_data.sh
proc_obspy_data.sh
mk_SPECFEM_STATIONS_dmt.py __OR__ mk_SPECFEM_STATIONS_obspy.py

### Generate Synthetics

navigate to `gcmt`

download_cmtsolutions.py

navigate to `specfem`

Par file
Station file
get_specfem_synthetics.sh
proc_specfem_synthetics.sh

### MTWSPy _Main_

navigate to `MTWSPy`

#### Contents:

* __params_in.yaml__ &rarr; Parameter file
* __v01_phasenames.py__ &rarr; Dictionary of phases desired
* __toolkit.py__ &rarr; Global functions
* __mk_events_csv.py__ &rarr; Make useable csv file of possible events
* __find_twin_obs.py__ &rarr; Find twins for observed data
* __find_twin_syn.py__ &rarr; Find twins for synth data
* __match_twin_files.py__ &rarr; Find matching twin files obs/synth 
* __phase_association_obs.py__ &rarr; Associate observed twin with a phase
* __phase_association_syn.py__ &rarr; Associate synth twin with a phase
* __correlate_twin.py__ &rarr; Correlate observed and synth time windows
* __MTWSPy_main.py__ &rarr; Execute main code

#### Description:

The overall philospy of the code is to have a setup where we can just hit `GO`.

To do this we rely heavily on an input parameter file `params_in.yaml` which must be studied and adapted for purpose before using the code. All parameters are described/commented and are grouped with each relevant part of the code. This file is loaded before the execution of any part of the code.

-> open params_in.yaml
-> Check home, data_loc, synth_loc, fmt_data_loc, cmt_infile, year, parallel, cores

Where possible the code also relies on a toolkit of global functions in `toolkit.py` these should not be changed but are crucial to the effective running of the code (e.g., initialising parallelisation of commands).

Although each step of the code can be launched individually (`python XXX.py`) recommended usage is carried out by executing the full code: `python MTWSPy.py` 

All steps of the code that utilise parallelisation (also written to work in serial) have a common input format:

* __main_function__ &rarr; Usually `process_one_event`, checks/prepares specific inputs for each function(_not to be changed_)
* __input_directory__ &rarr; Location of input data (_set based on parameter file_)
* __evt_id_tab__ &rarr; Loaded events (_pandas dataframe_)
* __functions__ &rarr; List of specifc functions for given code step (_unlikely to be changed_)
* __params_in__ &rarr; Loaded parameter file (_dict_)
* __phases__ &rarr; Loaded phase file (_dict_)



#### Detailed description:

__params_in.yaml__ &rarr; Parameter file containing all variables that can be changed within the code.

__v01_phasenames.py__ &rarr; Dictionary of phases desired, could be modified to search for alternative seismic phases

__toolkit.py__ &rarr; Global function definitions used throughout the code to speed up execution.

__mk_events_csv.py__ &rarr; Make useable csv file of possible earthquakes starting with gcmt catalog file at `gcmt/gcmt_catalog.csv`

__find_twin_obs.py__ &rarr; For observed data: get travel times, primary quality control (snr_pow, snr_amp), detect time windows, filter time windows, plot time windows, save to twin file

__find_twin_syn.py__ &rarr; For synth data: get travel times, detect time windows, filter time windows, plot time windows, save to twin file

__match_twin_files.py__ &rarr; Find matching twin files between observed and synth data, save to file

__phase_association_obs.py__ &rarr; Associate observed twin with a phase within a given limit of each time window, save to twin file

__phase_association_syn.py__ &rarr; Associate synth twin with a phase within a given limit of each time window, save to twin file

__correlate_twin.py__ &rarr; Correlate observed and synth time windows, select windows based on given quality criteria, save time delay files.

__MTWSPy_main.py__ &rarr; Execute main code in sequence.



Execute using batch scheduler? 1 node 96 cores...?




<!-- 
## Installation

Use the package manager [pip](https://pip.pypa.io/en/stable/) to install foobar.

```bash
pip install foobar
```

## Usage

```python
import foobar

# returns 'words'
foobar.pluralize('word')

# returns 'geese'
foobar.pluralize('goose')

# returns 'phenomenon'
foobar.singularize('phenomena')
``` -->

***

_________________
## Credits
_________________

The original serial code written in __matlab__ by Lei Li is available on request to stephanie.durand@ens-lyon.fr

***

_________________
## License
_________________

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