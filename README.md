# Replication code for "Differential privacy does not prevent lawful redistricting"

This repository contains code to:
1. Generate shapefiles for use in redistricting simulations.
2. Run MCMC redistricting chains using the ReCom algorithm (implemented in the GerryChainJulia software package [1]).
3. Reproduce results found in the paper "Differential privacy does not prevent lawful redistricting"

## Requirements

* Python 3.11
    * To install required packages: `pip install -r requirements.txt`
* Julia 1.8
    * Follow installation instructions at https://julialang.org/downloads/.
    * Install environment dependencies:
        * Run the Julia REPL while in the `redistricting_with_DP` directory.
        * Hit the `]` key to enter the Pkg REPL.
        * Run the following commands:
            * `(@v1.8) pkg> activate .`
            * `(redistricting_with_DP) pkg> instantiate`
        * This will automatically install all Julia dependencies. Further instructions on package management in Julia can be found at https://pkgdocs.julialang.org/.
 
## Preprocessing

The script `preprocessing/process_state.sh` takes in a state, a district type, a geolevel, and a data path, and downloads the corresponding TIGER/Line shapefile and 2020 DAS Demonstration dataset. It then calls `preprocessing/make_jl_shapefile.py`, which creates a new shapefile containing the Demonstration data.

## Ensemble Generation

The Julia script `julia/run_recom.jl` will sample an ensemble of a given length using the ReCom algorithm and save district tabulations to a csv file.

The Julia script `julia/run_recom_critical_offset.jl` will run a series of ReCom chains at varying offsets to determine the critical offset for a given geography (the minimum population tolerance offset for which less than 2% of plans exceed a given threshold).

The shell script `experiments.sh` contains 4 experiments which should be run before running the analysis notebooks.
* `./experiments.sh 1` will generate many ensembles for the Louisiana State Senate, using different population tolerance thresholds and standard redistricting constraints.
* `./experiments.sh 2` will generate one larger ensemble for the Louisiana State Senate, using a 5% population tolerance threshold and standard redistricting constraints.
* `./experiments.sh 3` will generate one ensemble for the Georgia State House, using a 5% population tolerance threshold and short burst optimization to increase the number of majority-Black districts in each sampled plan.
* `./experiments.sh 4` will compute the critical offset and no-offset discrepancy rate for the Georgia State House, using a 5% population tolerance threshold and standard redistricting constraints.

## Analysis

The `analysis` directory contains three Python notebooks to reproduce results from the paper:
* `LA_SS_offsets.ipynb` uses the ensembles generated by `./experiments.sh 1` to reproduce Figure 1 (Discrepancy with offsets in the Louisiana state senate).
* `LA_SS_pop_bal.ipynb` uses the ensemble generated by `./experiments.sh 2` to reproduce Figure 2 (The two components of population deviation in districts).
* `GA_SH_gingles.ipynb` uses the ensemble generated by `./experiments.sh 3` to reproduce Figures 3 (MMD discrepancies in the Georgia state lower house) and 4 (MMD discrepancy rate by BVAP margin for the Georgia state house).

## References

[1] bsuwal, Max Fan, and Matthew Sun. “Mggg/gerrychainjulia: Minor Fixes + Save as HDF5”. Zenodo, March 31, 2021. https://doi.org/10.5281/zenodo.4649464.