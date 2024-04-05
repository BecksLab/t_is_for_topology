# T is for Topology

Code associated with the topology generators manuscript. At a high level 
this code will do the relevant data cleaning and preparation as well as 
calculate the relevant benchmarking metrics.

## Repo Structure

### Folders

- `data/`: datasets used for analyses as well as processed data
- `lib/`: library of additional functions used in the analyses

### Files

Note scripts are numbered sequentially but are also self contained

- `00_mangal_networks.jl`: download networks from Mangal database,
summarises Mangal networks
- `01_code.jl`: implements the various network generating families, also
summarises their outputs

## Running the code

All dependencies are documented in `Project.toml` and it is possible
to activate this environment using `Pkg.activate topology_generators`