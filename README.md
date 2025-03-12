# MFT Dead Maps
A tool to visualize the content of MFT dead maps for a given run (or a group of runs)

## First steps:
1. Connect to LXPLUS: `ssh -X <user>@lxplus`
2. Load the ALICE environment. For example: `/cvmfs/alice.cern.ch/bin/alienv enter VO_ALICE@O2Physics::daily-20241102-0000-1` (this is a release from November 2024, any newer one can also be picked up). Alternatively: 
  - `export PATH=/cvmfs/alice.cern.ch/bin:$PATH`
  - `alienv enter VO_ALICE@O2Physics::daily-20241102-0000-1`
3. Initialize the GRID token: `alien-token-init <user>` (if not yet done, one first needs to export the GRID certificate into the `.globus` directory etc.: https://alice-doc.github.io/alice-analysis-tutorial/start/cert.html)

## The first time the script is used:
1. Clone the repo in LXPLUS: `git clone https://github.com/grundd/mft-deadmaps.git`
2. Enter the folder: `cd mft-deadmaps/`
3. Allow the Bash script: `chmod +x run_deadmaps.sh`

## How to run the script:
1. Log in to LXPLUS and go to the deadmaps folder: `cd mft-deadmaps/`
2. Either run the ROOT macro directly: `root -b -l -q mft_deadmaps.cxx\(<runNumber>\)`
3. Or use the Bash script: `./run_deadmaps.sh`. It also allows to automatically loop over many runs at the same time, but VS Code (or similar) is preferred to edit the script remotely using SSH.

## Output
Output plots will be stored in the subfolder `mft-deadmaps/<runNumber>/`. 

Example of a trending plot:
![trend](./readme/example_trend.png)

Example of a plot showing fractions of runtime when the chips were dead (only unmasked chips are shown):
![fractions](./readme/example_fractions.png)

To download the plots from LXPLUS to the local machine, either do it directly in VS Code, or:
1. Create a zip file on LXPLUS: `tar cvzf <runNumber>.zip <runNumber>/`
2. Open a different terminal window in the local machine and execute: `scp <user>@lxplus:/afs/cern.ch/user/<letter>/<user>/mft-deadmaps/<runNumber>.zip Downloads/`

