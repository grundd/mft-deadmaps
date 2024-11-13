# chmod +x run_deadmaps.sh
#!/bin/bash

runs=(
  559781
  559782
  559783
  559784
  559802
  559803
  559804
  559827
  559828
  559830
  559843
  559856
  559899
  559901
  559902
  559903
  559917
  559919
  559920
  559933
)
threshold=0. # already in percent
output_name=2024_PbPb_runs

for run in "${runs[@]}"
do
  echo "Run $run:"
  root -b -l -q 'mft_deadmaps.cxx('$run','$threshold')'
done

echo "Adding to .zip archive: ${runs[@]}"

tar cvzf $output_name.zip "${runs[@]}"
