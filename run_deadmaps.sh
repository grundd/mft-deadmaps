# chmod +x run_deadmaps.sh
#!/bin/bash

# 559781
# 559782
# 559783
# 559784
# 559802
# 559803
# 559804
# 559827
# 559828
# 559830
# 559843
# 559856
# 559899
# 559901
# 559902
# 559903
# 559917
# 559919
# 559920 
# 559933
# 559966
# 559968
# 559969
# 559970
# 559987
# 560012
# 560031
# 560033
# 560034
# 560049
# 560066
# 560067
# 560070
# 560089
# 560090

runs=(560105)

threshold=0. # fraction, not percent!
output_name=560105 # 2024_PbPb_runs

for run in "${runs[@]}"
do
  echo "Run $run:"
  root -b -l -q 'mft_deadmaps.cxx('$run','$threshold')'
done

#echo "Adding to .zip archive: ${runs[@]}"

#tar cvzf $output_name.zip "${runs[@]}"
