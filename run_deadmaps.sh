# chmod +x run_deadmaps.sh
#!/bin/bash

run=559211

root -b -l -q 'mft_deadmaps.cxx('$run')'
tar cvzf $run.zip $run/
