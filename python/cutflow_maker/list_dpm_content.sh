# properly lists a directory in ATLAS disk, giving a list of files
# author: Valerio Ippolito <valerio.ippolito@cern.ch>

#!/bin/bash

the_path="/dpm/roma1.infn.it/home/atlas/$1"

for the_dirname in `rfdir $the_path | awk '{print \$9}'`; do
  for the_filename in `rfdir $the_path/$the_dirname | awk '{print \$9}' | grep root`; do
    echo "root://grid-cert-03.roma1.infn.it/$the_path/$the_dirname/$the_filename"
  done
done
