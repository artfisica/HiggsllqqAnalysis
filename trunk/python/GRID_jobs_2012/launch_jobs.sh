#!/bin/bash
## author: arturos@cern.ch
## September 30th 2013

## User Options to Setup:
FILE="samples.txt"         ## The name of the text file where the list of files are saved.
USER="arturos"             ## User who will send the GRID jobs.
PRODUCTION="14"            ## Tag of the production version.
TAR="yes"                  ## To create the tar file.               Possible answer: "yes" or "not"
LAUNCH="not"               ## To really execute the production now. Possible answer: "yes" or "not"

## End of User Options, please do not chance the lines bellow if you are not sure of the procedures.


#############################################################################
echo "        ## ";date; 
echo "           "

./read_file.sh $FILE $USER $PRODUCTION $TAR > production_HZZllqq_$PRODUCTION;
cat                                           production_HZZllqq_$PRODUCTION;
chmod 755                                     production_HZZllqq_$PRODUCTION;
if [ $LAUNCH = "yes" ]; 
then
   nohup                                    ./production_HZZllqq_$PRODUCTION;
fi
rm                                                  Jobs_HZZllqq_$PRODUCTION;
mv      production_HZZllqq_$PRODUCTION              Jobs_HZZllqq_$PRODUCTION;
echo "        ## Done!"
#############################################################################

## Example of "samples.txt" file ##
    ##########################
    ## Signal_low_mass.txt  ##
    ## Signal_high_mass.txt ##
    ## Z_jets.txt           ##
    ## W_jets.txt           ##
    ## Diboson.txt          ##
    ## Top.txt              ##
    ## DY_jets.txt          ##
    ## Data.txt             ##
    ##(other sample?)       ##
    ##########################
