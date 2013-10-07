#!/bin/bash
## author: arturos@cern.ch
## Ceation: September 30th 2013
## Update:  October    7th 2013

## User Options to Setup:
FILE="samples.txt"         ## The name of the text file where the list of files are saved.
USER="arturos"             ## User who will send the GRID jobs.
PRODUCTION="13.2"          ## Tag of the production version.
TAR="yes"                  ## To create the tar file.               Possible answer: "yes" or "not"
LAUNCH="not"               ## To really execute the production now. Possible answer: "yes" or "not"
EXCLUDE="ANALY_SARA"       ## Site to exclude due to known problems. (In case of any, just write "")  Ex = ANALY_INFN-NAPOLI,ANALY_SARA

## End of User Options, please do not chance the lines bellow if you are not sure of the procedures.


#############################################################################################
echo "        ## ";date;                                                                   ##
echo "           "                                                                         ##
                                                                                           ##
./read_file.sh $FILE $USER $PRODUCTION $TAR $EXCLUDE >    production_HZZllqq_$PRODUCTION;  ##
## cat                                                    production_HZZllqq_$PRODUCTION;  ##
chmod 755                                                 production_HZZllqq_$PRODUCTION;  ##
if [ $LAUNCH = "yes" ];                                                                    ##
then                                                                                       ##
   nohup                                                ./production_HZZllqq_$PRODUCTION;  ##
fi                                                                                         ##
rm    xrootd_HZZllqq_$PRODUCTION    outs_HZZllqq_$PRODUCTION    Jobs_HZZllqq_$PRODUCTION;  ##
mv               production_HZZllqq_$PRODUCTION                 Jobs_HZZllqq_$PRODUCTION;  ##
                                                                                           ##
./create_outputs.sh   $FILE $USER $PRODUCTION $TAR $EXCLUDE >   outs_HZZllqq_$PRODUCTION;  ##
./find_xrootd_path.sh $FILE $USER $PRODUCTION $TAR $EXCLUDE > xrootd_HZZllqq_$PRODUCTION;  ##
chmod 755                                                       outs_HZZllqq_$PRODUCTION;  ##
chmod 755                                                     xrootd_HZZllqq_$PRODUCTION;  ##
echo "        ## Done!"                                                                    ##
#############################################################################################

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
