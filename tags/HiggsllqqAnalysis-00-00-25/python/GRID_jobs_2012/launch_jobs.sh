#!/bin/bash
## author: arturos@cern.ch
## Ceation: September 30th 2013
## Update:  November  28th 2013

## User Options to Setup:
FILE="SamplesD3PD.txt"     ## The name of the text file where the list of files are saved.
USER="arturos"             ## User who will send the GRID jobs.
PRODUCTION=$1              ## "16.0"                  ## Tag of the production version.
TAR="yes"                  ## To create the tar file.               Possible answer: "yes" or "not"
LAUNCH="not"               ## To really execute the production now. Possible answer: "yes" or "not"
EXCLUDE="ANALY_SARA"       ## Site to exclude due to known problems. (In case of any, just write "")  Ex = ANALY_INFN-NAPOLI,ANALY_SARA
SYST=$2                    ## NOW: 28thNov2013 => this will corespond to a subtag of the prod.
## End of User Options, please do not chance the lines bellow if you are not sure of the procedures.


#########################################################################################################
echo "        ## ";date;                                                                               ##
echo "           "                                                                                     ##
                                                                                                       ##
./read_file.sh $FILE $USER $PRODUCTION $TAR $EXCLUDE $SYST >    production_HZZllqq_$PRODUCTION.$SYST;  ##
## cat                                                          production_HZZllqq_$PRODUCTION.$SYST;  ##
chmod 755                                                       production_HZZllqq_$PRODUCTION.$SYST;  ##
if [ $LAUNCH = "yes" ];                                                                                ##
then                                                                                                   ##
   nohup                                                      ./production_HZZllqq_$PRODUCTION.$SYST;  ##
fi                                                                                                     ##
rm xrootd_HZZllqq_$PRODUCTION.$SYST       outs_HZZllqq_$PRODUCTION.$SYST;                              ##
rm   Jobs_HZZllqq_$PRODUCTION.$SYST   download_HZZllqq_$PRODUCTION.$SYST;                              ##
mv               production_HZZllqq_$PRODUCTION.$SYST                 Jobs_HZZllqq_$PRODUCTION.$SYST;  ##
                                                                                                       ##
./create_outputs.sh   $FILE $USER $PRODUCTION.$SYST $TAR $EXCLUDE >   outs_HZZllqq_$PRODUCTION.$SYST;  ##
./find_xrootd_path.sh $FILE $USER $PRODUCTION.$SYST $TAR $EXCLUDE > xrootd_HZZllqq_$PRODUCTION.$SYST;  ##
chmod 755                                                             outs_HZZllqq_$PRODUCTION.$SYST;  ##
chmod 755                                                           xrootd_HZZllqq_$PRODUCTION.$SYST;  ##
sed 's/user\./dq2-get\ user\./g' outs_HZZllqq_$PRODUCTION.$SYST > download_HZZllqq_$PRODUCTION.$SYST;  ##
chmod 755                                                         download_HZZllqq_$PRODUCTION.$SYST;  ##
echo "        ## Done!"                                                                                ##
ls -lhrt --col                                                                                         ##
#########################################################################################################

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
