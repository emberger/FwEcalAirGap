#!/bin/bash

# foldername=AngReso_IT20_OT20_Ogapfirst_25Inner_50Outer_lead1mm_Polystyrene5mm_20mmSteelVessel  #folder containing .root files
#
# FILES=/home/iwsatlas1/emberger/Geant4/Data/ForwardEcalWithPressureVessel/AngularResolution_OutlierRejection/$foldername
# ANALYSIS=$FILES/RejectionMethodAnalysis
#
# mkdir $ANALYSIS
#
# mkdir $ANALYSIS/ROOT
# mkdir $ANALYSIS/C
# mkdir $ANALYSIS/PDF
#
# counter=0
# for f in $FILES/*_.root
# do
#     #echo "$f"
#     FilePath=$f
#     filename="${FilePath##*/}"
#
#     foldername=${filename%_*}
#
#     # mkdir $AnaPath
#
#     ./Analysis 0.15 1000 $FilePath $ANALYSIS
#
# done
#
# # ./AngularResolution $FILES/RejectionMethodAnalysis/ROOT Res
# # ./PlotQuantiles $FILES/RejectionMethodAnalysis/ROOT Res
# #  ./EResAndShowerprofile 0.15 $FILES
#
# foldername=AngReso_IT20_OT20_Ogapfirst_25Inner_50Outer_lead2mm_Polystyrene5mm_20mmSteelVessel  #folder containing .root files
#
# FILES=/home/iwsatlas1/emberger/Geant4/Data/ForwardEcalWithPressureVessel/AngularResolution_OutlierRejection/$foldername
# ANALYSIS=$FILES/RejectionMethodAnalysis
#
# mkdir $ANALYSIS
#
# mkdir $ANALYSIS/ROOT
# mkdir $ANALYSIS/C
# mkdir $ANALYSIS/PDF
#
# counter=0
# for f in $FILES/*_.root
# do
#     #echo "$f"
#     FilePath=$f
#     filename="${FilePath##*/}"
#
#     foldername=${filename%_*}
#
#     # mkdir $AnaPath
#
#     ./Analysis 0.15 1000 $FilePath $ANALYSIS
#
# done
#
# foldername=AngReso_IT20_OT20_Ogapfirst_40Inner_50Outer_lead1mm_Polystyrene5mm_20mmSteelVessel  #folder containing .root files
#
# FILES=/home/iwsatlas1/emberger/Geant4/Data/ForwardEcalWithPressureVessel/AngularResolution_OutlierRejection/$foldername
# ANALYSIS=$FILES/RejectionMethodAnalysis
#
# mkdir $ANALYSIS
#
# mkdir $ANALYSIS/ROOT
# mkdir $ANALYSIS/C
# mkdir $ANALYSIS/PDF
#
# counter=0
# for f in $FILES/*_.root
# do
#     #echo "$f"
#     FilePath=$f
#     filename="${FilePath##*/}"
#
#     foldername=${filename%_*}
#
#     # mkdir $AnaPath
#
#     ./Analysis 0.15 1000 $FilePath $ANALYSIS
#
# done

foldername=AngReso_IT20_OT20_Ogapfirst_40Inner_50Outer_lead2mm_Polystyrene5mm_20mmSteelVessel  #folder containing .root files

FILES=/home/iwsatlas1/emberger/Geant4/Data/ForwardEcalWithPressureVessel/AngularResolution_OutlierRejection/$foldername
ANALYSIS=$FILES/RejectionMethodAnalysis

mkdir $ANALYSIS

mkdir $ANALYSIS/ROOT
mkdir $ANALYSIS/C
mkdir $ANALYSIS/PDF

counter=0
for f in $FILES/*_.root
do
    #echo "$f"
    FilePath=$f
    filename="${FilePath##*/}"

    foldername=${filename%_*}

    # mkdir $AnaPath

    ./Analysis 0.15 1000 $FilePath $ANALYSIS

done
