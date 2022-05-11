# refine.mcr: a YASARA macro energy minimize an AlphaFill model
# Runs the default YASARA energy minimization experiment. The returned model is protonated.  
#
# Version 0.01
#
# Minimum YASARA tier: YASARA dynamics 
#
# This script was created by Robbie P. Joosten (r.joosten@nki.nl)
# Reference: If you publish results (directly or indirectly) obtained by using this macro, please cite YASARA and (any of)
# these publications:
# 1) Maarten L. Hekkelman, Ida de Vries, Robbie P. Joosten, Anastassis Perrakis: "­AlphaFill: enriching AlphaFold models 
#    with  ligands and co-factors" (2022)
#
# Changelog:
# Version 0.01
# - First attempt.
#
# YASARA input variables:
# modelin  The mmCIF file of the Alphafill model
# modelout The minimised AlphaFill model
#
#Initialise
OnError Exit
Console off
Processors 8

#Check the YASARA version
if (Dynamics)==0
  Print "YASARA version too low, you need at least YASARA Dynamics"
  exit

#Load the model
LoadCIF (modelin), Center=No, Correct=No

#Setup the forcefield
Forcefield NOVA
Longrange None
Cutoff 10.5
Boundary Wall

#Do the energy minimization (default settings)
Experiment minimization
Experiment On

# Wait till end of experiment
Wait ExpEnd

#Write the coordinate file
SaveCIF 1, (modelout), Format=mmCIF, Transform=No

exit
