#!/bin/bash
#echo -e '\0033\0143'
#echo -e "Activating FEMfuns environment"

# This usually is not needed if already done before calling matlab
#source ~/.bashrc
#conda activate femfuns

#echo -e '\0033\0143'
echo "Converting mesh to FEMfuns format"
dolfin-convert $1 $2

python3.7 femfuns_caller.py $2