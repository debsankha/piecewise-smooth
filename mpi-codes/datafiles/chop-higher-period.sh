#!/bin/bash
#reads a datafile, prints only after >1 periods cease to appear
#Assumes datafiles are going from grazing to non-grazing 
#If otherwise, add a reverse

Start=`grep -nE "Period: (([2-9])|([0-9]{2,}))" $1 | awk '{FS=":"};{print $1}' | tail -n 1 `

tail -n +$Start $1 
