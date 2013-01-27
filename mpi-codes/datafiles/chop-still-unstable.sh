#!/bin/bash
#reads a datafile, prints only after >1 periods cease to appear
#Assumes datafiles are going from grazing to non-grazing 
#If otherwise, add a reverse

Start=`grep -nE "	1000$" $1 | awk '{FS=":"};{print $1}' | tail -n 1 `

tail -n +$((Start+1)) $1 | grep -v "#"
