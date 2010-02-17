#!/bin/bash
##############################################################################
# make_ui.sh
# Part of PAth Search Tool bAsed on Flexible Atomistic Reaction Image ANalysis
# (c) 2005-2008 by Jan M. Knaup <Knaup@bccms.uni-bremen.de>
# all rights reserved
##############################################################################
# Licensed under the Non-Profit Open Software License version 3.0
# see file LICENSE for details.
##############################################################################
# script file to compile all qt .ui files in cwd to pyqt3 python scripts
##############################################################################

#
# configuration options
#

UIC=/usr/bin/pyuic

echo "Trying to compile all .ui files to python scripts"

# get the basenames of .ui files in cwd
baseNames=""
for i in *.ui
do
  TEMP=`basename ${i} .ui`
  baseNames="${baseNames} ${TEMP}"
done

echo "found scripts: ${baseNames}"

# compile forms, make backup copies of scripts, if already present
for i in ${baseNames}
do
  echo -n "."
  if [ -f ${i}.py ]
  then
    mv ${i}.py ${i}.py.bak
  fi
  if ! $UIC ${i}.ui > ${i}.py
  then
    echo "error compiling form ${i}"
    exit 1
  fi
done
echo "done"
# finished
