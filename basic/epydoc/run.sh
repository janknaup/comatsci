#!/bin/bash

packAndRotate()
{
        local packAndRotateFile
        local packAndRotateI
        # iterate through list of arguments
        for packAndRotateFile in $@
        do
                # just skip nonexisting files
                if [[ -f ${packAndRotateFile} ]]
                then
                        # first check if the highest file numer we will rotate to already exists
                        # in that case abort
                        if [[ -f ${packAndRotateFile}.9.gz ]]
                        then
                                break
                        fi
                        # check backwards through preexisting rotated files
                        # and advance rotation
                        for ((packAndRotateI=8; packAndRotateI; packAndRotateI--))
                        do
                                if [[ -f ${packAndRotateFile}.${packAndRotateI}.gz ]]
                                then
                                        mv -fv ${packAndRotateFile}.${packAndRotateI}.gz ${packAndRotateFile}.$((${packAndRotateI}+1)).gz
                                fi
                        done
                        # rotate the unnumbered packed version if it exists
                        if [[ -f ${packAndRotateFile}.gz ]]
                        then
                                mv -fv ${packAndRotateFile}.gz ${packAndRotateFile}.1.gz
                        fi
                        # finally pack the target file
                        gzip -9 ${packAndRotateFile}
                fi
        done
        return
}

packAndRotate epydoc.log

epydoc --html --pdf --show-imports --graph=all -v ../../test_install/lib/python???/site-packages/comatsci ../../test_install/lib/python???/site-packages/geostatspack4 > epydoc.log

