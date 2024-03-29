#!/bin/bash

# Note: this script is should be sourced from the project root directory

set -e

# define some variables
yr=$(date +%Y)
dt=$(date +%Y-%m-%d)
cprt_R=misc/copyright.R
# cprt_cpp=misc/copyright.cpp
# citation=inst/CITATION
version=$(grep "Version" DESCRIPTION | awk '{print $NF}')

if [ "$(uname)" == "Darwin" ] || [ "$(uname)" == "Linux" ]; then
    printf "Updating date, version, and copyright year.\n"

    # update copyright year in the template headers
    regexp1="s/Copyright \(C\) [0-9]+/Copyright \(C\) $yr/"
    # regexp1="s/Copyright \(C\) 2020-[0-9]+/Copyright \(C\) 2020-$yr/"
    sed -i.bak -E "$regexp1" $cprt_R && rm $cprt_R.bak
    # sed "s_#_/_g" $cprt_R > $cprt_cpp

    # update copyright year in all R scripts
    for Rfile in R/*.R
    do
	if ! grep -q 'Copyright (C)' $Rfile; then
	    cat $cprt_R $Rfile > tmp
	    mv tmp $Rfile
	fi
	sed -i.bak -E "$regexp1" $Rfile && rm $Rfile.bak
    done

    # update copyright year in all C++ scripts
    # for cppfile in src/*.cpp inst/include/*.h inst/include/*/*.h
    # do
    #     if ! grep -q 'Copyright (C)' $cppfile; then
    #         cat $cprt_cpp $cppfile > tmp
    #         mv tmp $cppfile
    #     fi
    #     sed -i.bak -E "$regexp1" $cppfile && rm $cppfile.bak
    # done
    # rm $cprt_cpp

    # update date in DESCRIPTION
    regexp2="s/Date: [0-9]{4}-[0-9]{1,2}-[0-9]{1,2}/Date: $dt/"
    sed -i.bak -E "$regexp2" DESCRIPTION && rm DESCRIPTION.bak

    # update version and year in citation
    # regexp3="s/version ([0-9]+\.*)+/version $version/"
    # sed -i.bak -E "$regexp3" $citation && rm $citation.bak
    # restrict the search and only update the year of package
    # regexp4="/intsurv-package/,/^\)$/ s/20[0-9]{2}/$yr/"
    # sed -i.bak -E "$regexp4" $citation && rm $citation.bak

    # done
    printf "All updated.\n"
else
    printf "Remember to update date, version, and copyright year.\n"
fi
