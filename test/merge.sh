#!/bin/sh

if [ -z $ROOTSYS ]; then
    echo "ROOTSYS is not defined: source ROOT, or hadd won't work!"
    exit
fi

if ls ./Test_$1_p* &> /dev/null ; then
    ls ./Test_$1_p* 
    hadd -f Test_$1.root Test_$1_p*.root
    if ls Test_$1.root  &> /dev/null ; then
	rm Test_$1_p*.root
    fi
else
    echo "Pattern Test_$1 not found"
fi
