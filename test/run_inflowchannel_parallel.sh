#!/bin/bash

OFID=$1
BASHSCR=$2
DATADIR=$3
shift 3 # remove CMD line args, causes problems with some bashrcs

source ${BASHSCR}
if [ -d $OFID ]; then rm -rf $OFID; fi; mkdir $OFID && cd $OFID && (

analyze --workdir . $DATADIR/inflowchannel_parallel.ist

)
