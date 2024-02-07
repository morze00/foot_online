#!/bin/bash                                                                    
# Initialize variables for command line arguments
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/u/land/fake_cvmfs/10/cvmfs_fairroot_v18.8.0_fs_nov22p1_extra/ucesb/hbook
cp ../event_display_dev/pedestal.dat .
cd /u/land/r3broot/202402_s091_s118/R3BParams_s091_s118/macros/exp/online/foot/event_display/src
clear ; clear
make clean; make
./event_display
