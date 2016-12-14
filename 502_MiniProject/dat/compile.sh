#!/bin/bash
 
###########################################################
#
#
# M9 sim setup
#
#
###########################################################
 
# compile m9 sim
echo "Starting compilation"
echo "Reminder: Make sure root is sourced"
#rm mtest2

rm go
#macro2
g++ main.cpp vroot/root.cpp -o2 -o go `root-config --cflags --glibs` -std=c++0x -pthread
#macro1

#g++ muon.cpp gPT/GPT.cpp mROOT/mroot.cpp -o2 -o mugo `root-config --cflags --glibs` -std=c++0x -pthread

echo "done!"
 

