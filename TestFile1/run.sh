#!/bin/sh
# build the program!
# note: there will eventually be a separate build step for your bot, but for now it counts against your runtime.
LIBRARIES="-lutil -ldl -lrt -lpthread -lgcc_s -lc -lm -L/battlecode-c/lib/ -lbattlecode"
INCLUDES="-I/battlecode-c/include"
g++ main.cpp -o main -std=c++11  $LIBRARIES $INCLUDES

# run the program!
./main
