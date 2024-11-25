#!/bin/bash

alglib="alglib-cpp/src/"

for i in "$alglib"/*.cpp
do
g++ -std=c++0x -c -o "$i".o "$i"
done

for i in ./*.cpp
do
g++ -std=c++0x -c -o "$i".o "$i"
done


for i in ./*.c
do
gcc -std=c99 -c -o "$i".o "$i"
done


g++ -o PriorC ./*.o "$alglib"/*.o


rm -r ./*.o "$alglib"/*.o
