#!/bin/bash

touch icosphere_time.txt
touch brute_force_time.txt

orders=(1 2 3 4 5 6)

for order in "${orders}"; do
  echo "$order" | /usr/bin/time -v ./brute_force.exe >> brute_force_time.txt
  echo "$order" | /usr/bin/time -v ./icosphere.exe >> icosphere_time.txt
done

