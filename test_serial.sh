#!/bin/bash

echo "Serial testing"

for case in 32 64 128 256 512 1024 2048
do
    echo "Case $case x $case";
    time -f %e sh -c "cat ./test_cases/$case.txt | ./bin/serial > /dev/null";
done