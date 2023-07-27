#!/bin/bash
set -e

# cat in2.json | go run boundaries_setup.go -t 8 > out.txt 2>&1

# cat in2.json | go run boundaries_setup.go -nocload -t 8 > out2.txt 2>&1
cat in3.json | go run boundaries_setup.go -t 8 > out3.txt 2>&1
