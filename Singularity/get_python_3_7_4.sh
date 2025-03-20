#!/usr/bin/bash

## Downloads a python v3.7.4 image from https://depot.galaxyproject.org/singularity/
## Assume Singularity is already installed

singularity pull https://depot.galaxyproject.org/singularity/python:3.7

singularity build my_python3.7.sif python:3.7


