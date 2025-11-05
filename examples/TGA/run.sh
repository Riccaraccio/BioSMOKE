#!/bin/bash
if [ ! -f ../../build/source/BioSMOKEpp ]; then
    echo "Executable not found: Compile the project first"
    exit 1
fi

./../../build/source/BioSMOKEpp --input input.dic
