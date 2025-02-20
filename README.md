# Nanoflow analysis pipeline

Analysis pipeline for virus metagenomics using the sispa-protocol.

## Requirements

TODO

## Installation

Todo

## Run nanoflow 

Todo

## Options

TODO

## TODO

- write this README :)
- write tests in groovy at the beginning of the pipeline run, e.g.
    - are the required database files in the given paths? (in config)
    - is there enough disk space?
    - does the machine have given requirements (CPUs, RAM ..)
        - add a flag run_anyways in configs?
- gzip files when writing to results (e.g. cleaned fastq)
  - also in between to keep work small?

- assembly
  - exclude bubbles
  - input coverage parameter that determines the number of reads sampled