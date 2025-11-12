# Start an analysis with ViMOP using the command line

## Run ViMOP

1. Download our demo test set

```bash
wget https://opr.bnitm.de/example_data/vimop-demo.tar.gz
tar -xvzf vimop-demo.tar.gz
```

2. run ViMOP demo with default settings

```bash
nextflow run OPR-group-BNITM/vimop --fastq vimop-demo/lasv_simulated -o vimop-demo/output
```

To show all available options run:

```bash
nextflow run OPR-group-BNITM/vimop --help
```
