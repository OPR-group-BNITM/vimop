# Installation of ViMOP with command line
# Table of Contents

- [Prerequisites](#prerequisites)
  - [Docker installation](#docker-installation)
  - [Nextflow installation](#nextflow-installation)
- [Installation and setup of ViMOP](#installation-and-setup-of-vimop)
  - [Quick start](#quick-start)
  - [Step by step installation, set up and execution](#step-by-step-installation-set-up-and-execution)
    - [Clone Git repository](#clone-git-repository)
    - [Set up the database](#set-up-the-database)
  - [Run ViMOP](#run-vimop)
    - [Demo run](#demo-run)
    - [User input](#user-input)


# Prerequisites

ViMOP uses epi2me, nextflow and docker. All of these dependencies can be installed with and without the usage of the command line depending on the user's preference.

## Docker installation
You can find detailed installation steps for Docker in the [Docker docs](https://docs.docker.com/engine/install/). After following the steps on that web page for your operating system you should make Docker managable as non-root user and make it run on start up. For this, follow the steps on this website: [Post-install instructions](https://docs.docker.com/engine/install/linux-postinstall/).

## Nextflow installation
The detailed installation steps for the Nextflow installation can be found here: [nextflow documentation](https://nextflow.io/docs/latest/install.html#install-page).  

If you have conda installed, you can create a conda environment and install nextflow in there with:



```python
bash conda create -n nextflow nextflow
conda activate nextflow
```

# Installation and setup of ViMOP
To run ViMOP you need to install the application and set up the database.

## Quick start

If you want to run the pipeline without cloning our git repository you can just run the following steps
1. Setup database



```python
nextflow run OPR-group-BNITM/vimop --download_db_all
```


2. Run the pipeline with your own fastq and optionally your own configuration. To configure your run you can download the [nextflow.config](https://github.com/OPR-group-BNITM/vimop/blob/main/nextflow.config) in our GitHub repository and edit the parameters as needed.


To test the pipeline you can get some demo data [here](https://opr.bnitm.de/example_data/lasv_simulated_run.fastq)


```python
nextflow run OPR-group-BNITM/vimop --fastq /path/to/fastqfiles --out_dir /path/for/your/output [optional: -c /path/to/nextflow.config]
```

## Step by step installation, set up and execution

### Clone Git repository
To import the ViMOP workflow clone our GitHub repository.

1. Open the command line 
2. Clone the repository



```python
git clone https://github.com/OPR-group-BNITM/vimop.git
```

### Set up the database


1. Change directory into vimop folder


```python
cd vimop
```

2. Download the database


```python
nextflow run main.nf --download_db_all
```



If you want to replace an existing database with our latest version, add the option --download_db_update_existing  



```python
nextflow run main.nf --download_db_all --download_db_update_existing
```

## Run ViMOP

### Demo run

1. Download our demo test set


```python
wget https://opr.bnitm.de/example_data/vimop-demo.tar.gz
tar -xvzf vimop-demo.tar.gz
```

2. run ViMOP demo with default settings


```python
nextflow run main.nf --fastq vimop-demo/lasv_simulated -o vimop-demo/output
```

### User input

In general, you can run ViMOP with the following command:


```python
nextflow run main.nf  --fastq "/path/to/fastqfiles" --out_dir "/path/for/your/output"
```

To show all available options run:


```python
nextflow main.nf --help
```
