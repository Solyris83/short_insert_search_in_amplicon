# short_insert_search_in_amplicon
This repo is created for a custom analysis which involves sequencing target sequence which are very short (11-23bp long) being amplified by a set of primers from the system of interest (bacteria culture) and finding the expression level of the sequence of interest. 
HOW TO USE
1. Install miniconda/anaconda
2. Create conda environment with yaml file
    - conda env create --file environment.yml 
3. Change directory to work directory
    - cd /mnt/volume1/
4. Ensure data is found inside fastq folder ie /mnt/volume1/fastq
5. Run nextflow script
    - nextflow run custom_trimFASTP_mapPYTHON.nf -resume --outdir results -with-report log.html
