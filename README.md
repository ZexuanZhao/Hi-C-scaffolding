# HiFi-HiC-genome-assembly
 A snakemake workflow that assemble diploid genome using HiFi and HiC sequencing reads

# QC step

Quality control is not implemented. Usually adapters have been removed from HiFi reads. 
If they are not removed, they will be shown as overrepresented sequences in `fastQC`, 
along with telomeric repeats and other common repeats.

# To run
`snakemake --cores [cpu] --use-conda`
