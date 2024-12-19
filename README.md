## Introduction

**EMC/NBKsixteenS** is a bioinformatics pipeline designed for processing 16S sequences generated using Oxford Nanopore Technologies (ONT). The pipeline includes host depletion, making it suitable for experiments that use the Native Barcoding Kit (NBK) instead of ONT's dedicated 16S kit.

As input it requires a samplesheet with paths to long read, compressed, FASTQ files. The pipeline performs quality control and trimming on the reads, filters out reads mapping to a specified host reference genome and taxonomically classifies the remaining reads.

<!-- TODO nf-core: Include a figure that guides the user through the major workflow steps. Many nf-core
     workflows use the "tube map" design for that. See https://nf-co.re/docs/contributing/design_guidelines#examples for examples.   -->
<!-- TODO nf-core: Fill in short bullet-pointed list of the default steps in the pipeline -->

NBKsisteenS includes:
1. Raw read quality control ([`FastQC`](https://github.com/s-andrews/FastQC))
2. Adapter trimming ([`Porechop`](https://github.com/rrwick/Porechop))
3. Filtering by quality ([`Filtlong`](https://github.com/rrwick/Filtlong))
4. Visualize QC'ed data ([`NanoPlot`](https://github.com/wdecoster/NanoPlot))
5. Filter out reads mapping to a reference genome ([`minimap2`](https://github.com/lh3/minimap2) and [`Samtools`](https://www.htslib.org/doc/samtools-view.html))
6. Convert SAM file to FASTQ ([`Samtools`](https://www.htslib.org/doc/samtools-fasta.html))
7. Taxonomic classification ([`Kraken2`](https://github.com/DerrickWood/kraken2))
8. Generate report with quality metrics and used tools ([`MultiQC`](https://github.com/MultiQC/MultiQC))

To be added:
- Summarize mapping statistics ([`Samtools`](https://www.htslib.org/doc/samtools-flagstat.html))
- Visualize Kraken2 output with Krona ([`KrakenTools`](https://github.com/jenniferlu717/KrakenTools) and [`Krona`](https://github.com/marbl/Krona))
- Re-estimation of microbial abundances ([`Bracken`](https://github.com/jenniferlu717/Bracken))

## Usage
To use NBKsixteenS on your machine, follow the steps below:
1. Make sure you have correctly set-up Nextflow and it's dependencies
> [!NOTE]
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/usage/installation) on how to set-up Nextflow. 

2. Clone this GitHub repository
3. Prepare a samplesheet like the example below:
    `samplesheet.csv`:
    ```csv
    sample,fastq_1
    CONTROL_1,BR_PVP_0705.fastq.gz
    ```
    Each row represents a fastq file.
<!-- 
> [!TIP]
> If you don't have data available yet, or you want to test the pipeline first on a small dataset, use the [data that comes with this repo](https://github.com/BirgitRijvers/EMC-MetaMicrobes/tree/master/testdata). This data is subsampled from 3 RNA-seq samples with varying host contents, created by Marques *et al.* .
-->

> [!TIP]
> You can use the ["samplesheeter.py"](https://github.com/BirgitRijvers/EMC-MetaMicrobes/tree/master/samplesheeter) script that comes with another repo, a small command line tool that prepares the samplesheet for you based on a supplied data directory.

   <!-- TODO nf-core: Add documentation about samplesheeter and testdata -->
4. Download a FASTA file containing the reference genome you want to use for host depletion, for example [GRCh38](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001405.26/). 

Optionally, create a minimap2 index of this reference file and build your preferred Kraken2/Bracken database. If you don't supply these to the pipeline, NBKsixteenS will index your reference genome for you and build the Kraken2/Bracken standard database. 

4. Now, you can run the NBKsixteenS pipeline using:
   
    ```bash
    nextflow run <path/to/EMC-NBKsixteenS/directory/> \
       -profile <docker/singularity/conda/.../institute> \
       --input samplesheet.csv \
       --outdir <OUTDIR> \
       --fasta <path/to/reference_genome_fasta>
    ```
    
    If you have a pre-built minimap2 index or Kraken2/Bracken database, use a command like this:
      
      ```bash
      nextflow run <path/to/EMC-NBKsixteenS/directory/> \
         -profile <docker/singularity/conda/.../institute> \
         --input samplesheet.csv \
         --outdir <OUTDIR> \
         --fasta <path/to/reference_genome_fasta> \
         --minimap2_index <path/to/bwa_mem2_index> \
         --kraken2_db <path/to/kraken2_db> \
         --bracken_db <path/to/bracken_db>
      ```  

> [!WARNING]
> Please provide pipeline parameters via the CLI or Nextflow `-params-file` option.
> see [docs](https://nf-co.re/usage/configuration#custom-configuration-files).

## Credits

EMC/nbksixteens was originally written by BirgitRijvers.

## Support

Please open an issue in this repository if you experience problems or have development suggestions.

## Citations

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

This pipeline uses code and infrastructure developed and maintained by the [nf-core](https://nf-co.re) community, reused here under the [MIT license](https://github.com/nf-core/tools/blob/master/LICENSE).

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
