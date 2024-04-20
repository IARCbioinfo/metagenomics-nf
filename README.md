# metagenomics-nf
## Nextflow pipeline to perform metagenomic analyses from next-generation sequncing data

[![Docker Hub](https://img.shields.io/badge/docker-ready-blue.svg)](https://hub.docker.com/r/iarcbioinfo/template-nf/)

## Description
Nextflow pipeline running software centrifuge to detect reads mapping to microbial or viral references, and optionally software virusbreakend to detect viral integration

## Dependencies

1. This pipeline is based on [nextflow](https://www.nextflow.io). As we have several nextflow pipelines, we have centralized the common information in the [IARC-nf](https://github.com/IARCbioinfo/IARC-nf) repository. Please read it carefully as it contains essential information for the installation, basic usage and configuration of nextflow and our pipelines.
2. External software:
- [centrifuge](https://ccb.jhu.edu/software/centrifuge/manual.shtml)
- [samtools](https://www.htslib.org/doc/samtools.html)
- [virusbreakend](https://github.com/PapenfussLab/gridss/blob/master/VIRUSBreakend_Readme.md)

You can avoid installing all the external software by only installing Docker. See the [IARC-nf](https://github.com/IARCbioinfo/IARC-nf) repository for more information.

To use virusbreakend, you will need to download the database. See help at https://github.com/PapenfussLab/gridss/blob/master/VIRUSBreakend_Readme.md

## Input
  | Type      | Description     |
  |-----------|---------------|
  | input_folder   | Input folder with BAM or CRAM files |

## Parameters

  * #### Optional
| Name      | Default value | Description     |
|-----------|---------------|-----------------|
| --output_folder   |         | metagenomics-nf_results |
| --cpu    |            2 | Number of CPUs |
| --mem    |            8 | Memory in Gb |
|  --ref       |    | Host reference genome   |
| --virusbreakend_db | |  Virusbreakend database (e.g., virusbreakenddb_20210401 from https://github.com/PapenfussLab/gridss/blob/master/VIRUSBreakend_Readme.md) | 

  * #### Flags

Flags are special parameters without value. The virusbreakend process will be triggered automatically if the virusbreakend_db parameter is specified.

| Name      | Description     |
|-----------|-----------------|
| --help    | Display help |


## Usage
  ```
  nextflow run iarcbioinfo/metagenomics-nf --input_folder crams
  ```

## Output
  | Type      | Description     |
  |-----------|---------------|
  | *_centrifuge_report.tsv   | Centrifuge summary reports per sample |
  | *_centrifuge_results.txt    | Centrifuge output per sample |
  | *virusbreakend.vcf | Virusbreakend vcf file output |




## Contributions

  | Name      | Email | Description     |
  |-----------|---------------|-----------------|
  | Nicolas Alcala*    |       alcalan@iarc.who.int | Developer to contact for support |


## References 
[Kim D, Song L, Breitwieser FP, and Salzberg SL. Centrifuge: rapid and sensitive classification of metagenomic sequences. Genome Research 2016](https://genome.cshlp.org/content/early/2016/11/16/gr.210641.116.abstract)
