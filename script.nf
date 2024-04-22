#! /usr/bin/env nextflow

// Copyright (C) 2024 IARC/WHO

// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

//help function for the tool
def show_help (){
  log.info IARC_Header()
  log.info tool_header()
    log.info"""
    Usage:

    The typical command for running the pipeline is as follows:

    nextflow run iarcbioinfo/metagenomics-nf -singularity [OPTIONS]

    Mandatory arguments:
      --input_folder         [file] Folder with bam/cram files to be processed
    Optional arguments:
      --output_folder        [string] name of output folder
      --cpu                  [Integer]  Number of CPUs[def:2]
      --mem 		         [Integer] Max memory [def:8Gb]
      --ref               [file] Host reference genome fasta file (e.g., hg38.fa)
      --virusbreakend_db    [folder] Virusbreakend database (e.g., virusbreakenddb_20210401 from https://github.com/PapenfussLab/gridss/blob/master/VIRUSBreakend_Readme.md)
      """.stripIndent()
}

//extract unmapped reads and create fastq files
process bam2fq {
    cpus params.cpu
    memory params.mem+'G'
    tag { file_tag }
        
	input:
    path(bam)

    output:
    tuple val(file_tag), file("*_unmapped_{1,2}.fastq.gz")

    shell:
    file_tag = bam.baseName
    '''
    samtools collate -uOn 128 !{bam} tmp1_!{file_tag} | samtools fastq -f 4 -F8 --threads !{params.cpu} -1 !{file_tag}_unmapped1_1.fastq.gz -2 !{file_tag}_unmapped1_2.fastq.gz -s !{file_tag}_unmapped1_single.fastq.gz -N -
    samtools collate -uOn 128 !{bam} tmp2_!{file_tag} | samtools fastq -f 8 -F4 --threads !{params.cpu} -1 !{file_tag}_unmapped2_1.fastq.gz -2 !{file_tag}_unmapped2_2.fastq.gz -s !{file_tag}_unmapped2_single.fastq.gz -N -
    samtools collate -uOn 128 !{bam} tmp3_!{file_tag} | samtools fastq -f 12 --threads !{params.cpu} -1 !{file_tag}_unmappedpair_1.fastq.gz -2 !{file_tag}_unmappedpair_2.fastq.gz -s !{file_tag}_unmappedpair_single.fastq.gz -N -
    rm *single.fastq.gz
    cat !{file_tag}_*_1.fastq.gz > !{file_tag}_unmapped_1.fastq.gz 
    cat !{file_tag}_*_2.fastq.gz > !{file_tag}_unmapped_2.fastq.gz
    rm !{file_tag}_unmapped[12]_[12].fastq.gz
    rm !{file_tag}_unmappedpair_[12].fastq.gz
    '''
} //samtools view -h -u -f 12 !{bam} | samtools fastq -1 !{file_tag}_1.fastq.gz -2 !{file_tag}_2.fastq.gz -


//run centrifuge
process centrifuge {
    cpus params.cpu
    memory params.mem+'G'
    tag { ID }
        
	input:
    tuple val(ID), path(fastqpairs)

    output:
    path("${ID}_centrifuge*")
    publishDir params.output_folder, mode: 'copy'

    shell:
    '''
    centrifuge -x p_compressed+h+v -1 !{fastqpairs[0]} -2 !{fastqpairs[1]} -t -p !{params.cpu} --report-file !{ID}_centrifuge_report.tsv -S !{ID}_centrifuge_results.txt
    '''
}

//run viral integration detection with virusbreakend
process virusbreakend {
    cpus params.cpu
    memory params.mem+'G'
    tag { ID }
        
	input:
    path(bam)
    path(host_reference)

    output:
    path("${ID}_virusbreakend.vcf*")
    path("${ID}_tmp/*/*.kraken2.report.all.txt")
    publishDir params.output_folder+"/virusbreakend/", mode: 'copy'

    shell:
    ID = bam.baseName
    '''
    virusbreakend -t !{params.cpu} --gridssargs "--skipsoftcliprealignment" -w !{ID}_tmp \
    -r !{host_reference} \
    -o !{ID}_virusbreakend.vcf \
    --db !{params.virusbreakend_db} \
    !{bam}
    '''
}


// DSL2 workflow to run the processes
workflow{
    params.cpu = 2
    params.mem  = 8
    params.output_folder = "metagenomics-nf_results"
  //display help information
  if (params.help){ show_help(); exit 0;}
  //display the header of the tool
  log.info IARC_Header()
  log.info tool_header()
  //Check mandatory parameters
  assert (params.input_folder != null) : "please specify --input_folder"
  

  //channel with bam files
  bams = Channel.fromPath([params.input_folder+"/*.bam",params.input_folder+"/*.cram"])

  print_params()

  //run centrifuge
  bam2fq(bams) | centrifuge

  //if virusbreakend params given, run virusbreakend
  if(params.virusbreakend_db){
    virusbreakend(bams,params.ref)//,params.virusbreakend_db)
  }
}


// print the calling parameter to the log and a log file
def print_params () {
  //software versions
  def software_versions = ['centrifuge'   : '1.0.4']
  //we print the parameters
  log.info "\n"
  log.info "-\033[2m------------------Calling PARAMETERS--------------------\033[0m-"
  log.info params.collect{ k,v -> "${k.padRight(18)}: $v"}.join("\n")
  log.info "-\033[2m--------------------------------------------------------\033[0m-"
  log.info "\n"
  log.info "-\033[2m------------------Software versions--------------------\033[0m-"
  log.info software_versions.collect{ k,v -> "${k.padRight(18)}: $v"}.join("\n")
  log.info "-\033[2m--------------------------------------------------------\033[0m-"
  log.info "\n"


  //we print the parameters to a log file
   def output_d = new File("${params.output_folder}/nf-pipeline_info/")
   if (!output_d.exists()) {
       output_d.mkdirs()
   }
   def output_tf = new File(output_d, "run_parameters_report.txt")
   def  report_params="------------------Calling PARAMETERS--------------------\n"
        report_params+= params.collect{ k,v -> "${k.padRight(18)}: $v"}.join("\n")
        report_params+="\n--------------------------------------------------------\n"
        report_params+="\n------------------NEXTFLOW Metadata--------------------\n"
        report_params+="nextflow version : "+nextflow.version+"\n"
        report_params+="nextflow build   : "+nextflow.build+"\n"
        report_params+="Command line     : \n"+workflow.commandLine.split(" ").join(" \\\n")
        report_params+="\n--------------------------------------------------------\n"
        report_params+="-----------------Software versions--------------------\n"
        report_params+=software_versions.collect{ k,v -> "${k.padRight(18)}: $v"}.join("\n")
        report_params+="\n--------------------------------------------------------\n"

   output_tf.withWriter { w -> w << report_params}
}


//this use ANSI colors to make a short tool description
//useful url: http://www.lihaoyi.com/post/BuildyourownCommandLinewithANSIescapecodes.html
def tool_header (){
        return """
        metagenomics-nf: Pipeline to perform metagenomic analyses from next-generation sequencing data (${workflow.manifest.version})
        """
}

//header for the IARC tools
// the logo was generated using the following page
// http://patorjk.com/software/taag  (ANSI logo generator)
def IARC_Header (){
     return  """
#################################################################################
# ██╗ █████╗ ██████╗  ██████╗██████╗ ██╗ ██████╗ ██╗███╗   ██╗███████╗ ██████╗  #
# ██║██╔══██╗██╔══██╗██╔════╝██╔══██╗██║██╔═══██╗██║████╗  ██║██╔════╝██╔═══██╗ #
# ██║███████║██████╔╝██║     ██████╔╝██║██║   ██║██║██╔██╗ ██║█████╗  ██║   ██║ #
# ██║██╔══██║██╔══██╗██║     ██╔══██╗██║██║   ██║██║██║╚██╗██║██╔══╝  ██║   ██║ #
# ██║██║  ██║██║  ██║╚██████╗██████╔╝██║╚██████╔╝██║██║ ╚████║██║     ╚██████╔╝ #
# ╚═╝╚═╝  ╚═╝╚═╝  ╚═╝ ╚═════╝╚═════╝ ╚═╝ ╚═════╝ ╚═╝╚═╝  ╚═══╝╚═╝      ╚═════╝  #
# Nextflow pipelines for cancer genomics.########################################
"""
}
