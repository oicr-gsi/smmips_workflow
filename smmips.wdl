version 1.0

workflow smmips {
  input {
    File fastq1
    File fastq2
    File panel
    File reference
    String? outdir    
    String bwa    
    String prefix  
    Int? max_subs
    Int? upstream_nucleotides
    Int? umi_length  
    Int? match
    Int? mismatch
    Float? gap_opening
    Float? gap_extension  
    Int? alignment_overlap_threshold
    Float? matches_threshold  
    Boolean? remove
    Boolean? truncate
    String? stepper
    Int? max_depth
    Boolean? ignore_orphans
    String reference_name
    File cosmicfile
}


  parameter_meta {
    outdir: "Path to directory where directory structure is created"
    fastq1: "Path to Fastq1"
    fastq2: "Path to Fastq2"
    reference: "Path to the reference genome"
    prefix: "Prefix used to name the output files"
    remove: "Remove intermediate files if True"
    panel: "Path to file with smMIP information"
    upstream_nucleotides: "Maximum number of nucleotides upstream the UMI sequence"
    umi_length: "Length of the UMI"
    max_subs: "Maximum number of substitutions allowed in the probe sequence"
    match: "Score of identical characters"
    mismatch: "Score of non-identical characters"
    gap_opening: "Score for opening a gap"
    gap_extension: "Score for extending an open gap"
    alignment_overlap_threshold: "Cut-off value for the length of the de-gapped overlap between read1 and read2"
    matches_threshold: "Cut-off value for the number of matching pos"
    bwa: "Path to the bwa script"
    max_depth: "Maximum read depth. Default is 1000000"
    ignore_orphans: "Ignore orphans (paired reads that are not in a proper pair). Default is False"
    truncate: "Only pileup columns in the exact region specificied are returned. Default is False"
    stepper: "Filter or include reads in the pileup. See pysam doc for behavior of the all or nofilter options. Default is nofilter"
    reference_name: "Reference genome. Must be the same reference used in panel. Accepted values: 37 or 38"
    cosmicfile: "Tab separated table of all COSMIC coding point mutations from targeted and genome wide screens"
}

  meta {
    author: "Richard Jovelin"
    email: "richard.jovelin@oicr.on.ca"
    description: "Analysis of smMIP libraries"
  }

  call assign_smmips {
    input:
      fastq1 = fastq1,
      fastq2 = fastq2,
      panel = panel,
      reference = reference,
      outdir = outdir,
      bwa = bwa,    
      prefix = prefix,  
      max_subs = max_subs,
      upstream_nucleotides = upstream_nucleotides,
      umi_length = umi_length, 
      match = match,
      mismatch = mismatch,
      gap_opening = gap_opening,
      gap_extension = gap_extension,  
      alignment_overlap_threshold = alignment_overlap_threshold,
      matches_threshold = matches_threshold,  
      remove = remove
  }

  File sortedbam = assign_smmips.sortedbam 

  call count_variants {
    input: 
      sortedbam = sortedbam,
      panel = panel,
      outdir = outdir,
      prefix = prefix,  
      truncate = truncate,
      ignore_orphans = ignore_orphans,
      stepper = stepper,
      max_depth = max_depth,
      reference_name = reference_name,
      cosmicfile = cosmicfile
  }

  File count_table = count_variants.count_table
  
  output {
    Array[File] stats_files = assign_smmips.stats_files
    File sortedbam = assign_smmips.sortedbam
    File sortedbam_index = assign_smmips.sortedbam_index
    File assigned_bam = assign_smmips.assigned_bam
    File assigned_bam_index = assign_smmips.assigned_bam_index
    File unassigned_bam = assigned_smmips.unassigned_bam
    File unassigned_bam_index = assigned_smmips.unassigned_bam_index
    File empty_bam = assigned_smmips.empty_bam
    File empty_bam_index = assigned_smmips.empty_bam_index
    File count_table = count_variants.count_table 
  }
}


task count_variants {
  input {
    String modules = "smmips/1.0.0"
    File sortedbam
    File panel
    String? outdir
    String prefix  
    Boolean? truncate
    Boolean? ignore_orphans
    String? stepper
    Int? max_depth
    String reference_name
    File cosmicfile
    Int memory = 32
  }

  command <<<
    smmips variant -b ~{sortedbam} -p ~{panel} -m ${default=1000000 max_depth} \
    ${if outdir then "-o ~{outdir}" else ""} -stp ${default='nofilter' stepper} \
    -pf ~{prefix} -rf ~{reference_name} -c ~{cosmicfile} \
    ${if ignore_orphans then "-io" else ""} ${if truncate then "-t" else ""}
  >>>

  runtime {
    memory:  "~{memory} GB"
    modules: "~{modules}"
  }

  output {
  File count_table = "${outdir}/out/${prefix}_Variant_Counts.txt"
  }
}


task assign_smmips {
  input {
    String modules = "smmips/1.0.0"
    Int memory = 32
    File fastq1
    File fastq2
    File panel
    File reference
    String? outdir    
    String bwa    
    String prefix  
    Int? max_subs
    Int? upstream_nucleotides
    Int? umi_length  
    Int? match
    Int? mismatch
    Float? gap_opening
    Float? gap_extension  
    Int? alignment_overlap_threshold
    Float? matches_threshold  
    Boolean? remove
  }

  command <<<
    smmips assign -f1 ~{fastq1} -f2 ~{fastq2} -pa ~{panel} -r ~{reference} \
    -pf ~{prefix} -s ${default=0 max_subs} -up ${default=0 upstream_nucleotides} -umi ${default=4 umi_length} \ 
    -m ${default=2 match} -mm ${default=-1 mismatch} -go ${default=-5 gap_opening} -ge ${default=-1 gap_extension} \
    -ao ${default=60 alignment_overlap} -mt ${default=0.7 matches_threshold} ${if remove then "--remove" else ""} \
    ${if outdir then "-o ~{outdir}" else ""} -bwa ~{bwa}
  >>>

  runtime {
    memory:  "~{memory} GB"
    modules: "~{modules}"
  }

  output {
  Array[File] stats_files = glob("${outdir}/stats/*.json")
  File sortedbam = "${outdir}/out/${prefix}.sorted.bam"
  File sortedbam_index = "${outdir}/out/${prefix}.sorted.bam.bai"
  File assigned_bam = "${outdir}/out/${prefix}.assigned_reads.bam"
  File assigned_bam_index = "${outdir}/out/${prefix}.assigned_reads.bam.bai"
  File unassigned_bam = "${outdir}/out/${prefix}.unassigned_reads.bam"
  File unassigned_bam_index = "${outdir}/out/${prefix}.unassigned_reads.bam.bai"
  File empty_bam = "${outdir}/out/${prefix}.empty_reads.bam"
  File empty_bam_index = "${outdir}/out/${prefix}.empty_reads.bam.bai"
  }
}


task count_variants {
  input {
    String modules = "smmips/1.0.0"
    File sortedbam
    File panel
    String? outdir
    String prefix  
    Boolean? truncate
    Boolean? ignore_orphans
    String? stepper
    Int? max_depth
    String reference_name
    File cosmicfile
    Int memory = 32
  }

  command <<<
    smmips variant -b ~{sortedbam} -p ~{panel} -m ${default=1000000 max_depth} \
    ${if outdir then "-o ~{outdir}" else ""} -stp ${default='nofilter' stepper} \
    -pf ~{prefix} -rf ~{reference_name} -c ~{cosmicfile} \
    ${if ignore_orphans then "-io" else ""} ${if truncate then "-t" else ""}
  >>>

  runtime {
    memory:  "~{memory} GB"
    modules: "~{modules}"
  }

  output {
  File count_table = "${outdir}/out/${prefix}_Variant_Counts.txt"
  }
}

