version 1.0

workflow smmipsWorkflow {
  input {
    File fastq1
    File fastq2
    File panel
    File reference
    String? outdir    
    String bwa    
    String prefix  
    Int maxSubs = 0
    Int upstreamNucleotides = 0
    Int umiLength = 4  
    Int match = 2
    Int mismatch = -1
    Float gapOpening = -5
    Float gapExtension   = -1
    Int alignmentOverlapThreshold = 60
    Float matchesThreshold = 0.7  
    Boolean? remove
    Boolean? truncate
    String stepper = "nofilter"
    Int maxDepth = 1000000
    Boolean? ignoreOrphans
    String referenceName
    File cosmicFile
  }


  parameter_meta {
    outdir: "Path to directory where directory structure is created"
    fastq1: "Path to Fastq1"
    fastq2: "Path to Fastq2"
    reference: "Path to the reference genome"
    prefix: "Prefix used to name the output files"
    remove: "Remove intermediate files if True"
    panel: "Path to file with smMIP information"
    upstreamNucleotides: "Maximum number of nucleotides upstream the UMI sequence"
    umiLength: "Length of the UMI"
    maxSubs: "Maximum number of substitutions allowed in the probe sequence"
    match: "Score of identical characters"
    mismatch: "Score of non-identical characters"
    gapOpening: "Score for opening a gap"
    gapExtension: "Score for extending an open gap"
    alignmentOverlapThreshold: "Cut-off value for the length of the de-gapped overlap between read1 and read2"
    matchesThreshold: "Cut-off value for the number of matching pos"
    bwa: "Path to the bwa script"
    maxDepth: "Maximum read depth. Default is 1000000"
    ignoreOrphans: "Ignore orphans (paired reads that are not in a proper pair). Default is False"
    truncate: "Only pileup columns in the exact region specificied are returned. Default is False"
    stepper: "Filter or include reads in the pileup. See pysam doc for behavior of the all or nofilter options. Default is nofilter"
    referenceName: "Reference genome. Must be the same reference used in panel. Accepted values: 37 or 38"
    cosmicFile: "Tab separated table of all COSMIC coding point mutations from targeted and genome wide screens"
  }

  meta {
    author: "Richard Jovelin"
    email: "richard.jovelin@oicr.on.ca"
    description: "Analysis of smMIP libraries"
  }

  call assignSmmips {
    input:
      fastq1 = fastq1,
      fastq2 = fastq2,
      panel = panel,
      reference = reference,
      outdir = outdir,
      bwa = bwa,    
      prefix = prefix,  
      maxSubs = maxSubs,
      upstreamNucleotides = upstreamNucleotides,
      umiLength = umiLength, 
      match = match,
      mismatch = mismatch,
      gapOpening = gapOpening,
      gapExtension = gapExtension,  
      alignmentOverlapThreshold = alignmentOverlapThreshold,
      matchesThreshold = matchesThreshold,  
      remove = remove
  }

  #File sortedbam = assignSmmips.sortedbam 

  call countVariants {
    input: 
      sortedbam = assignSmmips.sortedbam,
      panel = panel,
      outdir = outdir,
      prefix = prefix,  
      truncate = truncate,
      ignoreOrphans = ignoreOrphans,
      stepper = stepper,
      maxDepth = maxDepth,
      referenceName = referenceName,
      cosmicFile = cosmicFile
  }

  output {
    Array[File] statsFiles = assignSmmips.statsFiles
    File sortedbam = assignSmmips.sortedbam
    File sortedbamIndex = assignSmmips.sortedbamIndex
    File assignedBam = assignSmmips.assignedBam
    File assignedBamIndex = assignSmmips.assignedBamIndex
    File unassignedBam = assignedSmmips.unassignedBam
    File unassignedBamIndex = assignedSmmips.unassignedBamIndex
    File emptyBam = assignedSmmips.emptyBam
    File emptyBamIndex = assignedSmmips.emptyBamIndex
    File countTable = countVariants.countTable 
  }
}


task assignSmmips {
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
    Int maxSubs = 0
    Int upstreamNucleotides = 0
    Int umiLength = 4 
    Int match = 2
    Int mismatch = -1
    Float gapOpening = -5
    Float gapExtension = -1  
    Int alignmentOverlapThreshold = 60
    Float matchesThreshold = 0.7  
    Boolean? remove
  }

  command <<<
    smmips assign -f1 ~{fastq1} -f2 ~{fastq2} -pa ~{panel} -r ~{reference} \
    -pf ~{prefix} -s ~{maxSubs} -up ~{upstreamNucleotides} -umi ~{umiLength} \ 
    -m ~{match} -mm ~{mismatch} -go ~{gapOpening} -ge ~{gapExtension} \
    -ao ~{alignmentOverlapThreshold} -mt ~{matchesThreshold} ${if remove then "--remove" else ""} \
    ${if outdir then "-o ~{outdir}" else ""} -bwa ~{bwa}
  >>>

  runtime {
    memory:  "~{memory} GB"
    modules: "~{modules}"
  }

  output {
  Array[File] statsFiles = glob("${outdir}/stats/*.json")
  File sortedbam = "${outdir}/out/${prefix}.sorted.bam"
  File sortedbamIndex = "${outdir}/out/${prefix}.sorted.bam.bai"
  File assignedBam = "${outdir}/out/${prefix}.assigned_reads.bam"
  File assignedBamIndex = "${outdir}/out/${prefix}.assigned_reads.bam.bai"
  File unassignedBam = "${outdir}/out/${prefix}.unassigned_reads.bam"
  File unassignedBamIndex = "${outdir}/out/${prefix}.unassigned_reads.bam.bai"
  File emptyBam = "${outdir}/out/${prefix}.empty_reads.bam"
  File emptyBamIndex = "${outdir}/out/${prefix}.empty_reads.bam.bai"
  }
}


task countVariants {
  input {
    String modules = "smmips/1.0.0"
    File sortedbam
    File panel
    String? outdir
    String prefix  
    Boolean? truncate
    Boolean? ignoreOrphans
    String stepper = "nofilter"
    Int maxDepth = 1000000
    String referenceName
    File cosmicFile
    Int memory = 32
  }

  command <<<
    smmips variant -b ~{sortedbam} -p ~{panel} -m ~{maxDepth} \
    ${if outdir then "-o ~{outdir}" else ""} -stp ~{stepper} \
    -pf ~{prefix} -rf ~{referenceName} -c ~{cosmicFile} \
    ${if ignoreOrphans then "-io" else ""} ${if truncate then "-t" else ""}
  >>>

  runtime {
    memory:  "~{memory} GB"
    modules: "~{modules}"
  }

  output {
  File countTable = "${outdir}/out/${prefix}_Variant_Counts.txt"
  }
}

