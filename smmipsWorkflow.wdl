version 1.0

workflow smmipsWorkflow {
  input {
    File fastq1
    File fastq2
    File panel
    String outdir = "./"    
    String outputFileNamePrefix  
    Int maxSubs = 0
    Int upstreamNucleotides = 0
    Int umiLength = 4  
    Int match = 2
    Int mismatch = -1
    Float gapOpening = -5
    Float gapExtension   = -1
    Int alignmentOverlapThreshold = 60
    Float matchesThreshold = 0.7  
    Boolean remove = false
    Boolean truncate = false
    String stepper = "nofilter"
    Int maxDepth = 1000000
    Boolean ignoreOrphans = false
    String referenceName = "37"
    File cosmicFile
  }


  parameter_meta {
    outdir: "Path to directory where directory structure is created"
    fastq1: "Path to Fastq1"
    fastq2: "Path to Fastq2"
    outputFileNamePrefix: "Prefix used to name the output files"
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
    dependencies: [
      {
        name: "bwa/0.7.12",
        url: "http://bio-bwa.sourceforge.net/"
      },
      {
        name: "python/3.6",
        url: "https://www.python.org/downloads/"
      },
      {
        name: "smmips/1.0.0",
        url: https://pypi.org/project/smmips/
      }
    ]
  }


  Boolean removeIntermediate = if remove then true else false

  call assignSmmips {
    input:
      fastq1 = fastq1,
      fastq2 = fastq2,
      panel = panel,
      outdir = outdir,
      outputFileNamePrefix = outputFileNamePrefix,  
      maxSubs = maxSubs,
      upstreamNucleotides = upstreamNucleotides,
      umiLength = umiLength, 
      match = match,
      mismatch = mismatch,
      gapOpening = gapOpening,
      gapExtension = gapExtension,  
      alignmentOverlapThreshold = alignmentOverlapThreshold,
      matchesThreshold = matchesThreshold,
      remove = removeIntermediate
  }

  File assignedBam = assignSmmips.assignedBam 
  File assignedBamIndex = assignSmmips.assignedBamIndex 

  Boolean truncateColumn = if truncate then true else false
  Boolean ignoreOrphanReads = if ignoreOrphans then true else false

  call countVariants {
    input: 
      assignedBam = assignedBam,
      assignedBamIndex = assignedBamIndex,
      panel = panel,
      outdir = outdir,
      outputFileNamePrefix = outputFileNamePrefix,  
      truncate = truncateColumn,
      ignoreOrphans = ignoreOrphanReads,
      stepper = stepper,
      maxDepth = maxDepth,
      referenceName = referenceName,
      cosmicFile = cosmicFile
  }

  output {
    File outputExtractionMetrics = assignSmmips.extractionMetrics
    File outputReadCounts = assignSmmips.readCounts
    File outputSortedbam = assignSmmips.sortedbam
    File outputSortedbamIndex = assignSmmips.sortedbamIndex
    File outputAssignedBam = assignSmmips.assignedBam
    File outputAssignedBamIndex = assignSmmips.assignedBamIndex
    File outputUnassignedBam = assignSmmips.unassignedBam
    File outputUnassignedBamIndex = assignSmmips.unassignedBamIndex
    File outputEmptyBam = assignSmmips.emptyBam
    File outputEmptyBamIndex = assignSmmips.emptyBamIndex
    File outputCountTable = countVariants.countTable 
  }
}


task assignSmmips {
  input {
    String modules = "smmips/1.0.0 hg19-bwa-index/0.7.12 bwa/0.7.12"
    Int memory = 32
    Int timeout = 36
    File fastq1
    File fastq2
    File panel
    String outdir = "./"    
    String outputFileNamePrefix  
    Int maxSubs = 0
    Int upstreamNucleotides = 0
    Int umiLength = 4 
    Int match = 2
    Int mismatch = -1
    Float gapOpening = -5
    Float gapExtension = -1  
    Int alignmentOverlapThreshold = 60
    Float matchesThreshold = 0.7  
    Boolean remove
    String refFasta = "$HG19_BWA_INDEX_ROOT/hg19_random.fa"
    String refFai = "$HG19_BWA_INDEX_ROOT/hg19_random.fa.fai"
    String refDict = "$HG19_BWA_INDEX_ROOT/hg19_random.dict"
    String bwa = "$BWA_ROOT/bin/bwa"
  }

  
    parameter_meta {
      modules: "Names and versions of modules to load"
      memory: "Memory allocated for this job"
      timeout: "Hours before task timeout"
      fastq1: "Path to Fastq1"
      fastq2: "Path to Fastq2"
      panel: "Path to file with smMIP information"
      outdir: "Path to directory where directory structure is created"
      outputFileNamePrefix: "Prefix used to name the output files"
      maxSubs: "Maximum number of substitutions allowed in the probe sequence"
      upstreamNucleotides: "Maximum number of nucleotides upstream the UMI sequence"
      umiLength: "Length of the UMI"
      match: "Score of identical characters"
      mismatch: "Score of non-identical characters"
      gapOpening: "Score for opening a gap"
      gapExtension: "Score for extending an open gap"
      alignmentOverlapThreshold: "Cut-off value for the length of the de-gapped overlap between read1 and read2"
      matchesThreshold: "Cut-off value for the number of matching pos"
      remove: "Remove intermediate files if True"
      refFasta: "Path to to the reference genome"
      refFai: "Path to the reference index"
      refDict: "Path to the reference dictionary"
      bwa: "Path to the bwa script"
  }

  String removeFlag = if remove then "--remove" else ""

  command <<<
    set -euo pipefail
    cp ~{refFai} .
    cp ~{refDict} .
    smmips assign -f1 ~{fastq1} -f2 ~{fastq2} -pa ~{panel} -r ~{refFasta} -pf ~{outputFileNamePrefix} -s ~{maxSubs} -up ~{upstreamNucleotides} -umi ~{umiLength}  -m ~{match} -mm ~{mismatch} -go ~{gapOpening} -ge ~{gapExtension}  -ao ~{alignmentOverlapThreshold} -mt ~{matchesThreshold} -o ~{outdir} -bwa ~{bwa} ~{removeFlag}
  >>>

  runtime {
    memory:  "~{memory} GB"
    modules: "~{modules}"
    timeout: "~{timeout}"
  }

  output {
  File extractionMetrics = "${outdir}/stats/${outputFileNamePrefix}_extraction_metrics.json"
  File readCounts = "${outdir}/stats/${outputFileNamePrefix}_smmip_counts.json"
  File sortedbam = "${outdir}/out/${outputFileNamePrefix}.sorted.bam"
  File sortedbamIndex = "${outdir}/out/${outputFileNamePrefix}.sorted.bam.bai"
  File assignedBam = "${outdir}/out/${outputFileNamePrefix}.assigned_reads.sorted.bam"
  File assignedBamIndex = "${outdir}/out/${outputFileNamePrefix}.assigned_reads.sorted.bam.bai"
  File unassignedBam = "${outdir}/out/${outputFileNamePrefix}.unassigned_reads.sorted.bam"
  File unassignedBamIndex = "${outdir}/out/${outputFileNamePrefix}.unassigned_reads.sorted.bam.bai"
  File emptyBam = "${outdir}/out/${outputFileNamePrefix}.empty_reads.sorted.bam"
  File emptyBamIndex = "${outdir}/out/${outputFileNamePrefix}.empty_reads.sorted.bam.bai"
  }

  meta {
    output_meta: {
      sortedbam: "Alignments of reads containing UMIs in paired input fastqs",
      sortedbamIndex: "Index file of aligned reads containing UMIs",
      assignedBam: "Alignment of reads assigned to smMIPs. Reads are tagged with smMIP and UMI",
      assignedBamIndex: "Index file of aligned and assigned reads",
      unassignedBam: "Alignment of reads that cannot be assigned to smMIPs",
      unassignedBamIndex: "Index file of aligned but unassigned reads",
      emptyBam: "Alignment of reads assigned to smMIPs but missing target capture",
      emptyBamIndex: "Index file of reads with empty smMIPs",
      extractionMetrics: "Json file with read counts",
      readCounts: "Json file with read counts with and without target for each smMIP in the panel"
    }
  }


task countVariants {
  input {
    String modules = "smmips/1.0.0"
    File assignedBam
    File assignedBamIndex
    File panel
    String outdir = "./"
    String outputFileNamePrefix  
    Boolean truncate
    Boolean ignoreOrphans
    String stepper = "nofilter"
    Int maxDepth = 1000000
    String referenceName = "37"
    File cosmicFile
    Int memory = 32
    Int timeout = 24
  }


  parameter_meta {
    modules: "Names and versions of modules to load"
    assignedBam: "Bam with UMI-ammotated reads" 
    panel: "Path to file with smMIP information"
    outdir: "Path to directory where directory structure is created"
    outputFileNamePrefix: "Prefix used to name the output files"
    truncate: "Only pileup columns in the exact region specificied are returned. Default is False"
    ignoreOrphans: "Ignore orphans (paired reads that are not in a proper pair). Default is False"
    stepper: "Filter or include reads in the pileup. See pysam doc for behavior of the all or nofilter options. Default is nofilter"
    maxDepth: "Maximum read depth. Default is 1000000"
    referenceName: "Reference genome. Must be the same reference used in panel. Accepted values: 37 or 38"
    cosmicFile: "Tab separated table of all COSMIC coding point mutations from targeted and genome wide screens"
    memory: "Memory allocated for this job"
    timeout: "Hours before task timeout"
  }


  String truncateFlag = if truncate then "-t" else ""  
  String ignoreOrphansFlag = if ignoreOrphans then "-io" else "" 

  command <<<
    set -euo pipefail
    cp ~{assignedBamIndex} .
    smmips variant -b ~{assignedBam} -p ~{panel} -m ~{maxDepth} -o ~{outdir} -stp ~{stepper} -pf ~{outputFileNamePrefix} -rf ~{referenceName} -c ~{cosmicFile} ~{ignoreOrphansFlag} ~{truncateFlag}
  >>>

  runtime {
    memory:  "~{memory} GB"
    modules: "~{modules}"
    timeout: "~{timeout}"
  }

  output {
  File countTable = "${outdir}/out/${outputFileNamePrefix}_Variant_Counts.txt"
  }

  meta {
    output_meta: {
      CountTable: "Table with variant counts at each smMIP position",
    }
  }
}

