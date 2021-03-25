version 1.0

workflow smmipsQC {
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
    Int distance
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
    distance: "Minimum distance between smmips for grouping"
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
        name: "smmips/1.0.3",
        url: "https://pypi.org/project/smmips/"
      },
      {
        name: "smmip-region-finder/1.0",
        url: "https://github.com/oicr-gsi/smmipRegionFinder"
      }
    ]
    
    output_meta: {
    outputExtractionMetrics: "Metrics file with extracted read counts",
    outputReadCounts: "Metric file with read counts for each smmip",
    }
  }


  Boolean removeIntermediate = if remove then true else false

  call align {
    input:
      fastq1 = fastq1,
      fastq2 = fastq2,
      outdir = outdir,
      outputFileNamePrefix = outputFileNamePrefix,
      remove = removeIntermediate
  }

   File sortedbam = align.sortedbam
   File sortedbamIndex = align.sortedbamIndex

  call assignSmmips {
    input:
      sortedbam = sortedbam,
      sortedbamIndex = sortedbamIndex,
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
      remove = removeIntermediate,
      chromosome = chromosome,
      start = start,
      end = end
  }

  output {
    File outputExtractionMetrics = assignSmmips.extractionMetrics
    File outputReadCounts = assignSmmips.readCounts
  }
}


task assignSmmips {
  input {
    String modules = "smmips/1.0.3"
    Int memory = 32
    Int timeout = 36
    File sortedbam
    File sortedbamIndex
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
    String chromosome
    Int start
    Int end
  }

  
  parameter_meta {
    modules: "Names and versions of modules to load"
    memory: "Memory allocated for this job"
    timeout: "Hours before task timeout"
    sortedbam: "Alignments of reads containing UMIs in paired input fastqs"
    sortedbamIndex: "Index file of aligned reads containing UMIs"
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
    chromosome: "Chromomosome where reads are mapped"
    start: "Start position of region on chromosome"
    end: "End position of region on chromosome"
  }

  String removeFlag = if remove then "--remove" else ""

  command <<<
    set -euo pipefail
    smmips assign -b ~{sortedbam} -pa ~{panel} -pf ~{outputFileNamePrefix} -ms ~{maxSubs} -up ~{upstreamNucleotides} -umi ~{umiLength}  -m ~{match} -mm ~{mismatch} -go ~{gapOpening} -ge ~{gapExtension}  -ao ~{alignmentOverlapThreshold} -mt ~{matchesThreshold} -o ~{outdir} -c ~{chromosome} -s ~{start} -e ~{end} ~{removeFlag}
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
}



task align {
  input {
    String modules = "smmips/1.0.3 hg19-bwa-index/0.7.12 bwa/0.7.12"
    Int memory = 32
    Int timeout = 36
    File fastq1
    File fastq2
    String outdir = "./"    
    String outputFileNamePrefix  
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
    outdir: "Path to directory where directory structure is created"
    outputFileNamePrefix: "Prefix used to name the output files"
    remove: "Remove intermediate files if True"
    refFasta: "Path to to the reference genome"
    refFai: "Path to the reference index"
    refDict: "Path to the reference dictionary"
    bwa: "Path to the bwa script"
  }

  String removeFlag = if remove then "--remove" else ""

  command <<<
    set -euo pipefail
    smmips align -bwa ~{bwa} -f1 ~{fastq1} -f2 ~{fastq2} -r ~{refFasta} -pf ~{outputFileNamePrefix} -o ~{outdir} ~{removeFlag}
  >>>

  runtime {
    memory:  "~{memory} GB"
    modules: "~{modules}"
    timeout: "~{timeout}"
  }

  output {
  File sortedbam = "${outdir}/out/${outputFileNamePrefix}.sorted.bam"
  File sortedbamIndex = "${outdir}/out/${outputFileNamePrefix}.sorted.bam.bai"
  }

  meta {
    output_meta: {
      sortedbam: "Alignments of reads containing UMIs in paired input fastqs",
      sortedbamIndex: "Index file of aligned reads containing UMIs"
    }
  }
}



task merge {
  input {
    String modules = "smmips/1.0.3"
    Int memory = 32
    Int timeout = 36
    String outdir = "./"    
    Boolean remove
    String outputFileNamePrefix
  }

  
  parameter_meta {
    modules: "Names and versions of modules to load"
    memory: "Memory allocated for this job"
    timeout: "Hours before task timeout"
    outdir: "Path to directory where directory structure is created"
    remove: "Remove intermediate files if True"
    outputFileNamePrefix: "Prefix used to name the output files"
  }

  String removeFlag = if remove then "--remove" else ""

  command <<<
    set -euo pipefail
    cp ~{refFai} .
    cp ~{refDict} .
    smmips merge -o ~{outdir} ~{removeFlag}
  >>>

  runtime {
    memory:  "~{memory} GB"
    modules: "~{modules}"
    timeout: "~{timeout}"
  }

  output {
  File extractionMetrics = "${outdir}/stats/${outputFileNamePrefix}_extraction_metrics.json"
  File readCounts = "${outdir}/stats/${outputFileNamePrefix}_smmip_counts.json"
  File assignedBam = "${outdir}/out/${outputFileNamePrefix}.assigned_reads.sorted.bam"
  File assignedBamIndex = "${outdir}/out/${outputFileNamePrefix}.assigned_reads.sorted.bam.bai"
  File emptyBam = "${outdir}/out/${outputFileNamePrefix}.empty_reads.sorted.bam"
  File emptyBamIndex = "${outdir}/out/${outputFileNamePrefix}.empty_reads.sorted.bam.bai"
  }

  meta {
    output_meta: {
      sortedbam: "Alignments of reads containing UMIs in paired input fastqs",
      sortedbamIndex: "Index file of aligned reads containing UMIs"
      extractionMetrics: "Metrics file with extracted read counts",
      readCounts: "Metric file with read counts for each smmip",
      emptyBam: "Alignment file with empty reads",
      emptyBamIndex: "Index file of the alignment file with empty reads",
    }
  }
}



task findRegions {
  input {
    String modules = "smmipRegionFinder/1.0"
    Int memory = 32
    Int timeout = 36
    String panel    
    Int distance
  }

  parameter_meta {
    modules: "Names and versions of modules to load"
    memory: "Memory allocated for this job"
    timeout: "Hours before task timeout"
    panel: "Path to file with smMIP information"
    distance: "Minimum distance between smmips for grouping"
  }

  command <<<
    set -euo pipefail
    smmipRegionFinder -p ~{panel} -d ~{distance}
  >>>

  runtime {
    memory:  "~{memory} GB"
    modules: "~{modules}"
    timeout: "~{timeout}"
  }

  output {
  Array coordinates
  }

  meta {
    output_meta: {
      coordinates: "Array with coordinates of smmip regions"
    }
  }
}


