version 1.0

workflow smmipsQC {
  input {
    File fastq1
    File fastq2
    File panel
    File smmipRegions
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
  }


  parameter_meta {
    outdir: "Path to directory where directory structure is created"
    fastq1: "Path to Fastq1"
    fastq2: "Path to Fastq2"
    outputFileNamePrefix: "Prefix used to name the output files"
    remove: "Remove intermediate files if True"
    panel: "Path to file with smMIP information"
    smmipRegions: "Path to bed file with smmip regions"
    upstreamNucleotides: "Maximum number of nucleotides upstream the UMI sequence"
    umiLength: "Length of the UMI"
    maxSubs: "Maximum number of substitutions allowed in the probe sequence"
    match: "Score of identical characters"
    mismatch: "Score of non-identical characters"
    gapOpening: "Score for opening a gap"
    gapExtension: "Score for extending an open gap"
    alignmentOverlapThreshold: "Cut-off value for the length of the de-gapped overlap between read1 and read2"
    matchesThreshold: "Cut-off value for the number of matching pos"
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
        name: "smmips/1.0.9",
        url: "https://pypi.org/project/smmips/"
      }
    ]
    
    output_meta: {
    outputExtractionMetrics: "Metrics file with extracted read counts",
    outputReadCounts: "Metric file with read counts for each smmip"
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


  call regionsToArray {
    input:
      regions = smmipRegions
  }

  Array[String] genomic_regions = regionsToArray.out

  scatter(region in genomic_regions) {
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
        region = region
    }
  }
  
  
  Array[File] extractionCounts = assignSmmips.extractionMetrics
  Array[File] readCounts = assignSmmips.readCounts
  
  
  call mergeExtraction {
    input:
      outdir = outdir,    
      remove = removeIntermediate,
      outputFileNamePrefix = outputFileNamePrefix,
      extractionCounts = extractionCounts
  }

  call mergeCounts {
    input:
      outdir = outdir,    
      remove = removeIntermediate,
      outputFileNamePrefix = outputFileNamePrefix,
      readCounts = readCounts
  }

  output {
    File outputExtractionMetrics = mergeExtraction.mergedExtractionMetrics
    File outputReadCounts = mergeCounts.mergedReadCounts
  }
}


task assignSmmips {
  input {
    String modules = "smmips/1.0.9"
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
    String region
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
    region: "Genomic region with grouped smmips"
  }

  String removeFlag = if remove then "--remove" else ""

  command <<<
    set -euo pipefail
    smmips assign -b ~{sortedbam} -pa ~{panel} -pf ~{outputFileNamePrefix} -ms ~{maxSubs} -up ~{upstreamNucleotides} -umi ~{umiLength}  -m ~{match} -mm ~{mismatch} -go ~{gapOpening} -ge ~{gapExtension}  -ao ~{alignmentOverlapThreshold} -mt ~{matchesThreshold} -o ~{outdir} -r ~{region} ~{removeFlag}
  >>>

  runtime {
    memory:  "~{memory} GB"
    modules: "~{modules}"
    timeout: "~{timeout}"
  }

  output {
  File extractionMetrics = "${outdir}/stats/${outputFileNamePrefix}_temp.${region}.extraction_metrics.json"
  File readCounts = "${outdir}/stats/${outputFileNamePrefix}_temp.${region}.smmip_counts.json"
  }

  meta {
    output_meta: {
      extractionMetrics: "Json file with read counts",
      readCounts: "Json file with read counts with and without target for each smMIP in the panel"
    }
  }
}



task align {
  input {
    String modules = "smmips/1.0.9 hg19-bwa-index/0.7.12 bwa/0.7.12"
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



task mergeExtraction {
  input {
    String modules = "smmips/1.0.9"
    Int memory = 32
    Int timeout = 36
    Boolean remove
    String outputFileNamePrefix
    Array[File] extractionCounts
    String outdir = "./"  
  }

  
  parameter_meta {
    modules: "Names and versions of modules to load"
    memory: "Memory allocated for this job"
    timeout: "Hours before task timeout"
    remove: "Remove intermediate files if True"
    outputFileNamePrefix: "Prefix used to name the output files"
    extractionCounts: "List of files to be merged"
    outdir: "Path to output directory"
  }

  String removeFlag = if remove then "--remove" else ""
  String outputName = basename(outputFileNamePrefix)

  command <<<
    set -euo pipefail
    smmips merge -pf ~{outputFileNamePrefix} -ft extraction -t ~{sep =" " extractionCounts}  ~{removeFlag}
  >>>

  runtime {
    memory:  "~{memory} GB"
    modules: "~{modules}"
    timeout: "~{timeout}"
  }

  output {
  File mergedExtractionMetrics = "${outputName}_extraction_metrics.json"
  }

  meta {
    output_meta: {
      mergedExtractionMetrics: "Metrics file with extracted read counts"
    }
  }
}



task mergeCounts {
  input {
    String modules = "smmips/1.0.9"
    Int memory = 32
    Int timeout = 36
    Boolean remove
    String outputFileNamePrefix
    Array[File] readCounts
    String outdir = "./"  
  }

  
  parameter_meta {
    modules: "Names and versions of modules to load"
    memory: "Memory allocated for this job"
    timeout: "Hours before task timeout"
    remove: "Remove intermediate files if True"
    outputFileNamePrefix: "Prefix used to name the output files"
    readCounts: "List of files to be merged"
    outdir: "Path to outout directory"
  }

  String removeFlag = if remove then "--remove" else ""
  String outputName = basename(outputFileNamePrefix)

  command <<<
    set -euo pipefail
    smmips merge -pf ~{outputFileNamePrefix} -ft counts -t ~{sep =" " readCounts} ~{removeFlag}
  >>>

  runtime {
    memory:  "~{memory} GB"
    modules: "~{modules}"
    timeout: "~{timeout}"
  }

  output {
  File mergedReadCounts = "${outputName}_smmip_counts.json"
  }

  meta {
    output_meta: {
      mergedReadCounts: "Metric file with read counts for each smmip"
    }
  }
}


task regionsToArray {
  input {
    File regions
    Int memory = 1
    Int timeout = 1
  }

  command <<<
    cat ~{regions} | sed 's/\t/./g'
  >>>

  output {
    Array[String] out = read_lines(stdout())
  }

  runtime {
    memory: "~{memory} GB"
    timeout: "~{timeout}"
  }

  parameter_meta {
    regions: "Bed file with smmip region coordinates"
    memory: "Memory allocated for this job"
    timeout: "Hours before task timeout"
  }
}