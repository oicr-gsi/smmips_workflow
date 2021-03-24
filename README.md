# smmipsWorkflow

Analysis of smMIP libraries

## Overview

## Dependencies

* [bwa 0.7.12](http://bio-bwa.sourceforge.net/)
* [python 3.6](https://www.python.org/downloads/)
* [smmips 1.0.0](https://pypi.org/project/smmips/)


## Usage

### Cromwell
```
java -jar cromwell.jar run smmipsWorkflow.wdl --inputs inputs.json
```

### Inputs

#### Required workflow parameters:
Parameter|Value|Description
---|---|---
`fastq1`|File|Path to Fastq1
`fastq2`|File|Path to Fastq2
`panel`|File|Path to file with smMIP information
`outputFileNamePrefix`|String|Prefix used to name the output files
`cosmicFile`|File|Tab separated table of all COSMIC coding point mutations from targeted and genome wide screens


#### Optional workflow parameters:
Parameter|Value|Default|Description
---|---|---|---
`outdir`|String|"./"|Path to directory where directory structure is created
`maxSubs`|Int|0|Maximum number of substitutions allowed in the probe sequence
`upstreamNucleotides`|Int|0|Maximum number of nucleotides upstream the UMI sequence
`umiLength`|Int|4|Length of the UMI
`match`|Int|2|Score of identical characters
`mismatch`|Int|-1|Score of non-identical characters
`gapOpening`|Float|-5|Score for opening a gap
`gapExtension`|Float|-1|Score for extending an open gap
`alignmentOverlapThreshold`|Int|60|Cut-off value for the length of the de-gapped overlap between read1 and read2
`matchesThreshold`|Float|0.7|Cut-off value for the number of matching pos
`remove`|Boolean|false|Remove intermediate files if True
`truncate`|Boolean|false|Only pileup columns in the exact region specificied are returned. Default is False
`stepper`|String|"nofilter"|Filter or include reads in the pileup. See pysam doc for behavior of the all or nofilter options. Default is nofilter
`maxDepth`|Int|1000000|Maximum read depth. Default is 1000000
`ignoreOrphans`|Boolean|false|Ignore orphans (paired reads that are not in a proper pair). Default is False
`referenceName`|String|"37"|Reference genome. Must be the same reference used in panel. Accepted values: 37 or 38


#### Optional task parameters:
Parameter|Value|Default|Description
---|---|---|---
`assignSmmips.modules`|String|"smmips/1.0.0 hg19-bwa-index/0.7.12 bwa/0.7.12"|Names and versions of modules to load
`assignSmmips.memory`|Int|32|Memory allocated for this job
`assignSmmips.timeout`|Int|36|Hours before task timeout
`assignSmmips.refFasta`|String|"$HG19_BWA_INDEX_ROOT/hg19_random.fa"|Path to to the reference genome
`assignSmmips.refFai`|String|"$HG19_BWA_INDEX_ROOT/hg19_random.fa.fai"|Path to the reference index
`assignSmmips.refDict`|String|"$HG19_BWA_INDEX_ROOT/hg19_random.dict"|Path to the reference dictionary
`assignSmmips.bwa`|String|"$BWA_ROOT/bin/bwa"|Path to the bwa script
`countVariants.modules`|String|"smmips/1.0.0"|Names and versions of modules to load
`countVariants.memory`|Int|32|Memory allocated for this job
`countVariants.timeout`|Int|24|Hours before task timeout


### Outputs

Output | Type | Description
---|---|---
`outputExtractionMetrics`|File|Metrics file with extracted read counts
`outputReadCounts`|File|Metric file with read counts for each smmip
`outputSortedbam`|File|Alignment file with all reads
`outputSortedbamIndex`|File|Index of the alignment file with all reads
`outputAssignedBam`|File|Alignment file with assigned reads to smmips
`outputAssignedBamIndex`|File|Index file of the alignment file with assigned reads to smmips
`outputUnassignedBam`|File|Alignment file with unassigned reads
`outputUnassignedBamIndex`|File|Index file of the alignment file with unassigned reds
`outputEmptyBam`|File|Alignment file with empty reads
`outputEmptyBamIndex`|File|Index file of the alignment file with empty reads
`outputCountTable`|File|Table with variant counts


## Niassa + Cromwell

This WDL workflow is wrapped in a Niassa workflow (https://github.com/oicr-gsi/pipedev/tree/master/pipedev-niassa-cromwell-workflow) so that it can used with the Niassa metadata tracking system (https://github.com/oicr-gsi/niassa).

* Building
```
mvn clean install
```

* Testing
```
mvn clean verify \
-Djava_opts="-Xmx1g -XX:+UseG1GC -XX:+UseStringDeduplication" \
-DrunTestThreads=2 \
-DskipITs=false \
-DskipRunITs=false \
-DworkingDirectory=/path/to/tmp/ \
-DschedulingHost=niassa_oozie_host \
-DwebserviceUrl=http://niassa-url:8080 \
-DwebserviceUser=niassa_user \
-DwebservicePassword=niassa_user_password \
-Dcromwell-host=http://cromwell-url:8000
```

## Support

For support, please file an issue on the [Github project](https://github.com/oicr-gsi) or send an email to gsi@oicr.on.ca .

_Generated with generate-markdown-readme (https://github.com/oicr-gsi/gsi-wdl-tools/)_
