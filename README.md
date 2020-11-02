# smmipsWorkflow #

smmips is tool for analysis Single Molecule Molecular Inversion Probes (smMIPs) libraries. It generates a table with variant counts (SNPs and indels) at each position covered in the input panel for each smMIP and UMI.
It can also be used for QC as it outputs read counts for each smMIP and estimates the number of empty smMIPs. 

## Dependencies ##
- [bwa/0.7.12](https://github.com/lh3/bwa)
- [smmips/1.0.0](https://github.com/oicr-gsi/pysmmips)
- [python3.6](https://www.python.org/downloads/)

## Usage ##

### Cromwell ###

```java -jar cromwell.jar run smmipsWorkflow.wdl --inputs inputs.json```

### Inputs ###

Required workflow parameters

| Parameter | Type | Description                                    |  Description                                    |
| ------- | ------- | ------------------------------------------ |  --------------                                    |
| ```fastq1``` | File | Path to Fastq1                |  Description                                    |
| ```fastq2``` | File | Path to Fastq2              |  Description                                    |
| ```panel```  | File | Path to file with smMIP information              |  Description                                    |
| ```prefix``` | String | Prefix used to name the output files              |  Description                                    |
| ```cosmicFile``` | File | Tab separated table of all COSMIC coding point mutations from targeted and genome wide screens              |  Description                                    |

Optional workflow parameters

| Parameter | Type | Default                                     |  Description                                    |
| ------- | ------- | ------------------------------------------ |  -----------                                    |
| ```remove``` | Boolean | false              |  Remove intermediate files                                    |
| ```truncate``` | Boolean | false              |  Only pileup columns in the exact region specificied are returned                                    |
| ```ignoreOrphans``` | Boolean | false              |  Ignore orphans (paired reads that are not in a proper pair)                                    |
| ```maxSubs``` | Int | 0              |  Maximum number of substitutions allowed in the probe sequence                                    |
| ```upstreamNucleotides``` | Int | 0              |  Maximum number of nucleotides upstream the UMI sequence                                    |
| ```umiLength``` | Int | 4              |  Length of the UMI                                    |
| ```match``` | Int | 2              |  Score of identical characters                                    |
| ```mismatch``` | Int | -1              |  Score of non-identical characters                                    |
| ```gapOpening``` | Float | -5              |  Score for opening a gap                                    |
| ```gapExtension``` | Float | -1              |  Score for extending an open gap                                    |
| ```alignmentOverlapThreshold``` | Int | 60              |  Cut-off value for the length of the de-gapped overlap between read1 and read2                                    |
| ```matchesThreshold``` | Float | 0.7              |  Cut-off value for the number of matching positions                                    |
| ```stepper``` | String | nofilter              |  Filter or include reads in the pileup                                    |
| ```maxDepth``` | Int | 1000000              |  Maximum read depth in the pileup                                    |
| ```referenceName``` | String | 37              |  Reference genome. Must be the same reference used in panel. Accepted values: 37 or 38                                    |
| ```outdir``` | String | "./"              |  Path to directory where directory structure is created                                    |


Optional task parameters

| Parameter | Type | Default                                     |  Description                                    |
| ------- | ------- | ------------------------------------------ |  -----------                                    |
| ```assignSmmips.memory``` | Int | 32              |  Memory allocated for job                                    |
| ```assignSmmips.timeout``` | Int | 36              |  Only pileup columns in the exact region specificied are returned                                    |
| ```assignSmmips.modules``` | String | smmips/1.0.0 hg19-bwa-index/0.7.12 bwa/0.7.12              |  Environment module names and version to load (space separated) before command execution                                    |
| ```assignSmmips.refFasta``` | Spring | "$HG19_BWA_INDEX_ROOT/hg19_random.fa"              |  Path to the reference fasta                                    |
| ```assignSmmips.refFai``` | Spring | "$HG19_BWA_INDEX_ROOT/hg19_random.fa.fai"              |  Index of the reference genome                                    |
| ```assignSmmips.refDict``` | Spring | "$HG19_BWA_INDEX_ROOT/hg19_random.dict"              |  Dictionary of the reference genome                                    |
| ```assignSmmips.bwa``` | Spring | "$BWA_ROOT/bin/bwa"              |  Path to the bwa script                                    |
| ```countVariants.modules``` | Spring | "smmips/1.0.0"              |  Environment module names and version to load (space separated) before command execution                                    |
| ```countVariants.memory``` | Int | 32              |  Memory allocated for job                                    |
| ```countVariants.timeout``` | Int | 24              |  Only pileup columns in the exact region specificied are returned                                    |


### Outputs ###


| Output | Type | Description                                    |
| ------- | ------- | ------------------------------------------ |
| ```outputSortedbam``` | File | Alignments of reads containing UMIs in paired input fastqs                |
| ```outputSortedbamIndex``` | File | Index file of aligned reads containing UMIs              |
| ```outputAssignedBam```  | File | Alignment of reads assigned to smMIPs. Reads are tagged with smMIP and UMI              |
| ```outputAssignedBamIndex``` | File | Index file of aligned and assigned reads              |
| ```outputUnassignedBam``` | File | Alignment of reads that cannot be assigned to smMIPs              |
| ```outputUnassignedBamIndex``` | File | Index file of aligned but unassigned reads              |
| ```outputEmptyBam``` | File | Alignment of reads assigned to smMIPs but missing target capture              |
| ```outputEmptyBamIndex``` | File | Index file of reads with empty smMIPs              |
| ```outputCountTable``` | File | Table with variant counts at each smMIP position              |
| ```outputExtractionMetrics``` | File | Json file with read counts              |
| ```outputReadCounts``` | File | Json file with read counts with and without target for each smMIP in the panel              |


## Niassa + Cromwell ##

This WDL workflow is wrapped in a Niassa workflow (https://github.com/oicr-gsi/pipedev/tree/master/pipedev-niassa-cromwell-workflow) so that it can used with the Niassa metadata tracking system (https://github.com/oicr-gsi/niassa).

- Building

```mvn clean install```

- Testing

```mvn clean verify \
-Djava_opts="-Xmx1g -XX:+UseG1GC -XX:+UseStringDeduplication" \
-DrunTestThreads=2 \
-DskipITs=false \
-DskipRunITs=false \
-DworkingDirectory=/path/to/tmp/ \
-DschedulingHost=niassa_oozie_host \
-DwebserviceUrl=http://niassa-url:8080 \
-DwebserviceUser=niassa_user \
-DwebservicePassword=niassa_user_password \
-Dcromwell-host=http://cromwell-url:8000```


## Support ##
For support, please file an issue on the [Github project](https://github.com/oicr-gsi/smmipsWorkflow) or send an email to gsi@oicr.on.ca .

Generated with wdl_doc_gen (https://github.com/oicr-gsi/wdl_doc_gen/)