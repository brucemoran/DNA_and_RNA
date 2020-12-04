#!/usr/bin/env nextflow

params.help = ""

if (params.help) {
  log.info ''
  log.info '------------------------------------------------------------'
  log.info 'NEXTFLOW 19.10 DNA AND RNA QC, TRIM, ALIGN, MUTATION CALLING'
  log.info '------------------------------------------------------------'
  log.info ''
  log.info 'Purpose: '
  log.info 'Analyse matched RNA- and DNA-seq, determine extent of same SNP, INDELs'
  log.info ''
  log.info 'Usage: '
  log.info 'nextflow run brucemoran/dna_and_rna'
  log.info ''
  log.info 'Mandatory arguments:'
  log.info '    -profile    Configuration profile (required: standard,singularity)'
  log.info '    --sampleCsv      STRING      CSV format, headers: type (either "DNA" or "RNA"),sampleID,/path/to/read1.fq.gz,/path/to/read2.fq.gz'
  log.info '    --refDir        STRING      dir in which reference data and required indices are held; if not specified, run refs.nf'
  log.info ''
  log.info 'Optional arguments:'
  log.info '    --outDir        STRING        path/to/output'
  log.info '    --multiqcConfig      STRING      config file for multiqc (default: bin/dna_and_rna.multiQC_config.yaml)'
  log.info '    --intList      STRING      path/to/interval.list in chr:start-end for interval calling (default: standardly made reference interval list for species)'
  log.info ''
  exit 1
}

// -2 Test if refDir is defined, if not run DNAseq_references pipeline under defaults
if(!params.refDir){
  exit 1, "Please run: nextflow run refs.nf --outDir work -profile standard,singularity, then specify: nextflow run brucemoran/dna_and_rna --refDir work/refs"
}

//Java task memory allocation via task.memory
javaTaskmem = { it.replace(" GB", "g") }

//Reference data as value channels and reusable therefore
ref = [
    fa: false,
    fai: false,
    dict: false,
    bwa: false,
    star: false,
    saf: false,
    intlist: false,
    refflat: false,
    rrna: false,
    vcf: false,
    vep: false,
    exome: false
]

ref.fa = Channel.value(file(params.references['genome'].fa))
ref.fai = Channel.value(file(params.references['genome'].fai))
ref.dict = Channel.value(file(params.references['genome'].dict))
ref.bwa = Channel.value(file(params.references['genome'].bwa))
ref.star_base = Channel.value(file(params.references['genome'].star_base))
ref.saf = Channel.value(file(params.references['genome'].saf))
ref.intlist = Channel.value(file(params.references['genome'].intlist))
ref.refflat = Channel.value(file(params.references['genome'].refflat))
ref.vep = Channel.value(file(params.references['genome'].vep))
ref.exome = Channel.value(file(params.references['genome'].exome))

//setting of intlist based on exome extant
ref.intlist = ref.exome == null ? ref.inlist : ref.exome

// 0.00: Input using sample.csv
Channel.fromPath("$params.sampleCsv")
       .splitCsv( header: true )
       .map { row -> [row.type, row.sampleID, file(row.read1, checkIfExists: true), file(row.read2, checkIfExists: true)] }
       .set { bbduk_in }

// 0.01: Input trimming
process bbduk {

  label 'medmem'

  publishDir path: "${params.outDir}/DNA_and_RNA/samples/$sampleID/bbduk", mode: "copy", pattern: "*.txt"

  input:
  tuple val(type), val(sampleID), file(read1), file(read2) from bbduk_in

  output:
  file('*.txt') into log_bbduk
  tuple val(type), val(sampleID), file("${sampleID}.bbduk.R1.fastq.gz"), file("${sampleID}.bbduk.R2.fastq.gz") into ( bwa_mem_in, star_in )
  tuple val(type), val(sampleID), file("${sampleID}.bbduk.R1.fastq.gz"), file("${sampleID}.bbduk.R2.fastq.gz"), file(read1), file(read2) into fastp_in

  script:
  taskmem = javaTaskmem("${task.memory}")
  """
  {
  ADAPTERS=\$(readlink -e /opt/miniconda/envs/dna_and_rna/opt/bbmap*/resources/adapters.fa)
  sh bbduk.sh -Xmx$taskmem \
    in1=$read1 \
    in2=$read2 \
    out1=$sampleID".bbduk.R1.fastq.gz" \
    out2=$sampleID".bbduk.R2.fastq.gz" \
    k=31 \
    mink=5 \
    hdist=1 \
    ktrim=r \
    trimq=20 \
    qtrim=rl \
    maq=20 \
    ref=\$ADAPTERS \
    tpe \
    tbo \
    stats=$sampleID".bbduk.adapterstats.txt" \
    overwrite=T
  } 2>&1 | tee > ${sampleID}.bbduk.runstats.txt
  """
}

/* 0.2: fastp QC of pre-, post-bbduk
*/
process fastp {

  label 'lowmem'
  publishDir "${params.outDir}/samples/$sampleID/fastp", mode: "copy", pattern: "*.html"

  input:
  tuple val(type), val(sampleID), file(preread1), file(preread2), file(postread1), file(postread2) from fastp_in

  output:
  file('*.html') into fastp_html
  file('*.json') into fastp_multiqc

  script:
  """
  fastp -w ${task.cpus} -h $sampleID"_pre.fastp.html" -j $sampleID"_pre.fastp.json" --in1 $preread1 --in2 $preread2

  fastp -w ${task.cpus} -h $sampleID"_post.fastp.html" -j $sampleID"_post.fastp.json" --in1 $postread1 --in2 $postread2
  """
}

/* 1.0: DNA alignment
*/
process bwamem {

  label 'highmem'

  input:
  tuple val(type), val(sampleID), file(read1), file(read2) from bwa_mem_in
  file(fa) from ref.fa
  file(bwa) from ref.bwa

  output:
  tuple val(type), val(sampleID), file('*.bam'), file('*.bai') into markdups_bwa_in

  when:
  type == "DNA"

  script:
  """
  DATE=\$(date +"%Y-%m-%dT%T")
  RGLINE="@RG\\tID:$sampleID\\tPL:ILLUMINA\\tSM:$sampleID\\tDS:$type\\tCN:UCD\\tLB:LANE_X\\tDT:\$DATE"

  bwa mem \
    -t${task.cpus} \
    -M \
    -R \$RGLINE \
    ${fa} \
    ${read1} ${read2} | \
    samtools sort -T "tmp."$sampleID -o $sampleID".sort.bam"
  samtools index $sampleID".sort.bam"
  """
}

/* 1.0: RNA alignment
*/
process star {

  label 'highmem'

  publishDir "${params.outDir}/DNA_and_RNA/samples/$sampleID/STAR", mode: "copy", pattern: "${sampleID}.*[!bam,!bai]"

  input:
  tuple val(type), val(sampleID), file(read1), file(read2) from star_in
  file(star_base) from ref.star

  output:
  file('*') into completedstar
  tuple val(type), val(sampleID), file("${sampleID}.Aligned.sortedByCoord.out.bam"), file("${sampleID}.Aligned.sortedByCoord.out.bam.bai") into markdups_star_in
  file("${sampleID}.Log.final.out") into starLOG_multiqc

  when:
  type == "RNA"

  script:
  starDir = "${star_base}/star_*"
  """
  BAMsortRAM=\$(echo ${task.memory} | sed 's/ G/000000000/' | \
    perl -ane '\$o=\$F[0]-1000000000; print "\$o\\n";')
  DATE=\$(date +"%Y-%m-%dT%T")

  STAR --runThreadN ${task.cpus} \
     --genomeDir $starDir \
     --readFilesIn $read1 $read2 \
     --readFilesCommand "zcat" \
     --outFileNamePrefix $sampleID"." \
     --outSAMmode Full \
     --outSAMstrandField intronMotif \
     --outSAMattributes All \
     --outFilterType BySJout \
     --seedSearchStartLmaxOverLread 0.5 \
     --outSAMtype BAM SortedByCoordinate \
     --quantMode TranscriptomeSAM GeneCounts \
     --twopassMode Basic \
     --limitBAMsortRAM \$BAMsortRAM \
     --outSAMattrRGline ID:$sampleID PL:ILLUMINA SM:$sampleID DS:$type CN:UCD LB:LANE_X DT:\$DATE
  samtools index "${sampleID}.Aligned.sortedByCoord.out.bam"
  """
}

/* 1.2: MarkDuplicates
*/
markdups_bwa_in.mix(markdups_star_in).set { markdups_in }
process mrkdup {

  label 'medmem'

  publishDir path: "${params.outDir}/DNA_and_RNA/samples/$sampleID/picard", mode: "copy", pattern: "*.txt"

  input:
  tuple val(type), val(sampleID), file(bam), file(bai) from markdups_in

  output:
  file('*.txt') into mrkdup_multiqc
  tuple val(type), val(sampleID), file('*.md.bam'), file('*.md.bam.bai') into gatk4recal_in

  script:
  taskmem = javaTaskmem("${task.memory}")
  """
  OUTBAM=\$(echo $bam | sed 's/bam/md.bam/')
  OUTMET=\$(echo $bam | sed 's/bam/md.metrics.txt/')
  {
  picard -Xmx$taskmem \
    MarkDuplicates \
    TMP_DIR=./ \
    INPUT=$bam \
    OUTPUT=/dev/stdout \
    COMPRESSION_LEVEL=0 \
    QUIET=TRUE \
    METRICS_FILE=\$OUTMET \
    REMOVE_DUPLICATES=FALSE \
    ASSUME_SORTED=TRUE \
    VALIDATION_STRINGENCY=LENIENT \
    VERBOSITY=ERROR | samtools view -Shb - > \$OUTBAM

  samtools index \$OUTBAM
  } 2>&1 | tee > ${sampleID}.picard_markDuplicates.log.txt
  """
}

/* 1.3: GATK4 BestPractices
*/
process gtkrcl {

  label 'medmem'

  publishDir path: "${params.outDir}/DNA_and_RNA/samples/$sampleID/gatk4", mode: "copy", pattern: "*.GATK4_BQSR.log.txt "

  input:
  tuple val(type), val(sampleID), file(bam), file(bai) from gatk4recal_in
  file(fa) from ref.fa
  file(vcf) from ref.vcf
  file(intlisd) from ref.intlist

  output:
  file('*.table') into gtkrcl_multiqc
  tuple val(type), val(sampleID), file('*.bqsr.bam'), file('*.bqsr.bam.bai') into (gatkhc_in, gmultimetric_in)
  file("${sampleID}.GATK4_BQSR.log.txt") into bqsr_log

  script:
  intlist = "${intlisd}/${params.vepGenome}_${params.vepVersion}.interval_list"
  """
  {
  INBAM=$bam
  if [[ $type == "RNA" ]];then
    mkdir tmp
    gatk SplitNCigarReads \
      -R $fa \
      -I $bam \
      --tmp-dir tmp \
      -O scn.bam
    samtools index scn.bam
    INBAM=scn.bam
    rm -rf tmp
  fi

  INTTEST=\$(grep -m1 "@SQ" $intlist | perl -ane 'print \$F[0];')
  if [[ ! \$INTTEST == "@SQ" ]]; then
    perl -ane 'print "\$F[0]\\n";' $intlist > int.list
    INTLIST=int.list
  else
    cp $intlist interval.list.interval_list
    INTLIST=interval.list.interval_list
  fi

  gatk BaseRecalibrator \
    -R $fa \
    -I \$INBAM \
    --known-sites $dbsnp \
    --use-original-qualities \
    -O ${sampleID}.recal_data.table \
    --disable-sequence-dictionary-validation true \
    -L \$INTLIST

  #ApplyBQSR
  OUTBAM=\$(echo $bam | sed 's/bam/bqsr.bam/')
  gatk ApplyBQSR \
    -R $fa \
    -I \$INBAM \
    --bqsr-recal-file ${sampleID}.recal_data.table \
    --add-output-sam-program-record \
    --use-original-qualities \
    -O \$OUTBAM \
    -L \$INTLIST

  samtools index \$OUTBAM \$OUTBAM".bai"
  } 2>&1 | tee >  ${sampleID}.GATK4_BQSR.log.txt
  """
}

/* 1.4: GATK4 Germline
*/
process gatkHC {

  label 'medmem'

  publishDir "${params.outDir}/DNA_and_RNA/samples/$sampleID/gatk4", mode: "copy", pattern: "*.log.txt"
  publishDir path: "${params.outDir}/DNA_and_RNA/output/vcf", mode: "copy", pattern: '*.vcf.*'

  input:
  tuple val(type), val(sampleID), file(bam), file(bai) from gatkhc_in
  file(fa) from ref.fa
  file(vcfd) from ref.vcf
  file(intlisd) from ref.intlist

  output:
  tuple val(type), val(sampleID), file("${sampleID}.${type}.vcf.gz"), file("${sampleID}.${type}.vcf.gz.tbi") into merge_in
  val(sampleID) into vcfGRaID
  file("${sampleID}.GATK4_HC.log.txt") into log_gatkgerm

  script:
  taskmem = javaTaskmem("${task.memory}")
  intlist = "${intlisd}/${params.vepGenome}_${params.vepVersion}.interval_list"
  vcf = "${vcfd}/${params.vepGenome}_${params.vepVersion}.vcf.gz"
  """
  {
  INTTEST=\$(grep -m1 "@SQ" $intlist | perl -ane 'print \$F[0];')
  if [[ ! \$INTTEST == "@SQ" ]]; then
    perl -ane 'print "\$F[0]\\n";' $intlist > int.list
    INTLIST=int.list
  else
    cp $intlist interval.list.interval_list
    INTLIST=interval.list.interval_list
  fi

  gatk --java-options -Xmx$taskmem HaplotypeCaller \
    -R $fa \
    -I $bam \
    -ERC NONE \
    --dont-use-soft-clipped-bases \
    --standard-min-confidence-threshold-for-calling 20 \
    --dbsnp $vcf \
    --pair-hmm-implementation FASTEST_AVAILABLE \
    --native-pair-hmm-threads ${task.cpus} \
    -O $sampleID"."$type".vcf" \
    --disable-sequence-dictionary-validation true \
    -L \$INTLIST

  bgzip $sampleID"."$type".vcf"
  tabix $sampleID"."$type".vcf.gz"

  } 2>&1 | tee > ${sampleID}.GATK4_HC.log.txt
  """
}

/* 1.25: CPSR annotation of GATK4 Germline
*/
//collect and map
merge_in
  .collect()
  .map { it -> tuple(it[0], it[1..-1].flatten()) }
  .set { merge_inc }

process mergevcfs {

  label 'highmem'

  publishDir "${params.outDir}/DNA_and_RNA/output/vcf", mode: "copy"

  input:
  tuple val(first), file(vcfs) from merge_inc

  output:
  file('all.dna_and_rna.merge.vcf') into vep_vcf

  script:
  taskmem = javaTaskmem("${task.memory}")
  """
  {
    VCFS=\$(ls *.vcf.gz)
    bcftools merge -O z \$VCFS > all.dna_and_rna.merge.vcf

  } 2>&1 | tee > mergevcfs.log.txt
  """
}

/* 2.1: PicardTools metrics suite for MultiQC HTML report
*/
process mltmet {

  label 'lowmem'

  publishDir "${params.outDir}/DNA_and_RNA/samples/$sampleID/metrics"

  input:
  tuple val(type), val(sampleID), file(bam), file(bai) from gmultimetric_in
  file(fa) from ref.fa
  file(intlisd) from ref.intlist

  output:
  file('*.txt') into multimetrics_multiqc

  when:
  params.multiMetrics

  script:
  taskmem = javaTaskmem("${task.memory}")
  intlist = "${intlisd}/${params.vepGenome}_${params.vepVersion}.interval_list"
  """
  {
  picard -Xmx$taskmem CollectHsMetrics \
    I=$bam \
    O=$sampleID".hs_metrics.txt" \
    TMP_DIR=./ \
    R=$fa \
    BAIT_INTERVALS=$intlist  \
    TARGET_INTERVALS=$intlist

  picard -Xmx$taskmem CollectAlignmentSummaryMetrics \
    I=$bam \
    O=$sampleID".AlignmentSummaryMetrics.txt" \
    TMP_DIR=./ \
    R=$fa

  picard -Xmx$taskmem CollectMultipleMetrics \
    I=$bam \
    O=$sampleID".CollectMultipleMetrics.txt" \
    TMP_DIR=./ \
    R=$fa

  picard -Xmx$taskmem CollectSequencingArtifactMetrics \
    I=$bam \
    O=$sampleID".artifact_metrics.txt" \
    TMP_DIR=./ \
    R=$fa

  picard -Xmx$taskmem CollectInsertSizeMetrics \
    I=$bam \
    O=$sampleID".insert_size_metrics.txt" \
    H=$bam".histogram.pdf" \
    TMP_DIR=./

  } 2>&1 | tee > ${sampleID}.picard.metrics.log
  """
}

/* 3.0 Run multiQC to finalise report
*/
process mltiQC {

  publishDir path: "${params.outDir}/DNA_and_RNA/reports", mode: "copy", pattern: "*html"

  input:
  file(fastps) from fastp_multiqc.collect()
  file(mrkdups) from mrkdup_multiqc.collect()
  file(starqcs) from starLOG_multiqc.collect()
  file(gtkrcls) from gtkrcl_multiqc.collect()
  file(multimetrics) from multimetrics_multiqc.collect()

  output:
  file('*') into completedmultiqc

  script:
  """
  OUTID=\$(basename ${workflow.launchDir})".dna_and_rna"
  multiqc . -i \$OUTID -f -c ${params.multiqcConfig}
  """
}

process vepann {

  publishDir path: "${params.outDir}/DNA_and_RNA/output/vcf", mode: "copy", pattern: '*.vcf'

  input:
  file(vcf) from vep_vcf
  file(fa) from ref.fa
  file(vepcache) from ref.vep

  output:
  file('*vep.vcf') into vcf_stats

  script:
  vcfanno = "${vcf}".replaceAll(".vcf", ".vep.vcf")
  """
  vep --dir_cache ${vepcache} \
    --offline \
    --assembly ${params.vepGenome} \
    --vcf_info_field ANN \
    --symbol \
    --species ${params.vepSpecies} \
    --check_existing \
    --cache \
    --fork ${task.cpus} \
    --vcf \
    --input_file ${vcf} \
    --output_file ${vcfanno} \
    --format "vcf" \
    --fasta ${fa} \
    --hgvs \
    --canonical \
    --ccds \
    --force_overwrite \
    --verbose
  """
}

// process vcf_report {
//
//   publishDir path: "${params.outDir}/DNA_and_RNA/output/vcf", mode: "copy", pattern: '*.vcf'
//
//   input:
//   file(vcf) from vcf_stats
//   file(intlist) from Channel.value([params.intlist])
//
//   output:
//   file('*.vcf') into vcf_stats
//
//   script:
//   """
//   TOTMB=\$(grep - "@" $intlist | awk '{s+=\$3-\$2}END{print s}')
//   perl ${workflow.projectDir}/bin/parse_merged_vcf.pl $vcf
//   TAB=\$(ls $vcf | sed 's/vcf\$/tab/')
//   Rscript --vanilla ${workflow.projectDir}/bin/summary_vcf_tab.R \$TOTMB ${workflow.projectDir}/bin/summary_vcf_tab.Rnw
//   """
// }
