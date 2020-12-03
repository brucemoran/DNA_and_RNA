#! /usr/bin/env nextflow

//script to download and parse references for DNA_and_RNA
//autoruns if --refDir is not supplied, or doesn't exist, or not all required files therein
//can be used alone to generate DNA, RNA references based on user supplying ENSEMBL URLs in config

params.help = ""

if (params.help) {
  log.info ''
  log.info '----------------------------------------------------------------------------'
  log.info 'NEXTFLOW: MAKE REFERENCE FILES AND INDICES FOR DNA_and_RNA NEXTFLOW PIPELINE'
  log.info '----------------------------------------------------------------------------'
  log.info ''
  log.info 'Usage: '
  log.info 'nextflow run refs.nf -profile standard,singularity'
  log.info ''
  log.info 'Optional arguments:'
  log.info '--exomeBed  STRING  path/to/exome.bed; exome corresponding to DNAseq data; determines regions included in analysis'
  log.info '--outDir  STRING  path/to/output/; where output is stored (default: ./refs/vepgenome)'
  log.info '--fa      STRING  URL to download fasta.gz'
  log.info '--gtf      STRING  URL to download gtf.gz'
  log.info '--cdna      STRING  URL to download cdna.gz'
  log.info '--vcf      STRING  URL to download known SNPs (vcf.gz)'
  log.info '--vepgenome     STRING  short identifier of genome assemby e.g. "Rnor_6.0"'
  log.info '--vepspecies      STRING  VEP species corresponding to vepgenome e.g. rattus_norveigcus'
  log.info '--vepversion      INT  VEP version e.g. 100'
  log.info '--sjdb      INT  length of reads in bp for STAR index (default: 100)'
  log.info ''
  exit 1
}

/* 0.0: Global Variables
*/
if(params.outDir){
  params.refDir = "$params.outDir"
}
if(!params.outDir){
  params.refDir = "refs/$params.vepgenome"
}

params.tag = Channel.from("$params.vepgenome").getVal()

/* 1.0: Download required files
*/
process downloads {

  publishDir path: "$params.refDir", mode: "copy"

  input:
  val(vepgenome) from Channel.value(params.vepgenome)

  output:
  file("${vepgenome}.fa") into (fa_dict, fa_bwa, fa_star, fa_rrna, fa_gatk4)
  file("${vepgenome}.gtf") into (gtf_star, gtf_rrna, gtf_gatk4, gtf_saf, gtf_refflat)
  file("${vepgenome}.vcf") into vcf_tabix

  """
  wget -O $vepgenome".fa.gz" ${params.ftpbase}/${params.fa}
  wget -O $vepgenome".gtf.gz" ${params.ftpbase}/${params.gtf}"."${params.vepversion}".gtf.gz"
  if [[ ${params.vcf} =~ ^http ]];then
    wget -O $vepgenome".vcf.gz" ${params.vcf}
  else
    wget -O $vepgenome".vcf.gz" ${params.ftpbase}/${params.vcf}
  fi
  gunzip $vepgenome".fa.gz"
  gunzip $vepgenome".gtf.gz"
  gunzip $vepgenome".vcf.gz"
  """
}

/* 1.1: dictionary
*/
process dictionary_pr {

  publishDir path: "$params.refDir", mode: "copy"

  input:
  file(fa) from fa_dict

  output:
  file('*.dict') into (dict_bwa, dict_window, dict_exome, dict_rrna, dict_gatk4)
  file('*') into completeddict

  """
  DICTO=\$(echo $fa | sed 's/fa/dict/')
  picard CreateSequenceDictionary \
    R=$fa \
    O=\$DICTO

  perl -ane 'if(\$F[0]=~m/SQ\$/){@sc=split(/:/,\$F[1]);@ss=split(/:/,\$F[2]); if(\$sc[1]!~m/[GLMT]/){ print "\$sc[1]\\t\$ss[1]\\n";}}' \$DICTO > seq.dict.chr-size

  bedtools makewindows -g seq.dict.chr-size -w 35000000 | perl -ane 'if(\$F[1]==0){\$F[1]++;};print "\$F[0]:\$F[1]-\$F[2]\n";' > 35MB-window.bed
  """
}

/* 2.0: Fasta processing
*/
process bwa_index {

  publishDir path: "$params.refDir", mode: "copy"

  input:
  file(fa) from fa_bwa
  file(dict) from dict_bwa

  output:
  file('*') into complete_fa

  """
  ##https://gatkforums.broadinstitute.org/gatk/discussion/2798/howto-prepare-a-reference-for-use-with-bwa-and-gatk
  samtools faidx $fa
  bwa index -a bwtsw $fa
  """
}


/* 2.0: Bed for exome
*/
process exomebed_pr {

  publishDir path: "$params.refDir", mode: "copy"

  input:
  file(dict) from dict_exome
  file(exome_bed) from Channel.fromPath("${params.exomeBed}").getVal()

  output:
  file('*') into allexome
  file("exome.bed") into exome_tabix

  when:
  params.exomeBed

  """
  ##test if all chr in fasta are in exome
  ##test if all regions are greater than length zero
  if [[ "$exome_bed" =~ ".gz" ]]; then
    gunzip -c $exome_bed > exome_bed;
  else
    cat $exome_bed > exome_bed
  fi
  perl -ane 'if(\$F[1] == \$F[2]){\$F[2]++;} if(\$F[0] !~m/^chrM/){print join("\\t", @F[0..\$#F]) . "\\n";}' exome_bed | grep -v chrM | sed 's/chr//g' > tmp.bed

  ##always make interval list so we are in line with fasta
  picard BedToIntervalList I=tmp.bed O=exome.bed.interval_list SD=$dict

  ##BedToIntervalList (reason unknown) makes 1bp interval to 0bp interval, replace with original
  perl -ane 'if(\$F[0]=~m/^@/){print \$_;next;} if(\$F[1] == \$F[2]){\$f=\$F[1]; \$f--; \$F[1]=\$f; print join("\\t", @F[0..\$#F]) . "\\n";} else{print \$_;}' exome.bed.interval_list > exome.bed.interval_list1
  mv exome.bed.interval_list1 exome.bed.interval_list

  ##output BED
  grep -v "@" exome.bed.interval_list | cut -f 1,2,3,5 > exome.bed
  rm tmp.bed
  """
}

/* 2.1: Tabix those requiring tabixing
*/
process vcf_pr {

  publishDir path: "$params.refDir", mode: "copy"

  input:
  file(tbtbx) from exome_tabix

  output:
  file('*') into tabixd

  when:
  params.exomeBed

  """
  bgzip $tbtbx
  tabix $tbtbx".gz"
  """
}

process indexfeaturefile_pr {

  publishDir path: "$params.refDir", mode: "copy"

  input:
  file(tbtbx) from vcf_tabix

  output:
  file('*') into indexfeatured

  """
  #https://gatkforums.broadinstitute.org/gatk/discussion/7020/error-malformed-vcf-empty-alleles-are-not-permitted-in-vcf-records
  gawk 'BEGIN{FS="\\t"; OFS="\\t"}{if (NF>1 && \$5=="") {\$5="."; print \$0} else print \$0}' $tbtbx > fortbx.vcf
  mv fortbx.vcf $tbtbx
  bgzip $tbtbx
  gatk IndexFeatureFile -I $tbtbx".gz"
  """
}

/* 4.0: STAR geneomeGenerate
 */
process star_index {

  publishDir "$params.refDir", mode: "copy"

  input:
  file(fa) from fa_star
  file(gtf) from gtf_star

  output:
  file('*') into completedstargg

  script:
  """
  #! /bin/bash
  let SJDB=${params.sjdb}-1
  mkdir -p STAR_\$SJDB

  RAM=\$(echo ${task.memory} | sed 's/ G/000000000/')
  STAR --runMode genomeGenerate \
    --genomeDir STAR_\$SJDB \
    --genomeFastaFiles $fa \
    --sjdbGTFfile $gtf \
    --sjdbOverhang \$SJDB \
    --runThreadN ${task.cpus} \
    --limitGenomeGenerateRAM \$RAM
  """
}

/* 5.0: refFlat conversion of GTF
 */
process reflat_pr {

  publishDir "$params.refDir", mode: "copy", pattern: "refFlat*"

  input:
  file(gtf) from gtf_refflat

  output:
  file('*') into completedrefflat

  script:
  """
  #! /bin/bash
  wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/gtfToGenePred
  chmod a+x ./gtfToGenePred
  mkdir -p refFlat
  cd refFlat
  ../gtfToGenePred -includeVersion ../$gtf $gtf".refFlat"
  perl -ane 'print "\$F[0]\\t\$_";' $gtf".refFlat" > 1
  mv 1 $gtf".refFlat"
  """
}

/* 3.0: rRNA_intervalList
 */
process rrna_pr {

  publishDir "$params.refDir", mode: "copy", pattern: "rRNA*"

  input:
  file(fa) from fa_rrna
  file(gtf) from gtf_rrna
  file(dict) from dict_rrna

  output:
  file('*') into completedrRNA

  script:
  """
  mkdir -p rRNA
  cd rRNA

  ##https://www.biostars.org/p/67079/
  cat ../$dict > $fa".dict.rRNA.interval_list"

  grep 'gene_biotype "rRNA"' ../$gtf | \
    awk '\$3 == "transcript"' | \
    cut -f1,4,5,7,9 | \
    perl -lane '/transcript_id "([^"]+)"/ or die "no transcript_id on \$.";
      print join "\\t", (@F[0,1,2,3], \$1);' | \
    sort -k1V -k2n -k3n >> $fa".dict.rRNA.interval_list"
  """
}

/* 4.0 GATK4 SNV calling
*/
process gatk4snv_pr {

  publishDir path: "$params.refDir", mode: 'copy'

  input:
  file(fa) from fa_gatk4
  file(gtf) from gtf_gatk4
  file(dict) from dict_gatk4

  output:
  file('*') into gatk4all

  """
  samtools faidx $fa

  ##https://gatkforums.broadinstitute.org/gatk/discussion/2798/howto-prepare-a-reference-for-use-with-bwa-and-
  perl -ane 'if(\$F[0]=~m/SQ\$/){@sc=split(/:/,\$F[1]);@ss=split(/:/,\$F[2]); if(\$sc[1]!~m/[GLMT]/){ print "\$sc[1]\\t\$ss[1]\\n";}}' $dict > seq.dict.chr-size.dict

  ##make GTF bed for 'gene'
  GENEINTLIST=\$(echo $gtf | sed 's/\\.gtf/\\.genes\\.interval_list/')
  cat $gtf | perl -ane 'if(\$F[2] eq "gene"){print "\$F[0]\\t\$F[3]\\t\$F[4]\\n";}' | \
    perl -ane 'if((\$F[0]=~m/^chr/) || (\$F[0]=~m/^[0-9]/)){\$_=~s/chr//g;print \$_;}' > genes.bed

  ##test for 'chr' in fasta dictionary
  ##always make interval list so we are in line with fasta
  CHRT=\$(cat seq.dict.chr-size.dict | head -n1 | \
    perl -ane '@s=split(/:/,\$F[1]);if(\$s[1]=~m/^chr/){print "CHR\\n"};')

  if [[ \$CHRT == "" ]];then
    cat seq.dict.chr-size.dict | perl -ane 'print "\$F[0]\\n";' | \
    while read CHR; do
      export CHR;
      sort -V genes.bed | perl -ane 'if(\$F[0] eq \$ENV{CHR}){print \$_;}else{next;}'
    done > 1
    picard BedToIntervalList I=1 O=\$GENEINTLIST SD=$dict
    rm 1
  else
    picard BedToIntervalList I=genes.bed O=\$GENEINTLIST SD=$dict
  fi

  ##tabix
  bgzip \$GENEINTLIST
  gunzip -c \$GENEINTLIST".gz" > \$GENEINTLIST
  tabix \$GENEINTLIST".gz"
  """
}

/* 5.0 SAF for featureCounts
*/
process gtfsaf_pr {

  publishDir path: "$params.refDir", mode: 'copy'

  input:
  file(gtf) from gtf_saf

  output:
  file('*.saf') into saf_out

  """
  ##make SAF format from gtf
  SAF=\$(echo $gtf | sed 's/\\.gtf/\\.saf/')
  echo -e "GeneID\\tChr\\tStart\\tEnd\\tStrand" > \$SAF
  grep -v "#" $gtf | perl -ane 'chomp;if(\$F[2] eq "gene"){\$g=\$F[9];\$g=~s/[";]//g;print "\$g\\t\$F[0]\\t\$F[3]\\t\$F[4]\\t\$F[6]\\n";}' >> \$SAF
  """
}

process vep_install {

  publishDir path: "$params.refDir", mode: 'copy'

  output:
  file(".vep") into vep_cache

  script:
  """
  {
  mkdir -p .vep && cd .vep
  wget ftp://ftp.ensembl.org/pub/release-${params.vepversion}/variation/indexed_vep_cache/${params.vepspecies}_vep_${params.vepversion}_${params.vepgenome}.tar.gz

  cd ../
  vep_install \
    --ASSEMBLY ${params.vepgenome} \
    --SPECIES ${params.vepspecies} \
    --VERSION ${params.vepversion} \
    --PLUGINS dbNSFP \
    --AUTO cp \
    --CACHEURL .vep \
    --CACHEDIR .vep \
    --NO_HTSLIB \
    --NO_TEST \
    --NO_UPDATE

  } 2>&1 | tee > vep_install.log
  """
}
