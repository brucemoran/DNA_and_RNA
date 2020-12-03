#! /usr/bin/env nextflow

//script to download and parse references for DNA_and_RNA
//autoruns if --refDir is not supplied, or doesn't exist, or not all required files therein
//can be used alone to generate DNA, RNA references based on user supplying ENSEMBL URLs in config

def helpMessage() {
  log.info"""
  ---------------------------------------------------------------------------
            REFERENCE FILES FOR 'DNA_and_RNA' NEXTFLOW PIPELINE
  ---------------------------------------------------------------------------
  Usage:

  nextflow run brucemoran/DNA_and_RNA/download-references.nf

  Mandatory arguments:
    -profile      [str]   standard,singularity
    --fa          [str]   URL to download fasta.gz
    --gtf         [str]   URL to download gtf.gz
    --cdna        [str]   URL to download cdna.gz
    --vcf         [str]   URL to download known SNPs (vcf.gz)
    --vepGenome   [str]   short identifier of genome assemby e.g. "Rnor_6.0"
    --vepSpecies  [str]   VEP species corresponding to vepGenome e.g. "rattus_norveigcus"
    --vepVersion  [int]   VEP version e.g. 100
    --sjdb        [int]   length of reads in bp for STAR index (default: 100)

  Optional arguments:
    --exomeBed    [str]   path/to/exome.bedcorresponding to DNAseq data; determines regions included in analysis
  """.stripIndet()
}

if (params.help) exit 0, helpMessage()

// 0.0: Global Variables
//refDir output, tag
params.refDir = "${params.vepGenome}_${params.vepVersion}"
params.tag = Channel.from("${params.vepGenome}_${params.vepVersion}").getVal()

// 1.0: Download required files
process downloads {

  publishDir path: "${params.refDir}/gtf", mode: "copy", pattern: "*gtf"
  publishDir path: "${params.refDir}/fa", mode: "copy", pattern: "*fa"
  publishDir path: "${params.refDir}/vcf", mode: "copy", pattern: "*vcf"

  input:
  val(vepGenome) from Channel.value(params.vepGenome)

  output:
  file("${params.vepGenome}_${params.vepVersion}.fa") into (fa_dict, fa_bwa, fa_star, fa_gatk4)
  file("${params.vepGenome}_${params.vepVersion}.gtf") into (gtf_star, gtf_rrna, gtf_gatk4, gtf_saf, gtf_refflat)
  file("${params.vepGenome}_${params.vepVersion}.vcf") into vcf_tabix

  """
  wget -O ${params.vepGenome}_${params.vepVersion}".fa.gz" ${params.fa}
  wget -O ${params.vepGenome}_${params.vepVersion}".gtf.gz" ${params.gtf}
  wget -O ${params.vepGenome}_${params.vepVersion}".vcf.gz" ${params.vcf}

  gunzip ${params.vepGenome}_${params.vepVersion}".fa.gz"
  gunzip ${params.vepGenome}_${params.vepVersion}".gtf.gz"
  gunzip ${params.vepGenome}_${params.vepVersion}".vcf.gz"
  """
}

// 1.1: Dictionary of fa
process dictionary_pr {

  publishDir path: "${params.refDir}/fa", mode: "copy"
  publishDir path: "${params.refDir}/bwa", mode: "symlink", pattern: "*[.fa, .dict]"

  input:
  file(fa) from fa_dict

  output:
  file("${dict}") into (dict_bwa, dict_window, dict_exome, dict_rrna, dict_gatk4)

  script:
  dict = "${fa}".replaceAll(".fa", ".dict")
  """
  picard CreateSequenceDictionary \
    R=${fa} \
    O=${dict}
  """
}

// 2.0: Fasta processing
process bwa_index {

  publishDir path: "${params.refDir}/bwa", mode: "copy", pattern: "*[!.fa, !.dict]"

  input:
  file(fa) from fa_bwa
  file(dict) from dict_bwa

  output:
  file('*') into complete_bwa

  """
  ##https://gatkforums.broadinstitute.org/gatk/discussion/2798/howto-prepare-a-reference-for-use-with-bwa-and-gatk
  samtools faidx ${fa}
  bwa index -a bwtsw ${fa}
  """
}

// 2.2: Bed for exome
process exomebed_pr {

  publishDir path: "${params.refDir}/exome", mode: "copy"

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
  if [[ "${exome_bed}" =~ ".gz" ]]; then
    gunzip -c ${exome_bed} > exome_bed;
  else
    cat ${exome_bed} > exome_bed
  fi
  perl -ane 'if(\$F[1] == \$F[2]){\$F[2]++;} if(\$F[0] !~m/^chrM/){print join("\\t", @F[0..\$#F]) . "\\n";}' exome_bed | grep -v chrM | sed 's/chr//g' > tmp.bed

  ##always make interval list so we are in line with fasta
  picard BedToIntervalList I=tmp.bed O=exome.bed.interval_list SD=${dict}

  ##BedToIntervalList (reason unknown) makes 1bp interval to 0bp interval, replace with original
  perl -ane 'if(\$F[0]=~m/^@/){print \$_;next;} if(\$F[1] == \$F[2]){\$f=\$F[1]; \$f--; \$F[1]=\$f; print join("\\t", @F[0..\$#F]) . "\\n";} else{print \$_;}' exome.bed.interval_list > exome.bed.interval_list1
  mv exome.bed.interval_list1 exome.bed.interval_list

  ##output BED
  grep -v "@" exome.bed.interval_list | cut -f 1,2,3,5 > exome.bed
  rm tmp.bed
  """
}

// 2.3: Tabix those requiring tabixing
process vcf_pr {

  publishDir path: "${params.refDir}/exome", mode: "copy"

  input:
  file(tbtbx) from exome_tabix

  output:
  file('*') into tabixd

  when:
  params.exomeBed

  """
  bgzip ${tbtbx}
  tabix ${tbtbx}".gz"
  """
}

// 2.4: Tabix those requiring tabixing
process indexfeaturefile_pr {

  publishDir path: "${params.refDir}/vcf", mode: "copy"

  input:
  file(tbtbx) from vcf_tabix

  output:
  file('*') into indexfeatured

  """
  #https://gatkforums.broadinstitute.org/gatk/discussion/7020/error-malformed-vcf-empty-alleles-are-not-permitted-in-vcf-records
  gawk 'BEGIN{FS="\\t"; OFS="\\t"}{if (NF>1 && \$5=="") {\$5="."; print \$0} else print \$0}' ${tbtbx} > fortbx.vcf
  mv fortbx.vcf ${tbtbx}
  bgzip ${tbtbx}
  gatk IndexFeatureFile -I ${tbtbx}".gz"
  """
}

// 3.0: STAR geneomeGenerate
process star_index {

  publishDir "${params.refDir}/star_${sjdbd}", mode: "copy", pattern: "*[!.fa, !.gtf]"

  input:
  file(fa) from fa_star
  file(gtf) from gtf_star

  output:
  file('*') into completedstargg

  script:
  sjdbd = "${params.sjdb - 1}"
  ram = "${task.memory}".replaceAll(" G", "000000000")
  """
  STAR --runMode genomeGenerate \
    --genomeDir ./ \
    --genomeFastaFiles ${fa} \
    --sjdbGTFfile ${gtf} \
    --sjdbOverhang ${sjdbd} \
    --runThreadN ${task.cpus} \
    --limitGenomeGenerateRAM ${ram}
  """
}

// 4.0: refFlat conversion of GTF
process reflat_pr {

  publishDir "${params.refDir}/refflat", mode: "copy"

  input:
  file(gtf) from gtf_refflat

  output:
  file("${reff_out}") into completedrefflat

  script:
  reff_out = "${gtf}".replaceAll(".gtf", ".refFlat")
  """
  #! /bin/bash
  wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/gtfToGenePred
  chmod a+x ./gtfToGenePred
  gtfToGenePred -includeVersion ${gtf} 1 && rm gtfToGenePred
  perl -ane 'print "\$F[0]\\t\$_";' 1 > ${reff_out}
  """
}

// 5.0: rRNA_intervalList
process rrna_pr {

  publishDir "${params.refDir}/rrna", mode: "copy"

  input:
  file(gtf) from gtf_rrna
  file(dict) from dict_rrna

  output:
  file("${rrna_out}") into completedrRNA

  script:
  rrna_out = "${dict}".replaceAll(".dict", ".rRNA.interval_list")
  """
  ##https://www.biostars.org/p/67079/
  cat ${dict} > ${rrna_out}

  grep 'gene_biotype "rRNA"' ${gtf} | \
    awk '\$3 == "transcript"' | \
    cut -f1,4,5,7,9 | \
    perl -lane '/transcript_id "([^"]+)"/ or die "no transcript_id on \$.";
      print join "\\t", (@F[0,1,2,3], \$1);' | \
    sort -k1V -k2n -k3n >> ${rrna_out}
  """
}

// 6.0 GATK4 SNV calling

process gatk4snv_pr {

  publishDir path: "${params.refDir}/intlist", mode: 'copy'

  input:
  file(fa) from fa_gatk4
  file(gtf) from gtf_gatk4
  file(dict) from dict_gatk4

  output:
  file('*interval_lis*') into gatk4all

  script:
  geneintlist = "${gtf}".replaceAll(".gtf", ".genes.interval_list")
  """
  samtools faidx ${fa}

  ##https://gatkforums.broadinstitute.org/gatk/discussion/2798/howto-prepare-a-reference-for-use-with-bwa-and-
  perl -ane 'if(\$F[0]=~m/SQ\$/){@sc=split(/:/,\$F[1]);@ss=split(/:/,\$F[2]); if(\$sc[1]!~m/[GLMT]/){ print "\$sc[1]\\t\$ss[1]\\n";}}' ${dict} > seq.dict.chr-size.dict

  ##make GTF bed for 'gene'
  cat ${gtf} | perl -ane 'if(\$F[2] eq "gene"){print "\$F[0]\\t\$F[3]\\t\$F[4]\\n";}' | \
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
    picard BedToIntervalList I=1 O=${geneintlist} SD=${dict}
    rm 1
  else
    picard BedToIntervalList I=genes.bed O=${geneintlist} SD=${dict}
  fi

  ##tabix
  bgzip ${geneintlist}
  gunzip -c ${geneintlist}.gz > ${geneintlist}
  tabix ${geneintlist}.gz
  """
}

// 7.0 SAF for featureCounts
process gtfsaf_pr {

  publishDir path: "${params.refDir}/saf", mode: 'copy'

  input:
  file(gtf) from gtf_saf

  output:
  file("${saf}") into saf_out

  script:
  saf = "${gtf}".replaceAll(".gtf", ".saf")
  """
  ##make SAF format from gtf
  echo -e "GeneID\\tChr\\tStart\\tEnd\\tStrand" > ${saf}
  grep -v "#" ${gtf} | perl -ane 'chomp;if(\$F[2] eq "gene"){\$g=\$F[9];\$g=~s/[";]//g;print "\$g\\t\$F[0]\\t\$F[3]\\t\$F[4]\\t\$F[6]\\n";}' >> ${saf}
  """
}

// 8.0 VEP
process vep_install {

  publishDir path: "${params.refDir}", mode: 'copy'

  output:
  file(".vep") into vep_cache

  script:
  """
  {
  mkdir -p .vep && cd .vep
  wget ftp://ftp.ensembl.org/pub/release-${params.vepVersion}/variation/indexed_vep_cache/${params.vepSpecies}_vep_${params.vepVersion}_${params.vepGenome}.tar.gz

  cd ../
  vep_install \
    --ASSEMBLY ${params.vepGenome} \
    --SPECIES ${params.vepSpecies} \
    --VERSION ${params.vepVersion} \
    --PLUGINS dbNSFP \
    --AUTO cp \
    --CACHEURL .vep \
    --CACHEDIR .vep \
    --NO_HTSLIB \
    --NO_TEST \
    --NO_UPDATE

  } 2>&1 | tee > ${params.vepGenome}_${params.vepVersion}.vep_install.log
  """
}
