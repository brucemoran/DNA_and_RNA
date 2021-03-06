#! /usr/bin/env perl

use strict;
use warnings;

##script to parse out XLSX-friendly table from VEP-annotated VCF
##also return a per-sample mutation count (for TMB estimation or such)
##single input: *vep.vcf merged VCF of all samples required

##input
my $VCF=$ARGV[0];
my $OUT=$VCF;
$OUT=~s/vcf/tab/;

##iterate over VCF
open(IN, $VCF);
my %PERSAMP;
my $tab;
my $tabo;
my $sampn=0;
while(<IN>){
  chomp;

  if($_=~m/^##/){
    next;
  }

  my @sp=split(/\t/);

  ##header line
  ##record info about variant
  ##unique ID e.g. chr1:100_G>A
  ##and canonical VEP result
  if($_=~m/^#CHROM/){
    $tab="chr:pos_REF>ALT\tConsequence\tIMPACT\tSYMBOL\tGene\tBIOTYPE;exon\tHGVSc\tHGVSp";
    for(my $i=9;$i<@sp;$i++){
      $tab.="\t" . $sp[$i];
      $sampn++;
    }
    $tab.="\tpop_freq\n";
    next;
  }

  my @info=split(/\|/, $sp[7]);
  if(scalar(@info)>=12){
    $tabo=$sp[0] . ":" . $sp[1] . "_" . $sp[3] . ">" . "$sp[4]\t$info[1]\t$info[2]\t$info[3]\t$info[4]\t$info[7];$info[8]\t$info[10]\t$info[11]";
  }
  if(scalar(@info)<12){
    ##intergenic MODIFIER with no annotation!
    ##why doesn't Perl know that the spaces are just space and not undef?!
    $tabo=$sp[0] . ":" . $sp[1] . "_" . $sp[3] . ">" . "$sp[4]\t$info[1]\t$info[2]\t-\t-\t-\t-\t-";
  }
  ##iterate over samples
  ##record population freq. of variant
  ##output VAF per sample
  my $pof=0;
  #print $_;
  for(my $i=9;$i<@sp;$i++){
    if($sp[$i]=~m/^\.\/\./){
      $tabo.="\t.";
      next;
    }
    else{
      $pof++;
      my @spp=split(/\:/,$sp[$i]);
      my @spr=split(/\,/, $spp[1]);
      ##issues with multiple variants at same allele
      #print "@spp\n\t@spr\n";
      ##impose filter, 10 total, 3 ALT
      if(($spr[-1] ne ".") && ($spr[-1]>=3) && (($spr[0]>=7) || ($spr[0]==0))){
        my $aft=$spr[0]+$spr[-1];
        my $vaf;
        if($aft==0){
          $vaf="multi";
        }
        else{
          $vaf=sprintf("%.3f", $spr[-1]/$aft);
        }
        $tabo.="\t$vaf";
        next;
      }
      if(($spr[-1] eq ".") && (($spr[0]>=7) || ($spr[0]==0))){
        ##case of multiple variants at position
        $tabo.="\tmulti";
        next;
      }
      else{
        if($spr[0]==0){
          ##not enough reads, but fixed
          $tabo.="\tlow_fixed";
          next;
        }
        else{
          $tabo.="\tlow";
        }
      }
    }
  }
  $tabo.="\t$pof" . "/$sampn\n";
  $tab.=$tabo;
}
close IN;
open(OT, ">$OUT");
print OT $tab;
close OT;
