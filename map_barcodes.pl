#!/usr/bin/perl

use strict;
use warnings;
use Text::LevenshteinXS qw(distance);
use Getopt::Std;

my %options=();
getopts('B:', \%options);

my $fasta = $ARGV[0];
my $read = $ARGV[1];
my $out = $ARGV[2];

open (FASTA, "$fasta") or die("ERROR: can not read file ($fasta): $!\n");
open (MATCH, ">$out".".match") or die("ERROR: can not create $out .matched: $!\n");
open (REJECT, ">$out".".reject") or die("ERROR: can not create $out .rejected: $!\n");

my $link_A_bc = "TCTAGA";
my $link_P_bc = "AGTCAG";
my $link_A_oligo = "AGTG";
my $link_P_oligo = "AGGG";
my $end_A_oligo = "CGTC";
my $end_P_oligo = "GCAA";
my $MIN_SEQ_SIZE = 100;

my $barcode_seq;
my $oligo_seq;
my $barcode_start;
my $oligo_start;
my $oligo_end;
my $oligo_length;
my $id;
my $r1;
my $revcomp;
my $link_index;
my $library;

while (<FASTA>){
# Extract Sequence ID
  chomp;
  $id = $_;
  $id =~ s/^>//;
	$id =~ s/\/1$//;

# Extract the sequence
  $r1 = <FASTA>;
  chomp $r1;
  next if(length($r1) < $MIN_SEQ_SIZE);

# If checking a Read 1 file take the reverse complement
  if($read == 1){
    $revcomp = reverse($r1);
  	$revcomp =~ tr/ACGTNacgtn/TGCANtgcan/;
    $r1 = $revcomp;
  }
# Check for presence of linker A
  if(substr($r1, 18, 10) =~ $link_A_bc){
    $link_index = index($link_A_bc, substr($r1, 18, 10));
    $link_index += 18;
    $library = "A";
  }
# Check for presence of linker P
  elsif(substr($r1, 8, 10) =~ $link_P_bc){
    $link_index = index($link_P_bc, substr($r1, 18, 10));
    $link_index += 8;
    $library = "P";
  }
# Check for presence of linker A for Duo
  elsif(substr($r1, 40, 12) =~ $link_A_bc){
    $link_index = index($link_A_bc, substr($r1, 18, 10));
    $link_index += 40;
    $library = "A";
  }
# Check for presence of linker P for Duo
  elsif(substr($r1, 36, 12) =~ $link_P_bc){
    $link_index = index($link_P_bc, substr($r1, 18, 10));
    $link_index += 36;
    $library = "P";
  }
  elsif(substr($r1, 18, 10) !~ $link_A_bc & substr($r1, 8, 10) !~ $link_P_bc & substr($r1, 40, 12) !~ $link_A_bc & substr($r1, 36, 12) =~ $link_P_bc){
    print REJECT "$id\n"
  }

# Match libraries to barcodes
  if($library eq "A"){
    if($link_index - 20 < 0){
      $barcode_start = 0;
      $barcode_seq = substr($r1, $barcode_start, 20);
    }
    if($link_index - 20 >= 0){
      $barcode_start = $link_index - 20;
      $barcode_seq = substr($r1, $barcode_start, 20);
    }
  }
  if($library eq "P"){
    if($link_index - 10 < 0){
      $barcode_start = 0;
      $barcode_seq = substr($r1, $barcode_start, 10);
    }
    if($link_index - 10 >= 0){
      $barcode_start = $link_index - 10;
      $barcode_seq = substr($r1, $barcode_start, 10);
    }
  }

# Find start of oligo
  if($library eq "A"){
    $oligo_start = index($link_A_oligo, substr($r1, $barcode_start+34, 6));
    $oligo_start += 4;
  }
  if($library eq "P"){
    $oligo_start = index($link_P_oligo, substr($r1, $barcode_start+24, 6));
    $oligo_start += 4;
  }

# Find the end of the oligo
  if($library eq "A"){
    $oligo_end = index($end_A_oligo, substr($r1, -18, 6));
    $oligo_end += -14;
  }
  if($library eq "P"){
    $oligo_end = index($end_P_oligo, substr($r1, -17, 6));
    $oligo_end += -13;
  }

# Define the substring that is the oligo
  $oligo_length = length($r1) + $oligo_end - $oligo_start;
  $oligo_seq = substr($r1, $oligo_start, $oligo_length);
  print MATCH join("\t", $id, $barcode_seq, $oligo_seq, $oligo_length)
}
