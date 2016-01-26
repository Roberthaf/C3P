###### Program description
#
# Authors:
#	 Jelena Calyseva
#	 Gad Hatem
#	 Robert Anton Hafthorsson
#	
# Description:
#	 This program improves gene prediction for complement C3 protein of 3 different gene prediction methods - 
#	 Augustus, FGENESH and GENSCAN. Taking prediction built by these prediction methods, this program 
#	 identifies the best predicted exons and combines them into new predicted gene, giving exon characteristics,
#	 coding DNA sequence and encoded protein sequence as output.
#	
# Modules:
#
#	fgenes, genscan, augustus - gather the information from
#		prediction software output
#
#	blast - aligns exon to all the exons of reference
#		and returns alignment score, e-value and the length of alignment
#		of the highest scoring alignment (references to arrays)
#
# Procedure:
#	1. Extract the information from Augustus, GENSCAN and FGENESH prediction outputs and store all of it in arrays.
#	2. Run quality check for every predicted exon:
#			1) Use BLAST algorithm to align it to all known exons and identify the highest scoring one;
#			2) Use conservation score, e-value and the length of alignment to identify the quality of
#																			the predicted exon;
#			3) Check the length of the exon comparing it to the known exons;
#			4) Store it as the best result for the exon of this number if the characteristics are the best.
#	3. Build a final set of exons in which every exon is the best exon of this number from the input.
#	4. Check if the phases of adjacent exons match and make adjustments if they don't.
#	5. Output the characteristics of final exons, DNA sequence and protein sequence in a "prediction.txt" file.
#
# usage: program.pl (input files separated by spaces)

use strict;
use Data::Dumper;
use genscan;
use aug;
use fgenes;
use blast;

# this hash holds the information on average lengths of exons (from the known ones)
my %length_ref = ( '1' => 77, '2' => 191, '3' => 88, '4' => 71, '5' => 95,
'6' => 83, '7' => 91,'8' => 103,'9' => 127,'10' => 116,'11' => 148,
'12' => 210,'13' => 207,'14' => 158,'15' => 130,'16' => 72,'17' => 199,
'18' => 111,'19' => 86,'20' => 143,'21' => 213,'22' => 67,'23' => 87,'24' => 204,
'25' => 76,'26' => 160,'27' => 99,'28' => 157,'29' => 164,'30' => 159,
'31' => 60,'32' => 91,'33' => 52,'34' => 88,'35' => 88,'36' => 106,
'37' => 90,'38' => 84,'39' => 84,'40' => 136,'41' => 142);
# this hash holds the information on the lowest alignment score for every exon (aligning known exons from different species to each other)
my %check_high_scores = ('1' => 80, '2' => 74, '3' => 83, '4' => 86, '5' => 83, '6', 88, '7' => 75, '8'=> 87, 
'9'=> 83, '10'=> 97, '11'=> 78, '12'=> 81, '13'=> 88, '14'=> 79, '15'=> 86, '16' => 90, '17'=> 77,
'18'=> 86, '19'=> 100, '20'=> 93, '21'=> 91, '22'=> 90, '23'=> 92, '24'=> 94, '25'=> 83, '26'=> 90, '27'=> 84,
'28'=> 75, '29'=> 90, '30'=> 83, '31'=> 90, '32'=> 86, '33'=> 87, '34'=> 85, '35'=> 90, '36'=> 94, '37'=> 82, 
'38'=> 81, '39'=> 70, '40'=> 84, '41'=> 91);
# all these arrays are created to store the predicted exons of the best quality outside of t"The Big Loop"
# they all are not empty, because we need 42 positions (position corresponds to the number of exon) which are 0 or less until defined
my @good_length = ( -41 .. 0 );
my @good_start = ( -41 .. 0 );
my @good_end = ( -41 .. 0 );
my @good_sequence_dna = ( -41 .. 0 );
my @good_sequence_prot = ( -41 .. 0 );
my @good_start_phase = ( -41 .. 0 );
my @good_end_phase = ( -41 .. 0 );
my @final_start_phase = ( -41 .. 0 );
my @final_end_phase = ( -41 .. 0 );
my @final_length = ( -41 .. 0 );
my @final_start = ( -41 .. 0 );
my @final_end = ( -41 .. 0 );
my @final_sequence_dna = ( -41 .. 0 );	
my @final_sequence_prot = ( -41 .. 0 );
my @bad_length = ( -41 .. 0 );
my @bad_start = ( -41 .. 0 );
my @bad_end = ( -41 .. 0 );
my @bad_sequence_dna = ( -41 .. 0 );
my @bad_sequence_prot = ( -41 .. 0 );
my @bad_start_phase = ( -41 .. 0 );
my @bad_end_phase = ( -41 .. 0 );
my @worst_length = ( -41 .. 0 );
my @worst_end = ( -41 .. 0 );
my @worst_start = ( -41 .. 0 );
my @worst_sequence_dna = ( -41 .. 0 );
my @worst_sequence_prot = ( -41 .. 0 );
my @worst_start_phase = ( -41 .. 0 );
my @worst_end_phase = ( -41 .. 0 );
# these arrays are the quality check arrays, for every exon (position in the array) only the best score will be stored
my @best = ( -41 .. 0 );
my @best_difference = ( -41 .. 0 );
my @worst = ( -41 .. 0 );
my @worst_difference = ( -41 .. 0 );
my @bad = ( -41 .. 0 );
my @bad_difference = ( -41 .. 0 );
# the hash of error messages, exons numbers being the keys, values being "1" if the warning message should be printed
# in the output, empty if shouldn't
my %error;
# the codon table for checking the codons which are split between 2 exons
my %aacode = (
  TTT => "F", TTC => "F", TTA => "L", TTG => "L",
  TCT => "S", TCC => "S", TCA => "S", TCG => "S",
  TAT => "Y", TAC => "Y", TAA => "STOP", TAG => "STOP",
  TGT => "C", TGC => "C", TGA => "STOP", TGG => "W",
  CTT => "L", CTC => "L", CTA => "L", CTG => "L",
  CCT => "P", CCC => "P", CCA => "P", CCG => "P",
  CAT => "H", CAC => "H", CAA => "Q", CAG => "Q",
  CGT => "R", CGC => "R", CGA => "R", CGG => "R",
  ATT => "I", ATC => "I", ATA => "I", ATG => "M",
  ACT => "T", ACC => "T", ACA => "T", ACG => "T",
  AAT => "N", AAC => "N", AAA => "K", AAG => "K",
  AGT => "S", AGC => "S", AGA => "R", AGG => "R",
  GTT => "V", GTC => "V", GTA => "V", GTG => "V",
  GCT => "A", GCC => "A", GCA => "A", GCG => "A",
  GAT => "D", GAC => "D", GAA => "E", GAG => "E",
  GGT => "G", GGC => "G", GGA => "G", GGG => "G",
);

# the references to be extracted from the prediction methods outputs
my ($gsstart_ref, $gsend_ref, $gslength_ref, $gsseqdna_ref, $gsseqprot_ref, $gsstartph_ref, $gsendph_ref);
my ($agstart_ref, $agend_ref, $aglength_ref, $agseqdna_ref, $agseqprot_ref, $agstartph_ref, $agendph_ref);
my ($fgstart_ref, $fgend_ref, $fglength_ref, $fgseqdna_ref, $fgseqprot_ref, $fgstartph_ref, $fgendph_ref);

# values to check if the input parsing has already been done or not
my $g = 0;
my $a = 0;
my $f = 0;

# collecting the information from the input
for ( my $i = 0; $i < 3; $i++ ) {
	open my $FH, '<', $ARGV[$i];
	while ( my $read = <$FH> ) {
		if ( ($read =~ /GENSCAN/) and ($g == 0) ){
		($gsstart_ref, $gsend_ref, $gslength_ref, $gsseqdna_ref, $gsseqprot_ref, $gsstartph_ref, $gsendph_ref) = genscan ($FH);
		$g = 1;
		}
		elsif ( ($read =~ /AUGUSTUS/) and ($a == 0)  ){ 
		($agstart_ref, $agend_ref, $aglength_ref, $agseqdna_ref, $agseqprot_ref, $agstartph_ref, $agendph_ref) = aug ($FH);
		$a = 1;
		}
		elsif ( ($read =~ /FGENESH/) and ($f == 0) ){
		($fgstart_ref, $fgend_ref, $fglength_ref, $fgseqdna_ref, $fgseqprot_ref, $fgstartph_ref, $fgendph_ref) = fgenes ($FH);
		$f = 1;
		}
	}
	close $FH;
} 

# print error message if the file is not found
if ($g == 0) { print "NO GENSCAN FILE\n" }
elsif ($a == 0) { print "NO AUGUSTUS FILE\n" }
elsif ($f == 0) { print "NO FGENESH FILE\n"}

# combine the input
my @start = (@$agstart_ref, @$gsstart_ref, @$fgstart_ref);
my @end = (@$agend_ref, @$gsend_ref, @$fgend_ref);
my @length = (@$aglength_ref, @$gslength_ref, @$fglength_ref);
my @sequence_dna = (@$agseqdna_ref, @$gsseqdna_ref, @$fgseqdna_ref);
my @sequence_prot = (@$agseqprot_ref, @$gsseqprot_ref, @$fgseqprot_ref);
my @start_phase = (@$agstartph_ref, @$gsstartph_ref, @$fgstartph_ref);
my @end_phase = (@$agendph_ref, @$gsendph_ref, @$fgendph_ref);

# the output file
my $file = 'prediction.txt';                
open my $OUT, '>', $file;

# The Big Loop starts
# in this loop, every predicted exon is checked for quality and stored if there are 
# no predicted exons of the same number with better score
for ( my $i = 0; $i < @length; $i++) {

	my $sequence_dna = $sequence_dna[$i];
	my $sequence_prot = $sequence_prot[$i];
	my $length_prot = length $sequence_prot[$i];
	my $length = $length[$i];
	my $start = $start[$i];
	my $end = $end[$i];
	
	# align to the known exons
	my $exon_nr = 0;
	my $score = 0;
	my $evalue = 1;
	my $length_hsp = 0;
	my $qseq = '';
	my $difference = 0;
	if ($length_prot > 0) {
		($exon_nr, $score, $evalue, $length_hsp, $qseq) = blast($sequence_prot);
		$difference = $length_hsp / $length_prot ;	# checking if the predicted exon is not too long comparing to the one it aligns to
	}
	if ( $length > 250 ) {
		my $position = index $sequence_prot, $qseq;
		my $count = ($position*3) + $start_phase[$i];
		if ( $count == 0 ) { $start = $start + $count } else { $start = $start + $count - 1 }
		$length = (length $qseq)*3;
		my $new_end = $end;
		$end = $start + $length - 1;
		my $new_start = $end + 1;
		my $new_length = $new_end - $new_start + 1;
		my $new_dna = substr $sequence_dna, (($position*3) + $start_phase[$i] + $length);
		$sequence_dna = substr $sequence_dna, (($position*3)+$start_phase[$i]), $length;
		my $new_prot = substr $sequence_prot, ($position + (length $qseq));
		$sequence_prot = $qseq;
		my $new_end_ph = $end_phase[$i];
		$start_phase[$i] = 0;
		$end_phase[$i] = 0;
		push @start, $new_start;
		push @end, $new_end;
		push @length, $new_length;
		push @sequence_dna, $new_dna;
		push @sequence_prot, $new_prot;
		push @start_phase, 0;
		push @end_phase, $new_end_ph;	
	}

	# check the alignment scores
	
	my $good_score = 0;
	my $bad_score = 0;
	
	# the score is defined as good if the alignment score is better than the one of reference for exon of this number
	if ( ($score >= $check_high_scores{$exon_nr}) and ( $evalue < 0.001 ) ){
		$good_score = $score;
	}	
	elsif ( $evalue < 0.001 ) {
		$bad_score = $score;
	}
	print $OUT "$i: $exon_nr exon, $evalue e-value, $bad_score bad score, $good_score good score\n";
	# the score is ignored if e-value is 0.001 or higher
	
	# check the length of the exon 														
	if ( $good_score > 0 ) {
		my $check = $length_ref{$exon_nr};		# from the hash of average lengths
		my $min = $check - 6;		 
		my $max = $check + 6;
		# if the length is in a reasonable range, we consider this exon of good quality and put
		# it in "good" or "first choice" arrays
		if ( ($length > $min) and ($length < $max) ) {
			# the exon is put in one of the final arrays only if the position is not filled with a better exon yet
			if ( ($good_score >= $best[$exon_nr]) and ($difference > $best_difference[$exon_nr])) {
				$good_length[$exon_nr] = $length;
				$good_start[$exon_nr] = $start;
				$good_end[$exon_nr] = $end;
				$good_sequence_dna[$exon_nr] = $sequence_dna;
				$good_sequence_prot[$exon_nr] = $sequence_prot;
				$good_start_phase[$exon_nr] = $start_phase[$i];
				$good_end_phase[$exon_nr] = $end_phase[$i];
				$best[$exon_nr] = $good_score;				# the score is considered the best score for this exon
				$best_difference[$exon_nr] = $difference;
			}
		}
		# if the score is good but the length is not in the range, the exon is put into a "bad", or "second choice",
		# set of arrays
		elsif ( ($good_score >= $bad[$exon_nr]) and ($difference > $bad_difference[$exon_nr]) ) {
				$bad_length[$exon_nr] = $length;
				$bad_start[$exon_nr] = $start;
				$bad_end[$exon_nr] = $end;
				$bad_sequence_dna[$exon_nr] = $sequence_dna;
				$bad_sequence_prot[$exon_nr] = $sequence_prot;
				$bad_start_phase[$exon_nr] = $start_phase[$i];
				$bad_end_phase[$exon_nr] = $end_phase[$i];
				$bad[$exon_nr] = $good_score;				# the score is considered the best score for the second choice exons
				$bad_difference[$exon_nr] = $difference;
		}
	}
	# if the score is off, but e-value is reasonable, the exon is put in the "worst", or
	# "last choice" set of arrays
	elsif ( ($bad_score >= $worst[$exon_nr]) and ($difference > $worst_difference[$exon_nr]) ) {
			$worst_length[$exon_nr] = $length;
			$worst_start[$exon_nr] = $start;
			$worst_end[$exon_nr] = $end;
			$worst_sequence_dna[$exon_nr] = $sequence_dna;
			$worst_sequence_prot[$exon_nr] = $sequence_prot;
			$worst_start_phase[$exon_nr] = $start_phase[$i];
			$worst_end_phase[$exon_nr] = $end_phase[$i];
			$worst[$exon_nr] = $bad_score;
			$worst_difference[$exon_nr] = $difference;
	}	
	# printing the percentage of the progress
	my $per = int(($i/@start)*100);
	print "\b\b\b\b$per%";
}
print "\b\b\b\bComputation is complete";


# create final set of exons ( prediction itself )
# in every array, position is the number of an exon
for ( my $o = 1; $o < 42; $o++) {
	# if the "first choice" array is not empty in that position, we take that exon as the final exon
	if ( $good_end[$o] > 0 ) {
		$final_end[$o] =  $good_end[$o];
		$final_length[$o] = $good_length[$o];
		$final_start[$o] = $good_start[$o];
		$final_sequence_dna[$o] = $good_sequence_dna[$o];
		$final_sequence_prot[$o] = $good_sequence_prot[$o];
		$final_start_phase[$o] = $good_start_phase[$o];
		$final_end_phase[$o] = $good_end_phase[$o];
	}
	# if the "first choice" array is empty, we take the "second choice"
	elsif ( $bad_end[$o] > 0 ) {
		$final_end[$o] =  $bad_end[$o];
		$final_length[$o] = $bad_length[$o];
		$final_start[$o] = $bad_start[$o];
		$final_sequence_dna[$o] = $bad_sequence_dna[$o];
		$final_sequence_prot[$o] = $bad_sequence_prot[$o];
		$final_start_phase[$o] = $bad_start_phase[$o];
		$final_end_phase[$o] = $bad_end_phase[$o];
	}
	# if the first and second choice arrays are empty in the position, we take the last choice
	elsif ( $worst_end[$o] > 0 ) {
		$final_end[$o] =  $worst_end[$o];
		$final_length[$o] = $worst_length[$o];
		$final_start[$o] = $worst_start[$o];
		$final_sequence_dna[$o] = $worst_sequence_dna[$o];
		$final_sequence_prot[$o] = $worst_sequence_prot[$o];
		$final_start_phase[$o] = $worst_start_phase[$o];
		$final_end_phase[$o] = $worst_end_phase[$o];
		$error{$o} = 1;			# notifying the program to print a warning message ( as the exon is of poor quality )
	}
	# if all 3 (good, bad and worst) sets of arrays are empty in the position, we leave the exon undefined
	else {
		$final_end[$o] = "NA";
		$final_length[$o] = "0";
		$final_start[$o] = "NA";
		$final_sequence_dna[$o] = "";
		$final_sequence_prot[$o] = "";
		$final_start_phase[$o] = "0";
		$final_end_phase[$o] = "0";
	}
}

# check for the phases and correct translation of the "middle" codons ( codons which start in the end
# of one exon and end in the beginning of the next one )
# loop through all codons, starting with the second, comparing phases of the start with the end of the previous one
for (my $v = 2; $v < @final_end; $v++) {
	if ( $final_start_phase[$v] == 0 ) {
		unless ( $final_end_phase[$v-1] == 0 ) { # in this case, we leave it alone, the phases match
			# if the end phase of a previous exon is not 0, the exon sequence is reduced
			my $phase = $final_end_phase[$v-1];
			$final_sequence_dna[$v-1] = substr $final_sequence_dna[$v-1], 0, -($phase);
			$final_end[$v-1] = $final_end[$v-1] - $phase;
			$final_length[$v-1] = $final_length[$v-1] - $phase;
			# if the end phase of a previous exon is 2, the amino acid sequence of the exon is also reduced
			if ( $final_end_phase[$v-1] == 2 ) {
				$final_sequence_prot[$v-1] = substr $final_sequence_prot[$v-1], 0, -1;
			}
			# the phase is then 0
			$final_end_phase[$v-1] = 0;
		}
	}
	elsif ( $final_start_phase[$v] == 1 ) {
		# if the phases don't match, we make adjustments
		if ( $final_end_phase[$v-1] != 2 ) {
			# changes in the later exon
			$final_sequence_dna[$v] = substr $final_sequence_dna[$v], 1;
			$final_start[$v] = $final_start[$v] + 1; 
			$final_length[$v] = $final_length[$v] - 1;
			$final_start_phase[$v] = 0;
			}
		if ( $final_end_phase[$v-1] == 1 ) {
			# changes in the previous exon, if the phase is 1
			$final_sequence_dna[$v-1] = substr $final_sequence_dna[$v-1], 0, -1;
			$final_end[$v-1] = $final_end[$v] - 1;
			$final_end_phase[$v-1] = 0;
			$final_length[$v-1] = $final_length[$v-1] - 1;
			}
		# if the exons match, we have to check if the codon encodes the right amino acid
		if ( $final_end_phase[$v-1] == 2 ) {
			# extract the codon
			my $start = substr $final_sequence_dna[$v-1], -2;
			my $end = substr $final_sequence_dna[$v], 0, 1;
			my $codon = $start . $end;
			# extract the amino acid
			my $seq = $final_sequence_prot[$v-1];
			my $aa = chop $seq;
			# check if they match, if not, make adjustments
			unless ($aacode{$codon} eq $aa) {
				$final_sequence_prot[$v-1] = $seq;
				$final_sequence_prot[$v-1] . $aacode{$codon};
			}
		}
	}
	elsif ( $final_start_phase[$v] == 2 ) {
		if ( $final_end_phase[$v-1] == 0 ) {
			$final_sequence_dna[$v] = substr $final_sequence_dna[$v], 2;
			$final_start[$v] = $final_start[$v] + 2;
			$final_start_phase[$v] = 0;
			$final_length[$v] = $final_length[$v] - 2;
		}
		elsif ( $final_end_phase[$v-1] == 1 ) {
			my $start = substr $final_sequence_dna[$v-1], -1;
			my $end = substr $final_sequence_dna[$v], 0, 2;
			my $codon = $start . $end;
			my $aa = substr $final_sequence_prot[$v], 0, 1;
			unless ($aacode{$codon} eq $aa) {
				my $seq = substr $final_sequence_prot[$v], 1;
				$final_sequence_prot[$v] = $aacode{$codon} . $seq;
			}
		}
		# in this case, the first nucleotide of the later exon is deleted
		elsif ( $final_end_phase[$v-1] == 2 ) {
			my $start = substr $final_sequence_dna[$v-1], -2;
			my $end = substr $final_sequence_dna[$v], 1, 2;
			my $codon = $start . $end;
			my $aa = substr $final_sequence_prot[$v-1], -1, 1;
			unless ($aacode{$codon} eq $aa) {
				my $seq = substr $final_sequence_prot[$v-1], -1;
				$final_sequence_prot[$v-1] = $seq . $aacode{$codon}
			}				
			$final_sequence_dna[$v] = substr $final_sequence_dna[$v], 1;
			$final_start[$v] = $final_start[$v] + 1;
			$final_start_phase[$v] = 1;
			$final_sequence_prot[$v] = substr $final_sequence_prot[$v], 1;
			$final_length[$v] = $final_length[$v] - 1;
		}
	}
}

# printing the output
print $OUT "Nr\tstart\tend\tlength\ts\te\tcomment\n";
my $dna_sequence;
my $protein_sequence;
for ( my $o = 1; $o < @final_sequence_prot; $o++ ) { # for every final exon
	print $OUT "$o\t$final_start[$o]\t$final_end[$o]\t$final_length[$o]\t\t$final_start_phase[$o]\t$final_end_phase[$o]\t";
	if ( $error{$o} ) { print $OUT "!"}		# printing the warning message
	print $OUT "\n";
	$dna_sequence .= $final_sequence_dna[$o];	# concatenating the sequence to print out the final one
	$protein_sequence .= $final_sequence_prot[$o];
}

print $OUT "DNA sequence:\n";
until ( (length $dna_sequence) == 0 ) {
	if ( (length $dna_sequence) < 60 ) { # the last line
		print $OUT "$dna_sequence\n";
		$dna_sequence = '';
	}
	else {
		my $line = substr $dna_sequence, 0, 60;	# new line every 60 letters
		$dna_sequence = substr $dna_sequence, 60;
		print $OUT "$line\n";
	}
}
print $OUT "Protein sequence:\n";
until ( (length $protein_sequence) == 0 ) {
	if ( (length $protein_sequence) < 60 ) {
		print $OUT "$protein_sequence\n";
		$protein_sequence = '';
	}
	else {
		my $line = substr $protein_sequence, 0, 60;
		$protein_sequence = substr $protein_sequence, 60;
		print $OUT "$line\n";
	}
}
close $OUT;