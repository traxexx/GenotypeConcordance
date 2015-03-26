#!/usr/bin/perl

# GT Concordance
# run on Indel only

use strict;

my $InVcf = $ARGV[0];
my $PairList = $ARGV[1];
my $Output = $ARGV[2];

# do pair list first
my @Pairs1;
my @Pairs2;
open( LIST, "<", $PairList ) or die $!;
while( <LIST> ) {
	chomp; my $line = $_;
	my @s = split "\t", $line;
	if ( scalar(@s) != 2 ) {
		die "List fields is not 2 at "."$."."\n";
	}
	push(@Pairs1, $s[0]);
	push(@Pairs2, $s[1]);
}
close LIST;

# open
my @vcf_name = split /\./, $InVcf;
if ( $vcf_name[-1] eq "gz" ) {
        open( IN, "zcat $InVcf|" ) or die $!;
}
else {
        open IN, "<", $InVcf or die $!;
}

# calculate
my $annotation_skipped = 0;
my $prev_line;
my @Index1;
my @Index2;
my $SampleStart;
my $GTindex;
my $PairCount = scalar(@Pairs1);

my @F1;
my @F2;
my @F3;
my @F4;
#my %hash = (0,0,1,0,2,0,3,0,4,0,5,0,6,0,7,0,8,0,9,0);

for( my $i = 0; $i < $PairCount; $i++ ) {
	push(@F1, {0=>0,1=>0,2=>0,3=>0,4=>0,5=>0,6=>0,7=>0,8=>0} );
	push(@F2, {0=>0,1=>0,2=>0,3=>0,4=>0,5=>0,6=>0,7=>0,8=>0} );
	push(@F3, {0=>0,1=>0,2=>0,3=>0,4=>0,5=>0,6=>0,7=>0,8=>0} );
	push(@F4, {0=>0,1=>0,2=>0,3=>0,4=>0,5=>0,6=>0,7=>0,8=>0} );
}


while( <IN> ) {
	chomp; my $line = $_;
	if ( $annotation_skipped == 0 ) {
		my @per_char_split = split "", $line, 2;
        if ( $per_char_split[0] eq '#' ) {
        	$prev_line = $line;
            next;
        }
        else { # prev line contains sample info
          # get GL index
            my @prev_fields = split "\t", $prev_line, 10;
            my @line_fields = split "\t", $line, 10;
            for (my $i=0; $i<10; $i++) {
            	if ( $prev_fields[$i] eq "FORMAT" ) {
            		$SampleStart = $i + 1;
            		my @info = split ":", $line_fields[$i];
            		for( my $ii=0; $ii<scalar(@info); $ii++ ) {
            			if ( $info[$ii] eq "GT" ) {
            				$GTindex = $ii;
            				last;
            			}
        			}
    			}
    		}
    		my @fields = split "\t", $prev_line; #split prev line here!
    	# get pair index
    		for( my $fc = 0; $fc < $PairCount; $fc++ ) {
    			for( my $i = $SampleStart; $i < scalar(@fields); $i++ ) {
    				if ( $fields[$i] eq $Pairs1[ $fc ] ) {
    					push( @Index1, ($i - $SampleStart) );
    				}
    				elsif( $fields[$i] eq $Pairs2[ $fc ] ) {
    					push( @Index2, ($i - $SampleStart) );
    				}
    			}
    		}
            # set parameters
            $annotation_skipped = 1;
        }
    }
# here get info of this line
    my @fields = split "\t", $line;
    my @gts;

    for( my $i = $SampleStart; $i < scalar(@fields); $i++ ) {
    	my @subfields = split ":", $fields[$i];
    # get gt
    	my @gt_str = split "", $subfields[$GTindex];
    	if ( $gt_str[0] > 1 ) {
    		$gt_str[0] = 1;}
    	if ( $gt_str[2] > 1 ) {
    		$gt_str[2] = 1;}
    	my $gt_number = $gt_str[0] + $gt_str[2];
    	push( @gts, $gt_number );
    }
# gt pr(g)
	my %Pg = ( 0,0,1,0,2,0 );
	my $maf;
	for( my $i = 0; $i < scalar(@gts); $i++ ) {
		my $gt_tmp = $gts[ $i] ;
		$Pg{ $gt_tmp } = $Pg{ $gt_tmp } + 1;
	}
# adjust to fraction
	my $sum = $Pg{'0'} + $Pg{'1'} + $Pg{'2'};
	$maf = 2 * $Pg{'2'} + $Pg{'1'};
	$maf = $maf / $sum / 2;
# get pair info
	for( my $pr = 0; $pr < scalar( @Index1 ); $pr++ ) {
		my $index1 = $Index1[$pr];
		my $index2 = $Index2[$pr];
		my $gt1 = $gts[ $index1 ];
		my $gt2 = $gts[ $index2 ];
		my $concur = $gt1 * 3 + $gt2;
		if ( $maf < 0.001) {
			$F1[$pr]{ $concur } = $F1[$pr]{ $concur } + 1;
		}
		elsif ( $maf < 0.01 ) {
			$F2[$pr]{ $concur } = $F2[$pr]{ $concur } + 1;
		}
		elsif ( $maf < 0.1 ) {
			$F3[$pr]{ $concur } = $F3[$pr]{ $concur } + 1;
		}
		else { # >= 10%
			$F4[$pr]{ $concur } = $F4[$pr]{ $concur } + 1;
		}
	}
}
close IN;

# print final
my $fname = "$Output.0";
open( OUT, ">", $fname ) or die $!;
for( my $i= 0; $i < $PairCount; $i++ ) {
	print OUT "$Pairs1[$i]-$Pairs2[$i]";
	for( my $j = 0; $j < 9; $j++ ) {
		print OUT "\t$F1[$i]{$j}";
	}
	print OUT "\n";
}	

my $fname = "$Output.1";
open( OUT, ">", $fname ) or die $!;
for( my $i= 0; $i < $PairCount; $i++ ) {
	print OUT "$Pairs1[$i]-$Pairs2[$i]";
	for( my $j = 0; $j < 9; $j++ ) {
		print OUT "\t$F2[$i]{$j}";
	}
	print OUT "\n";
}

my $fname = "$Output.2";
open( OUT, ">", $fname ) or die $!;
for( my $i= 0; $i < $PairCount; $i++ ) {
	print OUT "$Pairs1[$i]-$Pairs2[$i]";
	for( my $j = 0; $j < 9; $j++ ) {
		print OUT "\t$F3[$i]{$j}";
	}
	print OUT "\n";
}

my $fname = "$Output.3";
open( OUT, ">", $fname ) or die $!;
for( my $i= 0; $i < $PairCount; $i++ ) {
	print OUT "$Pairs1[$i]-$Pairs2[$i]";
	for( my $j = 0; $j < 9; $j++ ) {
		print OUT "\t$F4[$i]{$j}";
	}
	print OUT "\n";
}
































    