#!/usr/bin/perl

# GL Concordance
# run on Indel only

use strict;

my $InVcf = $ARGV[0];
my $PairList = $ARGV[1];
my $ExcludeList = $ARGV[2];

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

# read exclude list
my @Excludes;
open( EXC, "<", $ExcludeList ) or die $!;
while( <EXC> ) {
	chomp; my $line = $_;
	push( @Excludes, $line );
}
close EXC;

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
my $GLindex;
my $GTindex;
my $PairCount = scalar(@Pairs1);
my %ExcludeIndex;
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
            			}
            			elsif ( $info[$ii] eq "PL" ) {
            				$GLindex = $ii;
            			}
        			}
    			}
    		}
    		my @fields = split "\t", $prev_line; #split current line here!
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
    	# get exclude index
    		for( my $fc = 0; $fc < scalar(@Excludes); $fc++ ) {
    			for( my $i = $SampleStart; $i < scalar(@fields); $i++ ) {
    				if ( $fields[$i] eq $Excludes[ $fc ] ) {
    					$ExcludeIndex{ ($i - $SampleStart) } = 1;
    				}
    			}
    		}
            # set parameters & print header
            $annotation_skipped = 1;
            print "Position\tMinorAF";
            for( my $i = 0; $i < $PairCount; $i++ ) {
            	print "\t$Pairs1[$i]-$Pairs2[$i]";
            }
            print "\n";
        }
    }
# here get info of this line
    my @fields = split "\t", $line;
    my @gts;
    my @gls;
    for( my $i = $SampleStart; $i < scalar(@fields); $i++ ) {
    	my @subfields = split ":", $fields[$i];
    # get gt
    	my @gt_str = split "", $subfields[$GTindex];
    	my $gt_number = $gt_str[0] + $gt_str[2];
    	push( @gts, $gt_number );
    # get gl
    	push( @gls, $subfields[$GLindex] );
    }
# gt pr(g)
	my %Pg = ( 0,0,1,0,2,0 );
	my $maf;
	for( my $i = 0; $i < scalar(@gts); $i++ ) {
		if ( exists( $ExcludeIndex{$i} ) ) {
			next;
		}
		my $gt_tmp = $gts[ $i] ;
		$Pg{ $gt_tmp } = $Pg{ $gt_tmp } + 1;
	}
  # skip pairs with no 1|0 & 1|1
	if ( $Pg{'1'} == 0 && $Pg{'2'} == 0) {
		next;
	}
# adjust to fraction
	my $sum = $Pg{'0'} + $Pg{'1'} + $Pg{'2'};
	$maf = 2 * $Pg{'2'} + $Pg{'1'};
	$maf = $maf / $sum / 2;
	$maf = sprintf("%.4f", $maf);
	for( my $ele = 0; $ele <=2; $ele++) {
		if ( $Pg{$ele} == 0 ) {
			$Pg{ $ele } = 10 ** (-6);
		}
		else {
			$Pg{ $ele } = $Pg{ $ele } / $sum;
		}
	}
# get pair info
	for( my $pr = 0; $pr < scalar( @Index1 ); $pr++ ) {
		my $no_info = 0;
		my $index1 = $Index1[$pr];
		my $index2 = $Index2[$pr];
	# check if .
		if ($gls[ $index1 ] eq "." || $gls[ $index2 ] eq ".") {
			$no_info = 1;
		}
		
		my @gl1;
		my @gl2;
		if (  $no_info == 0) {		
			@gl1 = split ",", $gls[ $index1 ];
			@gl2 = split ",", $gls[ $index2 ];
	# check if all zero
			my $gl_sum = 0;
			foreach my $element(@gl1) {
				$gl_sum = $gl_sum + $element;
			}
			if ($gl_sum == 0 ) {
				$no_info = 1;
			}
			else {
				$gl_sum = 0;
				foreach my $element(@gl2) {
					$gl_sum = $gl_sum + $element;
				}
				if ( $gl_sum == 0 ) {
					$no_info = 1;
				}
			}
		}
		
		my $bayes;
		if ($no_info == 0) { # calculate bayes factor
			for( my $j = 0; $j <= 2; $j++ ) {
				if ( $gl1[ $j ] > 100 ) {
					$gl1[ $j ] = 100;
				}
				if ( $gl2[ $j ] > 100 ) {
					$gl2[ $j ] = 100;
				}
				$gl1[ $j ] = 10 ** ( -$gl1[ $j ] / 10 );
				$gl2[ $j ] = 10 ** ( -$gl2[ $j ] / 10 );
			}
		# h0
			my $ph0 = 0;
			for( my $cg = 0; $cg <=2; $cg++ ) {
				my $single = $gl1[ $cg ] * $gl2[ $cg ] * $Pg{ $cg };
				$ph0 = $ph0 + $single;
			}
		# h1
			my $ph1;
			my $pre1 = 0;
			for( my $g1 = 0; $g1 <= 2; $g1++ ) {
				$pre1 = $pre1 + $gl1[ $g1 ] * $Pg{ $g1 };
			}
			my $pre2 = 0;
			for( my $g2 = 0; $g2 <= 2; $g2++ ) {
				$pre2 = $pre2 + $gl2[ $g2 ] * $Pg{ $g2 };
			}
			$ph1 = $pre1 * $pre2;
			$bayes = log( $ph0 / $ph1 ) / log(10);
			$bayes = sprintf("%.3f", $bayes);
#print "\npres: $pre1\t$pre2\n";
#print "\nphs: $ph0\t$ph1\n";
		}
		else { # unable to calculate
			$bayes = "NA";
		}
	# final, save
		if ( $pr == 0 ) {
			print "$fields[1]\t$maf";
		}
		print "\t$bayes";
	}
	print "\n";
}

close IN;
close OUT;






























    