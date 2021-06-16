#!/usr/local/bin/perl

# This program is licensed to you under the Fred
# Hutchinos Cancer Research Center (FHCRC)
# NONCOMMERICAL LICENSE.  A copy of the license may be found at
# http://blocks.fhcrc.org/sift/license.html and should be attached
# to this software.
$file = @ARGV[0];
$queryseq = @ARGV[1];
$database = @ARGV[2];

print "hello $file query $queryseq data $database\n";

open (QUERYSEQ,">$queryseq");
open (DATABASE,">$database");

open (FILE,$file);
$count = 0;
$reached_first_seq = 0; # added to get rid of html crap 

while (<FILE>){
	if ($flag == 1){
		if ($prev ne ""){
			print DATABASE "$prev";
		}
		print DATABASE "$_";	
		$prev = "";
	}
	else{
		if ($_ =~ /\>/){
			$count ++;
			$reached_first_seq = 1;
			if ($count == 2){
				$prev = $_;
				$flag = 1;
				next;
			}
			else{
				$prev = "";
			}
		}	
		if ($reached_first_seq) {
			print QUERYSEQ "$_";
		}
	}
	
}

close (FILE);
close (QUERYSEQ);
close (DATABASE);
exit (0);
	
