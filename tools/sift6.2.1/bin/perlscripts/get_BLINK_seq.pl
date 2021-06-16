#!/usr/local/bin/perl 
#use strict;
use LWP::Simple;
use LWP::UserAgent;
use strict;

# This program is licensed to you under the Fred
# Hutchinos Cancer Research Center (FHCRC)
# NONCOMMERICAL LICENSE.  A copy of the license may be found at
# http://blocks.fhcrc.org/sift/license.html and should be attached
# to this software

# retrieves sequence from NCBI database

my $gi_number = $ARGV[0];
my $outfile = $ARGV[1];

my $best_hits_or_all_hits_option = "";
if (@ARGV >= 3) {
	$best_hits_or_all_hits_option = $ARGV[2];
}

my $URL = "https://www.ncbi.nlm.nih.gov/sutils/blink.cgi?pid=" . $gi_number;


#Dec 2016 - changed to 399 because original gi# will be added to 400
if ($best_hits_or_all_hits_option eq "BEST") {
#	$URL .= "&cut=100&org=1";
	$URL .= "&org=1&per_page=399";
##25/6/2013 HueyTyng modified - to recognise another option CG_ONLY
##CG_ONLY will take sequences from complete genomes db only + hide identical hits
} elsif ($best_hits_or_all_hits_option eq "CG_ONLY" ) {
        $URL .= "&org=3&set=4&per_page=399";
} else {
	$URL .= "&per_page=399";
}
## December 24,2009 edits to make timeout longer
my $browser = LWP::UserAgent->new();
$browser->timeout (60);
my $request = HTTP::Request->new (GET => $URL);
my $response = $browser->request ($request);

#my $NCBI_content = get ($URL);
#if (!defined ($NCBI_content) || $NCBI_content eq "") {
#	print "ERROR!  Unable to contact NCBI right now for $gi_number .  Please try later.\n";
#	exit (-1);
#}

if ($response->is_error()) { 
	print "ERROR! %s ERROR267\n", $response->status_line;
}

my $NCBI_content = $response->content();

if ($NCBI_content =~ /Could not retrieve/) {
	print "ERROR!  NCBI could not retrieve the required information. Please make sure your sequence is a PROTEIN sequence and then try <A HREF=\"http://sift-dna.org/www/SIFT_seq_submit2.html\">submitting the sequence</A> to SIFT.\n"; 
	exit (-1);
}

if ($NCBI_content =~ /No hits found/ ) {
	# ERROR
	print "ERROR!  No BLAST hits were found for gi:$gi_number by NCBI BLink .  Please check that this sequence is a protein.<BR>You can also submit your protein sequence directly to SIFT and have SIFT choose the protein sequences  <A HREF=\"http://sift-dna.org/www/SIFT_seq_submit2.html\">here.</A>\n";

	exit (-1);
}
if ($NCBI_content =~ /ERROR/) {
	print "ERROR trying to retrieve sequences from NCBI!  Please follow this <A HREF=$URL>link</A HREF> to NCBI to see the error. \n"; 
	exit (-1);
}

my @lines = split ("\n", $NCBI_content);
if (@lines == 0) {
	print "ERROR! Unable to retrieve sequences from NCBI. Please try later.\n";
	exit (-1); 

}
 
my $lineindex = 0;
my %seen;
my $gi_hits = $gi_number . ",";
while ($lineindex < @lines ) {
	my $loc = index ($lines[$lineindex], "name='gi");
#	print "location is " . $loc;
	while ($loc >= 0) {
		my $gi = substr ($lines[$lineindex], $loc + 8);
		$loc = index ($gi, "'");
		$gi = substr ($gi, 0, $loc);
		if (!exists ($seen{$gi})) {
			$gi_hits .= $gi . "," ;
			$seen{$gi} = 1;
		}
		$lines[$lineindex] = substr ($lines[$lineindex], $loc+1, 100000);
		$loc = index ($lines[$lineindex], "name='gi");
	}
	$lineindex++;
}
$gi_hits =~ s/\,$//;
print "GI hits" . $gi_hits;
#old 
#my $seq_URL = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&email=sift\@fhcrc.org&id=$acc_hits&rettype=fasta&retmode=text"; 
#new
my $seq_URL = "https://www.ncbi.nlm.nih.gov/sviewer/viewer.fcgi?db=protein&qty=1&c_start=1&list_uids=$gi_hits&uids=&dopt=fasta&dispmax=400&sendto=t&from=begin&to=end";
my $seq = get ($seq_URL);


# Pauline Sept 10,2010
# format seq because the NCBI name is gi|###|#### and psiblast will
# format it to gi|####
my @seq_lines = split (/\n/, $seq);
for (my $i= 0; $i < @seq_lines; $i++) {
	if ($seq_lines[$i] =~ /^\>gi\|/) {
		$_ = $seq_lines[$i];
		/^\>gi\|([0-9]+)\|(.*)/;
		$seq_lines[$i] = ">gi" . $1 . " " . $2;
#print "replacing with $seq_lines[$i]\n";
	}
} 
$seq = join ("\n", @seq_lines);

open (OUT, ">$outfile") || die "can't open $outfile";
print OUT $seq ; 
close (OUT);

exit (0);

 
