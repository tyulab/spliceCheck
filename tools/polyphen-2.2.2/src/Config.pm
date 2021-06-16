package PPH::Config;
use strict;
use warnings;

#DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD

=head1 NAME

PPH::Config - global variable container

=head1 DESCRIPTION

- config files from the subfolder config will be read
- contains also some fundamental values

=head1 AUTHOR

StS

=head1 SUBVERSION

 $LastChangedDate: 2012-03-04 19:09:14 -0500 (Sun, 04 Mar 2012) $
 $LastChangedRevision: 395 $
 $LastChangedBy: ivan $

=cut

#DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD

#BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB

our %CONFIG;
our @CONFKEYS;
our $DEBUG;
our @DEBUGGING_MODE;

BEGIN {
    use Exporter ();
    our ($VERSION, @ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);
    $VERSION = 0.1;
    @ISA           = qw(Exporter);
    @EXPORT        = qw($DEBUG @DEBUGGING_MODE %CONFIG @CONFKEYS);
    @EXPORT_OK     = qw();
    %EXPORT_TAGS = ();
}

our @EXPORT_OK;

#BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB

use Pod::Usage;
use File::Copy;
use File::stat;
use File::Temp;
use File::Basename;
use Cwd qw(abs_path);
use CGI;

# System architecture
our $ARCH = `uname -m`; chomp $ARCH;

#------------------------------------------------------------------------------
# Reading in the config files - stripping off blanks, comments, etc.
# Priority is: .pph:~/.pph:$ENV{PPH}/config
# Note that
#------------------------------------------------------------------------------
my $CNFPATH;
if      (-d '.pph') {
  $CNFPATH = '.pph';
} elsif (exists $ENV{'HOME'} && -d "$ENV{HOME}/.pph") {
  $CNFPATH = "$ENV{HOME}/.pph";
} elsif (exists $ENV{'PPH'} && -d "$ENV{PPH}/config") {
  $CNFPATH = "$ENV{PPH}/config";
} else {
  die "ERROR: Unable to locate configuration directory\n";
}
$CONFIG{'CNFPATH'} = $CNFPATH;
read_config() or die "ERROR: No valid configuration files found in: $CNFPATH\n";

#------------------------------------------------------------------------------

#==============================================================================
# TEMPORARY FILE EXTENSIONS
#------------------------------------------------------------------------------
# Files created in TMP_OUT need to be unique
#==============================================================================
$CONFIG{'PIDTIME'} = $$ . '_' . time();

#==============================================================================
# DEBUG flag
#==============================================================================
# Enables debugging, requires Smart::Comments installed
# DEBUGGING_MODE controls verbosity, 'perldoc Smart::Comments' for details
#----------------------------------------
$DEBUG = exists $CONFIG{DEBUG} ? $CONFIG{DEBUG} : 0;
if ($DEBUG) {
  # Enable Smart::Comments and set verbosity level
  for (1..$DEBUG) { push @DEBUGGING_MODE, '##'.('#' x $_); }
  $ENV{Smart_Comments} = join(':', @DEBUGGING_MODE);
  # Output of helper programs redirected to a file in the current working directory
  $CONFIG{LIMBO} = 'DEBUG.log';
} else {
  # Disable Smart::Comments
  $ENV{Smart_Comments} = 0;
  # Output of helper programs is normally ignored
  $CONFIG{LIMBO} = '/dev/null';
}
#==============================================================================

# Save the original %CONFIG before variable substitution
my %CONFIG_ORIG = %CONFIG;
eval_values();

#==============================================================================
# For more precise checking you can set the type of global variable here
#==============================================================================

#------------------------------------------------------------------------------
# Variables containing executables
#------------------------------------------------------------------------------
my %EXE;
@EXE{ qw(COILS FASTACMD BLAST DSSPCMBI HBPLUS BFACT DIST BLAT MSA
         CLUSPACK PSIC FETCHSEQ TWOBIT) } = ();

#------------------------------------------------------------------------------
# Variables containing location of files starting with a prefix or acutal file
#------------------------------------------------------------------------------
my %FILE;
@FILE{ qw(NRDB_BLAST NRDB_FASTA UNIPROT LIMBO WEKAMODEL) } = ();

#------------------------------------------------------------------------------
# Variables containing location of a dir
#------------------------------------------------------------------------------
my %DIR;
@DIR{ qw(PPH BIN MSA LEON WEKA PDB DSSP PDB2FASTA GOLDENPATH MATRIX
         MULTIZPATH PRECOMPATH SCRATCH TMP_OUT) } = ();

#==============================================================================

#------------------------------------------------------------------------------

=head2 read_config

  Function: Simple routine to read config files
  Params:   List of config filename tags (optional)
  Returns : List of full filenames read

=cut

#----------------------------------------
sub read_config {
  my @conftags = @_;
  my $confpath = $CONFIG{'CNFPATH'};
  die "Missing CNFPATH in config\n" unless length $confpath;
  # Compile a list of tags for a directory if it was not passed to us
  unless (@conftags) {
    # Compile a list of available config files, either .cnf or .cnf.dist
    my @conflist = glob "$confpath/*.{cnf,cnf.dist}";
    foreach my $cfnfile (@conflist) {
      my ($fname, $fpath, $fext) = fileparse($cfnfile, qw{ .cnf .cnf.dist });
      push @conftags, $fname;
    }
    # Ensure that the 'paths' config is always first on the list and 'options' if present is always the second one
    @conftags = sort { return -1 if $a eq 'options'; return 1 if $b eq 'options'; $a cmp $b } unique(@conftags);
    @conftags = sort { return -1 if $a eq 'paths'; return 1 if $b eq 'paths' } @conftags;
  }
  my @filesread;
  foreach my $f (@conftags) {
    # Read .cnf file if present, otherwise try reading .cnf.dist
    my $cnf = "$confpath/$f";
    $cnf .= '.cnf' unless $cnf =~ /\.cnf$/;
    unless (-e $cnf && -f $cnf && -s $cnf) {
      my $dist = $cnf . '.dist';
      if (-e $dist && -f $dist && -s $dist) { $cnf = $dist;}
      else { next; }
    }
    # Reading config file
    open my $F, '<', $cnf or die "Can't read config file: $cnf\n";
    while (<$F>) {
        chomp;
        s/^\s+//; s/\s+$//;
        next if /^#/;
        next unless length;
        my ($var, $value) = split(/\s*=\s*/, $_, 2);
        push @CONFKEYS, $var unless exists $CONFIG{$var};
        $CONFIG{$var} = $value;
    }
    close $F;
    push @filesread, $cnf;
  }
  return @filesread;
}
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------

=head2 eval_values

  Usage   : eval_values()
  Function: Routine to parse preferences recursively in order
            to allow variables in their definitions
  Returns : TRUE

=cut

#----------------------------------------
sub eval_values {
  my %quoted;
  # Recursion level is limited to 3!
  for (1..3) {
    foreach my $k (@CONFKEYS) {
      my $v = $CONFIG{$k};
      unless (defined $v) {
        die "CONFIGURATION FILE ERROR!\n" .
          "   Sorry, there seems to be an error in the configuration files!\n" .
          "   Please check the settings for '$k'\n";
      }
      # Do not attempt any interpolation if a value
      # string is quoted (as either "..." or '...'),
      # instead remove quotes surrounding the string
      # value and record the variable name to protect
      # it from further interpolation attempts.
      if ($v =~ /^(["'])/ && $v =~ /($1)$/) {
        my $qchar = $1;
        $v =~ s/^$qchar//; $v =~ s/$qchar$//;
        $CONFIG{$k} = $v;
        $quoted{$k} = 1;
      }
      # Enable warnings to catch possible references to undefined vars;
      eval qq{ use warnings; \$CONFIG{$k} = qq{$CONFIG{$k}} } unless exists $quoted{$k};
    }
  }
  return 1;
}
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------

=head2 test_path

 Usage   : test_patch($p)
 Function: Simply check if the path argument makes sense
           (i.e. hits directory or file)
 Returns : TRUE if ok / otherwise FALSE

=cut

#----------------------------------------
sub test_path {
    my $v = shift;
    my $f = shift;

    if (exists($EXE{$v})) {
        my $st = stat($f);
        return (-f $f) ? 1 : 0;
        return ($st->mode & 0111 ) ? 1 : 0;
    } elsif (exists($DIR{$v})) {
        # Check for 5 required subdirs
        if ($v eq 'SCRATCH' or $v eq 'PRECOMPATH') {
          return (
            -d "$f/alignments" and
            -d "$f/blastfiles" and
            -d "$f/profiles"   and
            -d "$f/structures" and
            -d "$f/lock"
          ) ? 1 : 0;
        }
        return (-d $f) ? 1 : 0;
        return 1;
    } elsif (exists($FILE{$v})) {
        my @f = glob("$f*");
        return scalar(@f);
    } else {
        warn "Unknown config option: $v\n";
        my @ok = glob("$f*");
        return scalar(@ok);
    }
}
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------

=head2 test_preferences()

  Usage   : PPH::Config->test_preferences();
  Function: Simple routine to test if the file setttings are o.k.
            It simply tests all variables starting with "/\S+" if
            glob("$var*") is ok
  Returns : TRUE if all ok

=cut

#----------------------------------------
sub test_preferences {
    ### Testing Preferences
    use Term::ANSIColor;

    foreach my $k (keys %CONFIG) {
        my $v = $CONFIG{$k};
        next unless $v =~ /^\/\S+/ and $k ne 'GRID_PATH';

        my @ok =  glob("$v*");

        unless ( test_path($k, $v) ) {
            print color("red");
            printf "? %25.25s = %s\n",  $k, $v;
            print color("reset")
        } else {
            print color("green"), "ok", color("reset");
            printf "%25.25s = %s\n",  $k, $v;
        }
    }
    return 1;
}
#------------------------------------------------------------------------------


#------------------------------------------------------------------------------

=head2 set_preferences

  Usage   : PPH::Config->>set_preferences($filename);
  Function: A routine to set variables and test them if the look ok
  Returns : TRUE

=cut

#----------------------------------------
sub set_preferences {
    my (undef, $file) = @_;
    $| = 1;
    $file =~ s/\s+$//; $file =~ s/^\s+//;

    #----------------------------------------
    # Reading the config file and put output
    # into the variable out
    #----------------------------------------
    open my $R, '<', $file or die "Can't open file: $file\n";
    my $out;
    while (<$R>) {
        #----------------------------------------
        # Blank or comment line
        #----------------------------------------
        (/^#/ or /^\s*$/) and do { $out .= $_; print; next };

        #----------------------------------------
        # Stripping the values
        #----------------------------------------
        chomp; s/#.*//; s/^\s+//; s/\s+$//;
        next unless length;
        my ($var, $value) = split(/\s*=\s*/, $_, 2);
        $value =~ s/^\s+//; $value =~ s/\s+$//;
        my $evalue = $value;
        for (1..3) { $evalue = eval qq{ qq{$evalue} }; }
        printf "%-16s= %s\n", $var, $evalue;

        #----------------------------------------
        # Now check paths and ask if ok
        #----------------------------------------
        while (1) {
            if ($evalue =~ /^\/\S+/ and $var ne 'GRID_PATH') {
                if (test_path($var, $evalue)) {
                    print color("green"),
                        "Path seems to be ok", color("reset"), "\n";
                } else {
                    if ($var eq 'SCRATCH' or $var eq 'PRECOMPATH') {
                      print color("red"),
                        "The path misses required subdirectories listed above!", color("reset"), "\n";
                    }
                    print color("red"),
                        "Is the path correct? ", color("reset"), "\n";
                }
            }
            print "Keep it (Y/n)? ";
            my $answer = <STDIN>;
            if ($answer =~ /^\s*[Yy]\s*$/ or $answer =~ /^\s*$/) {
                $value = $evalue if $var eq 'PPH';  # always expand PPH
                $out .= sprintf "%-16s= %s\n", $var, $value;
                $CONFIG{$var} = $value;
                last;
            } else {
                $evalue = $value = _change_value($var, $value);
            }
        } # end while (1)
    } # end while (<R>)
    close $R;

    #----------------------------------------
    # Alright - everything ok so far
    # Now create a backup file and overwrite
    # the config file.
    # IAA: if we are creating new config file from a default template
    # then we save it as a .cnf file and we do not make backups.
    #----------------------------------------
    my $save;
    if ($file =~ /\.cnf\.dist$/) {
      $file =~ s/\.dist$//;
    } else {
      my ($year, $month, $day) = (localtime(time()))[5,4,3];
      $year += 1900; $month++;
      $save = $file . ".$year.$month.$day.$$.bak";
      copy($file, $save) or die "Can't copy $file to $save\n";
    }

    open my $W, '>', $file or die "Can't write to $file";
    print $W "$out";
    close $W;
    print color('red'), "FINISHED\n", color('reset');
    if (defined $save) {
      print "Backup file stored: $save\n\n";
    } else {
      print "New config file created: $file\n\n";
    }
}

#------------------------------------------------------------------------------

#------------------------------------------------------------------------------

=head2 _change_value($$)

 Usage   : _change_value($var, $value)
 Function: Ask for the new value and return it if o.k.
 Input   : var   - simply the input text to display
           value - the standard value

=cut

#------------------------------------------------------------------------------
sub _change_value {
    my ($var, $value) = @_;
    printf "%-16s= ",  $var;
    $value = <STDIN>;
    $value =~ s/^\s+//; $value =~ s/\s+$//;
    return $value;
}
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
sub unique {
  my %seen;
  return grep { ! $seen{$_}++ } @_;
}
#------------------------------------------------------------------------------

sub new {
  my ($class, $opts) = @_;

  # Restore original %CONFIG before variable substitution
  %CONFIG = %CONFIG_ORIG;

  # Read configuration from user-specified directory
  # Note: overrides previously read configuration!
  if (my $confpath = $opts->{'c'}) {
    (-d $confpath) or pod2usage("Directory does not exist: $confpath");
    %CONFIG = ();
    $CONFIG{'CNFPATH'} = $confpath;
    read_config() or die "ERROR: No valid configuration files found in: $confpath\n";
  }

  if (my $queryfile = $opts->{'q'}) {
    my $abspath = abs_path($queryfile);
    die "ERROR: Input file does not exist: $queryfile\n"
      unless defined $abspath && length $abspath && -e $abspath;
    $queryfile = $abspath;
    # Input file is a CGI.pm query in Boulder format
    open(INP, $queryfile) or die "ERROR: Can't open query input file: $queryfile\n";
    my $query = new CGI(\*INP);
    close(INP);
    # --- Mandatory query params ---
    # accid seqres seqpos seqvar1 seqvar2
    # --- Advanced (structural) query params, optional ---
    # SORTBYIDE MAP3D2MISMATCH STRUCTALLHITS CONTALLHITS
    # MIN3DHITLEN MIN3DHITIDE MAX3DHITGAPS CONTTHRESH
    my @qpars_mandatory = qw( accid seqres seqpos seqvar1 seqvar2 );
    my @qpars_optional  = qw{ SORTBYIDE MAP3D2MISMATCH STRUCTALLHITS CONTALLHITS MIN3DHITLEN MIN3DHITIDE MAX3DHITGAPS CONTTHRESH };
    # Clean up parameters that have been input by user
    foreach my $name (qw{ accid seqres seqpos }) {
      my $value = $query->param($name);
      if (defined $value && length $value) {
        my $value = $query->param($name);
        $value =~ s/^\s+//s; $value =~ s/\s+$//s;
      } else {
        $value = '';
      }
      $query->param($name, $value);
    }
    # Validate accid & seqres params (mutually exclusive)
    my ($accid, $seqres, $seqfile);
    $accid  = $query->param('accid');
    $seqres = $query->param('seqres');
    if (length $accid && length $seqres) {
      die "Either protein identifier or dbSNP rsID, or a protein sequence in FASTA format should be specified, but not both\n";
    } elsif (!(length $accid || length $seqres)) {
      die "Neither protein identifier nor protein sequence in FASTA format were specified\n";
    # Ensure that we have only a single protein ID submitted
    } elsif (length $accid) {
      die "Invalid protein identifier specified: $accid\n"
        if scalar(my @a = split(' ', $accid)) > 1;
    }
    # Validate rest of mandatory parameters
    my $dbsnp = 0;
    foreach my $name (@qpars_mandatory) {
      # dbSNP rsID instead of a protein ID is supported since v2.0.20
      $dbsnp = 1, last if $name eq 'accid' && $query->param($name) =~ /^rs\d+$/;
      next if $name eq 'accid' || $name eq 'seqres';
      die "Mandatory query parameter missing: $name\n" unless defined $query->param($name);
      die "Mandatory query parameter value missing for: $name\n" unless length $query->param($name);
    }
    # Parse FASTA sequence for ID, then save sequence to a temp file
    if (length $seqres) {
      my ($defline) = $seqres =~ /^>(\S.+)(\n|$)/;
      die "ERROR: Submitted sequence have a missing or malformed FASTA definition line\n" unless defined $defline && length $defline;
      require PPH::Align; # we can't do 'use PPH::Align' since Align.pm does 'use PPH::Config'
      $accid = (PPH::Align::_parse_defline($defline))[2];
      $seqfile = File::Temp->new( UNLINK=>0 );  # FIXME
      print $seqfile $seqres;
      $seqfile->flush;  # flush buffers (Perl 5.8.0 and above)
      $CONFIG{'SEQFILE'} = $seqfile;
    }
    # Set the rest of parameters only for param names defined in CONFIG
    foreach my $name (@qpars_optional) {
      my $value = $query->param($name);
      $CONFIG{$name} = $value if defined $value && exists $CONFIG{$name};
    }
    # Construct query string
    if ($dbsnp) {
      $CONFIG{'QUERYSTRING'} = sprintf "%s\n", $accid;
    } else {
      $CONFIG{'QUERYSTRING'} = sprintf "%s\t%s\t%s\t%s\n",
        $accid, $query->param('seqpos'), $query->param('seqvar1'), $query->param('seqvar2');
    }
  }

  # Now process rest of command-line options (will override options set via -q query_file)
  #----------------------------------------
  # Debug mode
  #----------------------------------------
  if (defined $opts->{'v'}) {
    pod2usage("Unsupported verbosity level: $opts->{v}") if $opts->{'v'} =~ /[^\d]/;
    $CONFIG{'DEBUG'} = $DEBUG = $opts->{'v'};
    if ($DEBUG) {
      # Enable Smart::Comments and set verbosity level
      for (1..$DEBUG) { push @DEBUGGING_MODE, '##'.('#' x $_); }
      $ENV{Smart_Comments} = join(':', @DEBUGGING_MODE);
      # Output of helper programs redirected to a file in the current working directory
      $CONFIG{LIMBO} = 'DEBUG.log';
    } else {
      # Disable Smart::Comments
      $ENV{Smart_Comments} = 0;
      # Output of helper programs is normally ignored
      $CONFIG{LIMBO} = '/dev/null';
    }
  }
  #----------------------------------------

  #-----------------------------------------------
  # Set scratch directory (create it if necessary)
  #-----------------------------------------------
  if (defined $opts->{'d'}) {
    unless (-d $opts->{'d'}) {
      unless (mkdir $opts->{'d'}) {
        my $msg = "$!";
        # Ignore "file exists" errors due to possible race conditions
        die "ERROR: Can't create scratch directory $opts->{d}: $msg\n"
          unless $msg =~ /File exists/i;
      }
    }
    foreach my $subdir (qw{ alignments blastfiles profiles structures lock }) {
      unless (-d "$opts->{d}/$subdir") {
        unless (mkdir "$opts->{d}/$subdir") {
          my $msg = "$!";
          # Ignore "file exists" errors due to possible race conditions
          die "ERROR: Can't create directory $opts->{d}/$subdir: $msg\n"
            unless $msg =~ /File exists/i;
        }
      }
    }
    $CONFIG{'SCRATCH'} = abs_path($opts->{'d'});
  }
  #----------------------------------------

  #----------------------------------------
  # Precomputed only?
  #----------------------------------------
  if ($opts->{'p'}) {
    $CONFIG{'PRFPRECOMPONLY'} = 1;
    $CONFIG{'PRECOMPONLY3D'}  = 1;
  }
  #----------------------------------------

  #----------------------------------------
  # Old polyPhen-1 MSA pipeline?
  #----------------------------------------
  if ($opts->{'t'}) {
      $CONFIG{'USEBLASTMSA'} = 1;
  }
  #----------------------------------------

  #----------------------------------------
  # Use full length sequences for MSA?
  #----------------------------------------
  if ($opts->{'f'}) {
      $CONFIG{'USEFULLSEQ'} = 1;
  }
  #----------------------------------------

  #----------------------------------------
  # Set genome assembly version
  #----------------------------------------
  if ($opts->{'g'}) {
    my $geneset = $opts->{'g'};
    if (grep { $_ eq $geneset; } qw{ hg19 hg18 }) {
      $CONFIG{'GENESET'} = $geneset;
    } else {
      pod2usage("Unsupported genome assembly version specified: $geneset");
    }
  }
  #----------------------------------------

  #----------------------------------------
  # Map protein sequences to genes?
  #----------------------------------------
  if ($opts->{'m'}) {
      $CONFIG{'MAPGENE'} = 1;
  }
  #----------------------------------------

  #----------------------------------------
  # Set BLAST E-value cutoff
  #----------------------------------------
  if (exists $opts->{'e'} && length $opts->{'e'}) {
      $CONFIG{'BLAST_EVALUE'} = $opts->{'e'};
  }
  #----------------------------------------

  #----------------------------------------
  # Output saved in XML file?
  #----------------------------------------
  if ($opts->{'x'}) {
      $CONFIG{'XMLOUTPUT'} = 1;
  }
  #----------------------------------------

  #----------------------------------------
  # Enable/disable job array mode?
  #----------------------------------------
  if ($opts->{'a'}) {
      $CONFIG{'GRID_ARRAY'} = 1;
  }
  if ($opts->{'A'}) {
      $CONFIG{'GRID_ARRAY'} = 0;
  }
  #----------------------------------------

  #----------------------------------------
  # Set BLAST nr database pathname
  #----------------------------------------
  if ($opts->{'b'}) {
      my ($name, $path) = fileparse($opts->{'b'});
      (-d $path) or pod2usage("Directory does not exist: $path");
      $path = abs_path($path);
      $CONFIG{'NRDB_BLAST'} = "$path/$name";
      $CONFIG{'NRDB_FASTA'} = $CONFIG{'NRDB_BLAST'};
  }
  #----------------------------------------

  #----------------------------------------
  # Sequence file?
  # Note: if -q option is specified then -s will be ignored
  #       since the sequence will be read from a query file
  #----------------------------------------
  if (defined $opts->{'s'} && !$opts->{'q'}) {
    pod2usage("Sequence file does not exist: $opts->{s}") unless -e $opts->{'s'};
    $CONFIG{'SEQFILE'} = abs_path($opts->{'s'});
  }
  #----------------------------------------

  eval_values();
}

1;
