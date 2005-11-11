#!/usr/local/bin/perl5.6.0 -w
#
# make_acefiles 
# dl1
#

#################################################################################
# Variables                                                                     #
#################################################################################

$|=1;
use lib "/wormsrv2/scripts/";   
use strict;
use Wormbase;
use IPC::Open2;
use IO::Handle;
use Getopt::Std;
use vars qw( $opt_n $opt_s $opt_e $opt_d $opt_h);

 ##############################
 # Paths etc                  #
 ##############################

my $autoace_config = "/wormsrv2/autoace_config/autoace.config";
my $basedir        = "/wormsrv2/wormbase";
my $miscdir        = "$basedir/misc";
my $tace           = "/nfs/disk100/acedb/RELEASE.DEVELOPMENT/bin.ALPHA_4/tace";
my $autodir        = "/wormsrv2/autoace";

 ##############################
 # command-line options       #
 ##############################

getopts ('ndh');

# help page
&usage(0) if ($opt_h);

# no option supplied
&usage(1) if (!$opt_n);


#################################################################################
# Main Loop                                                                     #
#################################################################################

&mknewacefiles   if ($opt_n);

close STDERR;
close STDOUT;

exit (0);


#################################################################################
# Erases old acefiles and make new ones                                         #
#################################################################################

sub mknewacefiles {
    
    my ($dbname,$dbdir,$targetdir,$exe,$object,$criteria,$criterianoasterisk);
    my ($filename,$command,$outfile,$tag,@deletes);

    local (*CONFIG);

    open (CONFIG, "$autoace_config");
    while (<CONFIG>) {
	s/^\s+//;s/\s+$//;    
			
	# next if comment line
	next if (/^\#/  || /^$/);

        # single database mode
	if ($opt_s) {
	    next unless (/$opt_s/);
	}

	if (/^P\s+(\S+)\s+(\S+)$/) {
	    $dbname    = $1;
	    $dbdir     = $2;
	    $targetdir = "$basedir"."/$dbname";
	    $exe       = "$tace $dbdir";
	    next;
	}
	next if ($dbname =~ /misc/);
	# single database mode

	# nuts and bolts of acefile generation
	# format:  database object criteria

	# clear out the @deletes array eachtime
	undef (@deletes);

	# deal with complete queries
	unless (/\[(\S+.+)\]$/) { 
	    if (/^\S+\s+(\S+)$/) {
		$object=$1; 
		$criteria=""; 
		$criterianoasterisk="";
	    }
	    elsif (/^\S+\s+(\S+)\s+(\S+.+)$/) {
		$object   = $1; 
		$criteria = $2;
		$criterianoasterisk = $2;
		chop($criterianoasterisk) if ($criteria =~ /\*$/);
	    } 
	}
	# deal with queries requiring tag deletions
	else {
	    (/\[(\S+.+)\]$/);
	    @deletes = split (/\s/,$1);
	    my $report = join ' + ', @deletes;
	
	    if (/^\S+\s+(\S+)\s+\[/) {  
		$object=$1; 
		$criteria=""; 
		$criterianoasterisk="";
	    } 
	    elsif (/^\S+\s+(\S+)\s+(\S+.+)\s+\[/) {
		$object   = $1; 
		$criteria = $2;
		$criterianoasterisk = $2;
		chop($criterianoasterisk) if ($criteria =~ /\*$/);
	    }
	}

	# write tace command for different options:
	# Simple   :- whole ?Class
	# Advanced :- ?Class with tag deletions based on @deletes
	# DNA      :- ?Sequence with Species, follow DNA

#	if ($object eq "DNA") {
#	    ($command,$filename) = &make_command_DNA($object,$dbname,$targetdir,$criteria);
#	}
	if ($criteria ne "") {
	    ($command,$filename) = &make_command_advanced($object,$dbname,$targetdir,$criteria,$criterianoasterisk);
	}
#	else {
#	    ($command,$filename) = &make_command_simple($object,$dbname,$targetdir,$criteria,$criterianoasterisk);
#	}
	print "Command is $command\n" if ($command);
	print "Criteria is *$criteria*\n" if ($criteria);

	
	# dump out from ACEDB
#	open (TACE,"| $exe");
#	print TACE $command;
#	close TACE;

	
	# delete things
	if (scalar (@deletes) > 0) {
	    print "Going to delete some tags ....\n";
#	    &strip_object($filename,@deletes);
	}

	# remove ghost objects
        # don't run on LongText files as these may look like ghosts

#	unless ($filename =~ /LongText/) {
#	    &strip_ghosts($filename,$object);
#	}

	# process database dumps to add database names into time stamps
	# rather than keeping the existing 'wormpub', 'lstein' etc.  name
#	&process_ace_file($filename,$dbname);

	next;
    }
    die "It's the end of the world as we know it\n";

    close CONFIG;
}

sub make_command_simple {

   my ($object,$dbname,$targetdir,$criteria,$criterianoasterisk) = @_;

   my $filename="$targetdir/".$dbname."_".$object.".ace";
   my $command;
   $command=<<END;
query find $object
show -a -T -f $filename
quit
END

   return ($command,$filename);
}

sub make_command_advanced {

   my ($object,$dbname,$targetdir,$criteria,$criterianoasterisk) = @_;

   my $filename="$targetdir/".$dbname."_".$object.".ace";
   my $command;
   $command=<<END;
query find $object where ($criteria)
show -a -T -f $filename
quit
END
   return ($command,$filename);
}

sub make_command_DNA {

    my ($object,$dbname,$targetdir,$criteria) = @_;    
    my $filename = "$targetdir/" . $dbname . "_DNA.ace";
    my $command;

    if(($dbname eq "stlace") || ($dbname eq "briggsae")){
      $command=<<END;
query find Sequence where ($criteria)
follow DNA
show -a -T -f $filename
quit
END
}
    else{
      $command=<<END;
find Sequence
follow DNA
show -a -T -f $filename
quit
END

    }
    return ($command,$filename);

}

#######################################################################
# Strip unwanted lines from the acefiles                              #
#######################################################################

sub strip_object {
    my ($filename,@deletes) = @_;
    my $strip = "";
    local (*OUT,*FILE);

    foreach $strip (@deletes) {

	open (OUT, ">${filename}.rm");
	open (FILE, "<$filename") || warn "can't do this chappie\n";
	while (<FILE>) {
	    next if (/$strip/);
	    print OUT "$_";
	}
	close FILE;
	close OUT;
	system ("mv -f ${filename}.rm $filename");
    }
}

sub strip_ghosts {
    my ($filename,$object) = @_;
    my ($object_header) = "";
    my $object_look = 0;

    local (*OUT,*FILE);

    open (OUT, ">${filename}.rm");
    open (FILE, "<$filename") || warn "can't do this either chappie\n";
    while (<FILE>) { 
	chomp;
	if (/^$object \: \"\S+.+\"$/) {
	    $object_header = $_;
	    $object_look = 1;
	    next;
	}
	
	if (($object_look == 1) && ($_ ne "")) {
	    print OUT "$object_header\n";
	    print OUT "$_\n";
	    $object_look = 0;
	    next;
	}
	elsif (($object_look == 1) && ($_ eq "")) {
	    $object_header = "";
	    $object_look = 0;
	    next;
	}

	print OUT "$_\n";
    }
    close FILE;
    close OUT;
    system ("mv -f ${filename}.rm $filename");
}



########################################################################
# Process ace files to change 'wormpub' in timestamps to database name #
########################################################################

sub process_ace_file{
    my $filename = shift;
    my $database = shift;

    
    open (INPUT, "<$filename") || die "Couldn't open $filename for reading\n";
    open (OUTPUT, ">${filename}.tmpfile") || die "Couldn't write to tmpfile\n";

    
    while (<INPUT>) {

      if (/\"\d{4}-\d{2}-\d{2}_\d{2}:\d{2}:.*?_.*?\"/){
	s/(\"\d{4}-\d{2}-\d{2}_\d{2}:\d{2}:.*?)_.*?\"/$1_$database\"/g;
	print OUTPUT "$_";
      }
      else{
	print OUTPUT "$_";
      }
    }
    close(INPUT);
    close(OUTPUT);
    system ("mv -f ${filename}.tmpfile $filename");
}


#######################################################################
# Help and error trap outputs                                         #
#######################################################################
 


sub usage {
    my $error = shift;
    
    if ($error == 1) {
	# No command-line options
	print "\nNo command line options given\n\n";
	print "Usage: make_acefiles [-options]\n";
	print " -n : Write .ace files\n";
	print " -e : Write non-elegans nematode EST data sets\n";
	exit(0);
    }
    elsif($error == 2) {
	# Single database mode: database name not recognised
	print "\nNo database of this name exists ('$opt_s')\n";
	print "Check database names in config file\n\n";
	exit(0);
    }
    elsif ($error == 0) {
	# Normal help menu
	system ('perldoc',$0);
	exit (0);
    }
}

__END__



