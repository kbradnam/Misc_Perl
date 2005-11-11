#!/usr/local/bin/perl

use Socket;

my $querycontent="%5Blibs%3D%7Bembl_SP_emblnew%7D-Organism%3A%22Caenorhabditis%20elegans%22%5D%26%28%5Blibs-FtQualifier%3Agene%5D%3Eparent%29%26%5Blibs-AllText%3Acosmid%5D+-f+div+-f+key+-f+org+-f+prd+-f+ftd+-f+ftl+-f+date";

$querycontent = "-e+[{swall}-Organism:%20\"Caenorhabditis%20elegans\"-AccNumber:%20Q22242]";

#$querycontent = "-e+[swall-AccNumber:%20Q222*]%20&%20[swall-Organism:%20%22Caenorhabditis%20elegans%22]";

$querycontent ="-e+[embl-Organism:Caenorhabditis%20elegans*]";

my $request = "/srs6135bin/cgi-bin/wgetz?$querycontent";

my $server = "www.sanger.ac.uk";
if (!defined(open_TCP(*F,$server,80))) {
        print "Error connecting to server at \n";
        exit(-1);
        }
print F "GET $request HTTP/1.0\n\n";
print F "Accept: */*\n";
print F "User-Agent: socketsrs/1.0\n\n";

my $GENE ; my $GENE_1 ; my $SPT ; my $SWISS; my %SPTREMBL ; my %SWISS_PROT ; my %PROTEIN;


# Parsing annotation

my $genes_parsed=0;

while (my $return_line=<F>) {

 print "$return_line";
 if ($return_line =~ /^<pre>ID/){
   $return_line =~ s/<pre>//;
   print "\n$return_line";
 }
 if ($return_line =~ /^(DE|DT|AC|GN|DR|KW)/) {
    print "$return_line";
  }
}
close F;



sub open_TCP 
{
        my ($FS,$dest,$port) = @_;
        my $proto = getprotobyname ('tcp');
        socket ($FS,PF_INET,SOCK_STREAM,$proto);
        my $sin = sockaddr_in($port,inet_aton($dest));
        connect($FS,$sin) || return undef;
        my $old_fh = select($FS);
        $| = 1;
        select($old_fh);

}
