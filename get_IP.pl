#!/usr/bin/perl -w

use strict;
use MIME::Lite;

# 1) Get IP address
my $ip = `curl -s http://checkip.dyndns.org`;
$ip =~ s/.*?(\d+\.\d+\.\d+\.\d+).*/$1/s;

# 2) Look up where that IP address corresponds to
my $address = `curl -s http://www.antionline.com/tools-and-toys/ip-locate/`;
$address =~ m/are located in (.* United States)/;
my $place = $1;


# 3) Find out which machine is this being run on?
my $machine = `uname -n`;
$machine =~ s/\.local$//;

#print "$machine\n";
#print "$place\n";
#print "IP is $ip\n\n";


# 4) Create screenshot
`screencapture /tmp/aaa.png`;

# 5) convert to jpg.
`sips -s format jpeg -s formatOptions 25% --resampleWidth 1024 /tmp/aaa.png --out /tmp/aaa.jpg`;

### Create a new multipart message:

my $date = `date`;

my $msg = MIME::Lite->new( 
    From    =>'kbradnam@mac.com',
    To      =>'kbradnam@mac.com',
    Subject =>"MyIP: $machine - $ip - $place\n",
    Type    =>'multipart/mixed');
        
### Add parts (each "attach" has same arguments as "new"):
$msg->attach(Type     =>'TEXT',   
             Data     =>"Sent on $date\n"
             );  

# Only attach picture if on mosquito
#if($machine eq "Mosquito"){
    $msg->attach(Type     =>'image/gif',
		 Path     =>'/tmp/aaa.jpg',
		 Filename =>'aaa.jpg',
		 Disposition => 'attachment'
		 );
#}

### use Net:SMTP to do the sending
$msg->send('smtp','localhost', Debug=>0 );



# tidy up
#unlink '/tmp/aaa.jpg','/tmp/aaa.png';

exit;
