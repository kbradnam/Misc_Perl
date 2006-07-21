#!/usr/bin/perl -w
#
# get_IP.pl
#
# a script to email your IP address and include a screenshot
#
# by Keith Bradnam
# July 2006
#
###########################################################

use strict;
use MIME::Lite;

# NOTES
# i)  You must install MIME::Lite manually, it is not installed by default
# ii) You must enable the admin group to run the screencapture command without
#     requiring a password (use the visudo command to do this)


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


# 4) Create screenshot
system("sudo /usr/sbin/screencapture -m /tmp/screengrab.png");


# 5) convert to jpg.

# set jpg compression level (1=lowest quality, 100=highest quality)
my $compression = 25;

system("/usr/bin/sips -s format jpeg -s formatOptions ${compression}% --resampleWidth 1024 /tmp/screengrab.png --out /tmp/screengrab.jpg");


# 6) create email

my $msg = MIME::Lite->new( 
			   From    =>'kbradnam@mac.com',
			   To      =>'kbradnam@mac.com',
			   Subject =>"MyIP: $machine - $ip - $place\n",
			   Type    =>'image/jpg',
			   Encoding => 'base64',
			   Path => '/tmp/screengrab.jpg');
        

# 7) send email
$msg->send('smtp','localhost', Debug=>0 );


# 8) tidy up
unlink '/tmp/aaa.jpg','/tmp/aaa.png';

exit;
