#!/usr/bin/perl -w
#
# get_desktop.pl
#
# a script to email a screenshot of your computers desktop
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

# 1) Find out which machine is this being run on?
my $machine = `uname -n`;
$machine =~ s/\.local$//;
chomp($machine);

# 2) Create screenshot
# -x to turn sound off
# -m to capture main monitor
system("sudo /usr/sbin/screencapture -x -m /tmp/screengrab.png");


# 3) convert to jpg.
# set jpg compression level (1=lowest quality, 100=highest quality)
my $compression = 25;
system("/usr/bin/sips -s format jpeg -s formatOptions ${compression}% --resampleWidth 1024 /tmp/screengrab.png --out /tmp/screengrab.jpg");

# 4) create email

my $msg = MIME::Lite->new( 
			   From    =>'keithwho@mac.com',
			   To      =>'keithwho@mac.com',
			   Subject =>"SCREENSHOT: $machine\n",
			   Type    =>'image/jpg',
			   Encoding => 'base64',
			   Path => '/tmp/screengrab.jpg');
        
# 5) send email
$msg->send('smtp','localhost', Debug=>0 );


# 6) tidy up
unlink '/tmp/screengrab.jpg','/tmp/screengrab.png';

exit;
