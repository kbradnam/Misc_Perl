#!/usr/bin/perl -w
#
# get_isight.pl
#
# a script to email a picture from your isight camera
#
# by Keith Bradnam
# July 2006
#
###########################################################

use strict;
use MIME::Lite;

# NOTES
# i)  You must install MIME::Lite manually to send email attachments, it is not installed by default
# ii) You must also install the command line tool 'isightcapture' (http://www.intergalactic.de/hacks.html)


# 1) Find out which machine is this being run
my $machine = `uname -n`;
$machine =~ s/\.local$//;


# 2) Take picture
system("/usr/bin/isightcapture /tmp/isight.jpg");

# 3) create email
my $msg = MIME::Lite->new( 
			   From    =>'kbradnam@mac.com',
			   To      =>'kbradnam@mac.com',
			   Subject =>"ISIGHT: $machine\n",
			   Type    =>'image/jpg',
			   Encoding => 'base64',
			   Path => '/tmp/isight.jpg');
        
# 5) send email
$msg->send('smtp','localhost', Debug=>0 );


# 6) tidy up
unlink '/tmp/isight.jpg';

exit;
