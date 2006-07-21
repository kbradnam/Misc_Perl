#!/usr/bin/perl -w

use strict;
use MIME::Lite;

# 1) Get IP address
my $ip = `curl -s http://checkip.dyndns.org`;
$ip =~ s/.*?(\d+\.\d+\.\d+\.\d+).*/$1/s;


# 4) Create screenshot
system("sudo /usr/sbin/screencapture -m /tmp/screengrab.png");

# 5) convert to jpg.
system("/usr/bin/sips -s format jpeg -s formatOptions 25% --resampleWidth 1024 /tmp/screengrab.png --out /tmp/screengrab.jpg");

### Create a new multipart message:

my $msg = MIME::Lite->new( 
			   From    =>'kbradnam@mac.com',
			   To      =>'kbradnam@mac.com',
			   Subject =>"TEST IP: $ip\n",
			   Type    =>'image/jpg',
			   Encoding => 'base64',
			   Path => '/tmp/screengrab.jpg');
        

### use Net:SMTP to do the sending
$msg->send('smtp','localhost', Debug=>0 );



# tidy up
unlink '/tmp/screengrab.jpg','/tmp/screengrab.png';

exit;
