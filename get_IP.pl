#!/usr/bin/perl -w
#
# get_IP.pl
#
# a script to email your computers IP address and possible geographic location
#
# by Keith Bradnam
# July 2006
#
################################################################################

use strict;
use Net::SMTP;


# 1) Get IP address
my $ip = `curl -s http://checkip.dyndns.org`;
$ip =~ s/.*?(\d+\.\d+\.\d+\.\d+).*/$1/s;


# 2) Look up where that IP address corresponds to
my $address = `curl -s http://www.antionline.com/tools-and-toys/ip-locate/`;
$address =~ m/are located in (.* United States)/;
my $place = $1;


# 3) Find out which machine is this being run on (remove any .local suffix)
my $machine = `uname -n`;
$machine =~ s/\.local$//;


# 4) Send email using Net::STMTP
my $smtp = Net::SMTP->new("mailhost") || die "Couldn't connect to mailhost\n";

my $from = "kbradnam\@mac.com";
my $to   = "kbradnam\@mac.com";

$smtp->mail($from);
$smtp->to($to);
$smtp->data();
$smtp->datasend("Subject: MyIP: $machine - $ip - $place\n");
$smtp->dataend();
$smtp->quit();


exit;

