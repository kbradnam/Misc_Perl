#!/usr/bin/perl -w

use strict;
use Net::SMTP;

my $mailprog = '/usr/sbin/sendmail';
my $ip = `curl -s http://checkip.dyndns.org`;
$ip =~ s/.*?(\d+\.\d+\.\d+\.\d+).*/$1/s;

my $address = `curl -s http://www.antionline.com/tools-and-toys/ip-locate/`;
$address =~ m/are located in (.* United States)/;
my $place = $1;

my $date = `date`;

# which machine is this being run on?
my $machine = `uname -n`;
$machine =~ s/\.local$//;

#print "$machine\n";
#print "$place\n";
#print "IP is $ip\n\n";


my $smtp = Net::SMTP->new("localhost");

die "Couldn't connect to SMTP server\n" unless $smtp;

my $from = "kbradnam\@mac.com";
my $to = "kbradnam\@mac.com";

$smtp->mail($from);
$smtp->to($to);
$smtp->data();
$smtp->datasend("Subject: MyIP: $machine - $ip - $place\n");
$smtp->datasend("Sent on $date\n");
$smtp->dataend();
$smtp->quit();


exit;
