#!/usr/bin/perl -w
#
# Birthday reminder script
# by Keith Bradnam

use strict;
use Net::SMTP;

my @date = localtime;
my $this_day = $date[3];
my $this_month = $date[4];
$this_month++;

my $text = "Keith's Birthday Reminder\n";
my $path = glob("~keith/Work/Code/my_scripts/birthdays.inf");

open(DATA,"$path") || die "Couldn't open file\n";

while(<DATA>){
  chomp;

  my ($birth_day,$birth_month,$name) = split(/,/);
  $name =~ s/:/, /g;
  if ($birth_month == $this_month){
    if((($birth_day - $this_day) < 7) && (($birth_day - $this_day) > 0) ){
      $text .=  "$birth_day/$birth_month\t$name \n"; 
    }
  }
  # check for birthdays in next month if near the end of the present month
  if ($this_day > 24){
    if (($birth_month == $this_month+1) && ($birth_day < 7)){
      $text .=  "$birth_day/$birth_month\t$name \n"; 	
    } 
  }
}
close(DATA);



# send the actual mail if there is data to send
if($text){

    # set up mail connection
#    my $mail_server = "smtp.ucdavis.edu";
    
    # Connect to the server 
    my $smtp = Net::SMTP->new("smtp.ucdavis.edu"); 
    die "Couldn't connect to server" unless $smtp;
    
    my $MailFrom = "kbradnam\@mac.com"; 
    my $MailTo   = "kbradnam\@mac.com"; 
    
    $smtp->mail($MailFrom); 
    $smtp->to($MailTo ); 

    $smtp->data();
    $smtp->datasend("Subject: Birthday reminder\n");		   
    $smtp->datasend("$text\n");
    $smtp->dataend();
    $smtp->quit();
}
