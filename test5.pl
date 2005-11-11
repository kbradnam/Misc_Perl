#!/usr/local/bin/perl -w

($date1,$date2) = get_wormbase_release_date("both");

#print "Date is $date\n\n";
print "Dates are $date1, $date2\n";
sub get_wormbase_release_date{

  my $format = shift;
  print "Format is *$format*\n\n";

  if(!(defined($format))){
    $format = "long";
  }
  elsif($format eq "short"){
    $format = "short";
  }
  elsif($format eq "both"){
    $format = "both";
  }
  else{$format = "long";}

  print "Format is *$format*\n\n";
  
  my $line = `ls -l /nfs/disk69/ftp/pub/wormbase/current_release/letter.WS??`;
  my @split = split(/\s+/,$line);

  my $month = $split[5];
  my $month2;

  if($month eq "Jan"){ $month = "January";   $month2 = "01";}
  if($month eq "Feb"){ $month = "February";  $month2 = "02";}
  if($month eq "Mar"){ $month = "March";     $month2 = "03";}
  if($month eq "Apr"){ $month = "April";     $month2 = "04";}
  if($month eq "May"){ $month = "May";       $month2 = "05";}
  if($month eq "Jun"){ $month = "June";      $month2 = "06";}
  if($month eq "Jul"){ $month = "July";      $month2 = "07";}
  if($month eq "Aug"){ $month = "August";    $month2 = "08";}
  if($month eq "Sep"){ $month = "September"; $month2 = "09";}
  if($month eq "Oct"){ $month = "October";   $month2 = "10";}
  if($month eq "Nov"){ $month = "November";  $month2 = "11";}
  if($month eq "Dec"){ $month = "December";  $month2 = "12";}

  my $day = $split[6];

  my $day2;
  if (length($day) == 1){$day2 = "0".$day;}
  else{$day2 = $day;}

  if    ($day eq "1") {$day .= "st";}
  elsif ($day eq "2") {$day .= "nd";}
  elsif ($day eq "3") {$day .= "rd";}
  elsif ($day eq "21"){$day .= "st";}
  elsif ($day eq "22"){$day .= "nd";}
  elsif ($day eq "23"){$day .= "rd";}
  elsif ($day eq "31"){$day .= "st";}
  else                {$day .= "th";} 



  my $year = `date`;
  $year = substr($year,-3,2);

  # make a text style date
  my $date = $day." ".$month;

  # make a regular xx/xx/xx date
  my $date2 = $day2."/".$month2."/".$year;

  return($date)        if ($format eq "long");
  return($date2)       if ($format eq "short");
  return($date2,$date) if ($format eq "both"); 
}
