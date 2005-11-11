#!/usr/local/bin/perl5.6.0 -w

use Ace;


$|=1;

# open a local database connection
$db = Ace->connect(-path  =>  '/nfs/disk100/wormpub/DATABASES/current_DB') || die "Can't connect\n";


$person = $db->fetch(Person => "WBPerson1971");
print "Person is $person\n";
($lab) = $person->at("Address.Street_address");
print "Lab is $lab\n";

($country) = $person->at("Address.Country");
print "Country is $country\n";

@address = $person->at("Address")->col;
print "Address is @address\n\n";

@street =  $person->at("Address.$address[0]")->col;

print "@street\n\n";

$db->close;
exit;





