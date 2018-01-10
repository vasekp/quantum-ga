use strict;
use warnings;
use utf8;

my $ma = 0;
my $mb = 0;
my $mc = 0;
my $last = "";

$| = 1;

while(<>) {
  last if(/{0\.000.,0\.000.,2,15,2}/);
  if(/{0\.0...,0\.0...,2,15,2}/ && !$ma) {
    print $last;
    print;
    $ma = 1;
  } elsif(/{0\.00..,0\.00..,/ && !$mb) {
    print $last;
    print;
    $mb = 1;
  } elsif(/{0\.000.,0\.000.,/ && !$mc) {
    print $last;
    print;
    $mc = 1;
  } else {
    last if(/^Gen 3000/);
    s/\n/\r/;
    print;
  }
  $last = $_;
}
print $last;
print;
