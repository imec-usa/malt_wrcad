#!/usr/bin/perl -w

# Search each file for lines with unmatched *alloc calls and print them out.

my %allocs = ();
my %frees = ();

while (<>) {
  if (/(?:\bi?vector|\bmatrix|resprintf|strdup|(?:m|re|c)alloc)\s*\(/) {
    if (m{// ([a-z]+-[a-z]+-[a-z]+)$}) {
      $allocs{$1} = [$ARGV, $., $_];
    } else {
      print "(unmatched) $ARGV:$. $_";
    }
  }
  if (m{free\w*\(.*// ([a-z]+-[a-z]+-[a-z]+)$}) {
    $frees{$1} = 1;
  }
} continue {
  close ARGV if eof;
}

foreach (keys %allocs) {
  if (!exists $frees{$_}) {
    my($file, $line_no, $line) = @{$allocs{$_}};
    print "$file:$line_no $line";
  }
}
