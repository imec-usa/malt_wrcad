#!/usr/bin/perl -w

# Search each file for lines with unmatched *alloc calls and print them out.

my %allocs = ();
my %frees = ();

while (<>) {
  if (/(?<alloc>\bi?vector|\bmatrix|resprintf|strdup|(?:m|re|c)alloc)\s*\(/) {
    if (m{// mem:([a-z]+)$}) {
      if (defined $allocs{$1}) {
        my ($file, $line_no, $_line) = @{$allocs{$1}};
        print STDERR "(warning) duplicate label '$1': $ARGV:$. (originally defined at $file:$line_no)\n";
      } else {
        $allocs{$1} = [$ARGV, $., $_];
      }
    } else {
      print STDERR "(warning) unlabeled `$+{alloc}`: $ARGV:$. $_";
    }
  }
  if (m{free\w*\(.*// mem:([a-z]+)$}) {
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
