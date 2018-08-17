global env;
set outfile [open charmmdir.txt w];
puts $outfile "$env(CHARMMPARDIR)";
puts $outfile "$env(CHARMMTOPDIR)";
close $outfile
exit;