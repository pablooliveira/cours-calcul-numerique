# Support for pythontex in v. 0.16 or higher, with latexmk 4.62 or higher
#
# What these definitions provide/do:
# 1. Variable specifying command string for invoking pythontex
# 2. Addition to %extra_rule_spec of template for pythontex rule. This
#    tells latexmkrc to create the rule when it is initializing for
#    processing a TeX file.
# 3. A subroutine mypythontex that the pythontex rule is defined to
#    call. This runs pythontex and then sets dependency information.
# 4. Settings for the files generated by the pythontex package and the
#    pythontex program so that the files are deleted in a clean-up
#    operation.

$clean_ext .= " pythontex-files-%R/* pythontex-files-%R";
push @generated_exts, 'pytxcode';

$pdf_mode = 5;

$pythontex = 'pythontex --interpreter python:python3 %O %S';
$extra_rule_spec{'pythontex'}  = [ 'internal', '', 'mypythontex', "%Y%R.pytxcode",  "%Ypythontex-files-%R/%R.pytxmcr",    "%R", 1 ];

sub mypythontex {
   my $result_dir = $aux_dir1."pythontex-files-$$Pbase";
   my $ret = Run_subst( $pythontex, 2 );
   rdb_add_generated( glob "$result_dir/*" );
   my $fh = new FileHandle $$Pdest, "r";
   if ($fh) {
      while (<$fh>) {
         if ( /^%PythonTeX dependency:\s+'([^']+)';/ ) {
	     print "Found pythontex dependency '$1'\n";
             rdb_ensure_file( $rule, $aux_dir1.$1 );
	 }
      }
      undef $fh;
   }
   else {
       warn "mypythontex: I could not read '$$Pdest'\n",
            "  to check dependencies\n";
   }
   return $ret;
}
