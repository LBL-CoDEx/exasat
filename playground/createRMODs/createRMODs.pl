#!/usr/bin/perl                                                                        

use strict;                                                                      
use Cwd;

my %not_compiled=();
my %compiled=() ;

#my $filename = $ARGV[$#ARGV];                  
#my $directory="/global/homes/d/dunat/hopper/s3d/source/modules";
my $directory = "/global/homes/d/dunat/hopper/BoxLib/Src/F_BaseLib";
my $rmod_dir = "./rmod";
my $include_dir = "-I$rmod_dir -I/usr/common/usg/openmpi/1.4.5_mom/gcc/include ";
my $options="-fno-range-check";
#if($#ARGV  < 0 || !($directory) )                                                     
#{                                                                                     
#    print " no input file  \n";                                                       
#    exit(-1);                                                                         
#}                                                                                     

sub collectModuleNames()
{
    my @files = <$directory/*.f90>;
    my $file;

    foreach $file(@files)
    {
        my $my_module_name ;

        open(IN, "+<$file");
        my $line;

        while($line = <IN>)
        {
            if($line =~ /^\s*(module|program)\s+(\w+)/)
            {
                if($line !~ /procedure/){
                    $not_compiled{$2} = "1";
                }
            }
        }
    }
}

sub compileModules()
{
    my @files = <$directory/*.f90>;
    my $file;
    my $at_least_one = 1;

    while( (keys %not_compiled ne 0) and ($at_least_one eq  1))
    {
        $at_least_one = 0;

        foreach $file(@files)
	{
            my @my_module_name ;

            open(IN, "+<$file");
            my $line;

            my $found = 0;
            my $num_modules = 0;

            while($line = <IN>)
            {
                if($line =~ /^\s*(module|program)\s+(\w+)/)
		{
                    if($line !~ /procedure/)
                    {
			if(!$compiled{$2}){
			    push(@my_module_name, $2);
			}
                    }
                }
                if($line =~ /^\s*(use)\s+(\w+)/ and scalar(@my_module_name)> 0)
                {
                    my $module_name = $2;
                    #search this module in the modules we already compiled             

                    $num_modules = $num_modules + 1;

                    if($compiled{$module_name})
                    {
			$found = $found + 1;
                    }
                }
            } # end of while                                                           

	    if ( ($num_modules eq $found) and scalar(@my_module_name) > 0  ) # ready to compile  
	    {
		system("./identityTrans -rose:skipfinalCompileStep $include_dir $options $file");
		system("mv *.mod $rmod_dir");
		system("mv *.rmod $rmod_dir");

		foreach my $my_mod(@my_module_name)
		{
		    print "$file  $my_mod \n"; 
		    # compile and move that to compiled list                               
		    $compiled{$my_mod} = "1";
		    
		    my $numkeys =  keys(%not_compiled);
		    delete($not_compiled{$my_mod});
		    
		    if($numkeys ne keys(%not_compiled))
		    {
			$at_least_one =1;
		    }
		}
		
	    }
	} # end of files                                                               
    }# end of while                                                                    
}

system("mkdir -p $rmod_dir");

collectModuleNames();

print keys(%not_compiled)."\n";

for my $key( keys %not_compiled)
{
#    print $key."\n";
}

compileModules();

print keys(%compiled)."\n";

for my $key( keys %compiled)
{
 #   print $key."\n";
}
