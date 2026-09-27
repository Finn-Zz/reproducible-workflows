use strict;
use MaterialsScript qw(:all);

# Define the base directory
my $cwd = "D:\\MS\\MIL_Files\\Documents";

# Concatenate the subdirectory path
my $input_directory = $cwd . "\\NEW\\stp_test\\all";
# Output directory to save the optimized CIF files
my $output_directory = $cwd . "\\NEW\\stp_test\\opt_all";

# Get all CIF files in the directory
opendir(DIR, $input_directory) or die "Cannot open directory: $!";
my @files = grep(/\.xsd$/, readdir(DIR));
closedir(DIR);

foreach my $file (@files) {
    my ($file_name) = $file =~ /([^\.]+)/;
    my $input_xsd_name = $file_name . ".xsd";
    # Full path to the CIF file
    my $input_file_path = "$input_directory\\$input_xsd_name";
    
    # Step 1: Read the CIF file
    my $doc = Documents->Import($input_file_path);
    if (!$doc) {
        print "Failed to open CIF file: $input_file_path\n";
        next;
    }
    
    #my $num_of_atoms = $doc->UnitCell->Atoms->Count;
    #my $calc_quality
    #if ($num_of_atoms>2000) { $calc_quality = "Fine" }
    #else { $calc_quality = "Ultra-fine" };
    
    # Step 2: Set up the Forcite optimization task with the DREIDING force field
    my $forcite_opt = Modules->Forcite->GeometryOptimization;
    
    # Change Forcite settings
    Modules->Forcite->ChangeSettings([
        OptimizeCell       => 1,
        CurrentForcefield  => "Universal",
        Quality            => "Ultra-fine"
        ,
    ]);
    
    eval {
        # Step 3: Run the optimization
        my $results = $forcite_opt->Run($doc, Settings(OptimizeStructure => "Yes", WriteLevel => "Silent"));
        if (!$results) {
            print "Optimization failed for: $input_file_path\n";
            $doc->Close();
            next;
        }

        # Step 4: Save the optimized structure as a new XSD file
	my $output_xsd_name = $file_name . ".xsd";
	my $optimized_xsd_path = "$output_directory\\$output_xsd_name";
	$doc->Export($optimized_xsd_path);

    };
    if ($@) {
        print "Failed to optimize file: $file_name. Error: $@\n";
        next;
    };
    
    # Step 5: Close the document
    $doc->Close();
    $doc->Delete;
    
}