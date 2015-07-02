#!/usr/bin/perl -w
#
#   Perl equivalent of the popular dos2unix utility:
#
#   Convert DOS line endings to Unix line endings:
#   works in bulk, safely updates files in place.
#
my  ($filename, $line, $count);
$count = 0;

#   If no arguments, print an error message
if( $#ARGV < 0 ) {
	    print "Usage: $0 filenames\n";
		    print ";Replace DOS line endings with Unix line endings\n";
			    exit(5);
				    }

#   Loop through each given filename
foreach $filename (@ARGV)
{
	    if( -e "$filename.bak" ) {
			        printf "Skipping $filename.bak - it already exists\n";
					    }
		    elsif(!( -f $filename && -r $filename && -w $filename  )) {
				        printf "Skipping $filename - not a regular writable file\n";
						    }
			    else {
					        rename("$filename","$filename.bak");
							        open INPUT, "$filename.bak";
									        open OUTPUT, ">$filename";

											        while( <INPUT> ) {
														            s/\r\n$/\n/;     # convert CR LF to LF
																		            print OUTPUT $_;
																	        }

													        close INPUT;
															        close OUTPUT;
																	        unlink("$filename.bak");
																			        $count++;
																					    }
}
printf "Processed $count files.\n";
