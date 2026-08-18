import sys
import os

import calculate_gc_from_genome

if len( sys.argv ) < 4:
	print( "Error: expected at least 3 arguments" )
	print( f"{sys.argv[ 0 ]} sort/nosort calculated_gc.txt bacteria_1.genomic.fna [bacteria_2.genomic.fna, bacteria_3.genomic.fna, ...]" )
	sys.exit( 1 )

sort_flag = sys.argv[ 1 ]
store_file = sys.argv[ 2 ]
genome_list = sys.argv[ 3: ]

bacteria_gc_list = []

for dna_filename in genome_list:
	group_name = os.path.basename( os.path.dirname( dna_filename ) )

	filename_only = os.path.basename( dna_filename )
	bacteria_name = filename_only.split( "." )[ 0 ]
	bacteria_name = bacteria_name.strip().replace( "_", " " )

	gc_ratio = calculate_gc_from_genome.calculate_gc( dna_filename )
	result = f"{group_name}, {bacteria_name}, {gc_ratio}"

	bacteria_gc_list.append( result )

if sort_flag == "sort":
	bacteria_gc_list.sort( key = lambda item: float( item.split( "," )[ -1 ] ) )

##print( bacteria_gc_list )
with open( store_file, "w" ) as file:
	for item in bacteria_gc_list:
		file.write( item + "\n" )

