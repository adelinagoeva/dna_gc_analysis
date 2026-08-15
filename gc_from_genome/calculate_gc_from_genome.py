##import sys
import os

# Import python script and use this big function:
def calculate_gc( dna_filename ):
	# variables.
	A_count = 0
	T_count = 0
	G_count = 0
	C_count = 0

	# inner function.
	def count_bases( dna ):
		nonlocal A_count, T_count, G_count, C_count

		for base in dna:
			match base:
				case "A":
					A_count +=1
				case "T":
					T_count +=1
				case "G":
					G_count +=1
				case "C":
					C_count +=1

	# open file.
	dna_file = open( dna_filename, "rt" ) # rt to work with text
	print( f"Calculate GC ratio from genomic file: {dna_filename}" )

	# read huge DNA file on chunks.
	dna_chunk = dna_file.readline()
	while len( dna_chunk ) > 0:
		if dna_chunk[ 0 ] == ">":
			dna_chunk = dna_file.readline()
			continue

		# operate over DNA data line.
		count_bases( dna_chunk )

		dna_chunk = dna_file.readline()

	dna_file.close()

	# Calculate GC ratio.
	GC_percent = 100 * ( ( G_count + C_count ) / ( A_count + T_count + G_count + C_count ) )

	return( f"{GC_percent:.3f}" )

