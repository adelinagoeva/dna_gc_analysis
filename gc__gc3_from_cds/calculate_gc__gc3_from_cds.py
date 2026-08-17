from Bio import SeqIO

def calculate_gc3( cds_filename ):
	total_gc3_g_count = 0
	total_gc3_c_count = 0
	total_gc3_bases_count = 0

	total_g_count = 0
	total_c_count = 0
	total_bases_count = 0

	print( f"Calculate GC & GC3 ratio for CDS file: {cds_filename}" )
	for record in SeqIO.parse( cds_filename, "fasta" ):
		if "[pseudo=true]" in record.description:
			continue

		# requires the sequence length to be a multiple of 3 (coding DNA)
		if len( record.seq ) % 3 != 0:
			continue

		sequence = record.seq
		gc3_seq = sequence[ 2::3 ]

		total_gc3_g_count += gc3_seq.count( "G" )
		total_gc3_c_count += gc3_seq.count( "C" )
		total_gc3_bases_count += len( gc3_seq )

		total_g_count += sequence.count( "G" )
		total_c_count += sequence.count( "C" )
		total_bases_count += len( sequence )


	total_gc3_ratio = (
		100 * ( total_gc3_g_count + total_gc3_c_count )
		/ total_gc3_bases_count
	)

	total_gc_ratio = (
		100 * ( total_g_count + total_c_count )
		/ total_bases_count
	)

	return( f"{total_gc_ratio:.3f}, {total_gc3_ratio:.3f}" ) 

