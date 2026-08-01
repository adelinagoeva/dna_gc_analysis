################################################################################
# slide_window_of_file : 21.JUL.2026 : Read DNA.fasta using Slide Window algo  #
# ====================                                                         #
# Implement Slide Window algorithm for fasta DNA files.                        #
#                                                                              #
#                                                                              #
################################################################################

import sys

from typing import NamedTuple

# Counter for bases in calculated slice
class BasesCount( NamedTuple ):
    A_cnt: int
    T_cnt: int
    G_cnt: int
    C_cnt: int

# List of calculated bases per slide
records_list = []

# add below your own function(s).
def count_bases( dna ) -> BasesCount:
    A_count = 0
    T_count = 0
    G_count = 0
    C_count = 0

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

    return BasesCount( A_cnt = A_count, T_cnt = T_count, G_cnt = G_count, C_cnt = C_count )


# check input arguments.
if len( sys.argv ) != 4:
    print( f"Incorrect number of arguments!\nUsage:\n{sys.argv[ 0 ]} DNA_filename Window_size Slide_length" )
    sys.exit( 1 )

# read program arguments.
dna_filename = sys.argv[ 1 ] # the 1st argument of our program is DNA filename.
try:
    window_size  = int( sys.argv[ 2 ] ) # the 2nd argument of our program is the size of the Window.
    slide_length = int( sys.argv[ 3 ] ) # the 3rd argument of our program is the length of the Slide.

except ValueError:
    print( "ERROR: the argument is not a number\nUsage:\n{sys.argv[ 0 ]} DNA_filename Window_size Slide_length" )
    sys.exit( 1 )


print( f"DNA filename : {dna_filename}" )
print( f"Window size  : {window_size}" )
print( f"Slide length : {slide_length}" )

if 2 * slide_length != window_size:
    print( "Window size shall be twice larger than a Slde length." )
    sys.exit( 2 )

# open file.
try:
    dna_file = open( dna_filename, "rt" ) # rt to work with text
    print( f"Successfully opened file: {dna_filename}" )

except FileNotFoundError:
    print( f"ERROR: file {dna_filename} does not exsist." )
    sys.exit( 2 )

except PermissionError:
    print( f"ERROR: insufficient permisions to open file {dna_filename} or directory to that file." )
    sys.exit( 2 )

except Exception as e:
    print( f"ERROR: unexpected error with file {dna_filename} - {e}." )
    sys.exit( 2 )

# read the header line of the huge DNA file
# and validate it.
dna_chunk = dna_file.readline()
if len( dna_chunk ) == 0:
    print( "Zero size file." )
    sys.exit( 3 )

if dna_chunk[ 0 ] != ">":
    print( "Wrong header format." )
    sys.exit( 3 )

not_last_chunk_flag = True
while not_last_chunk_flag:
    # read twice as much data as the window size,
    # then clean up from new lines.
    dna_chunk = dna_file.read( 30 )
    if len( dna_chunk ) < 30:
        not_last_chunk_flag = False
        print( "Last chuk read." )

    dna_chunk = dna_chunk.replace("\r", "").replace("\n", "")
    chunk_length = len( dna_chunk )
    print( f"Cleaned chunk has length = {chunk_length}" )

    position = 0
    while position < chunk_length:
        if position + window_size > chunk_length:
            print( f"Leftover is smaller : position = {position}" )
            data_chunk_one = dna_chunk[ position : position + int( chunk_length / 2 ) ]
            position += int( chunk_length / 2 )
            data_chunk_two = dna_chunk[ position : ] # get bytes from position to the end of the dna_chunk.
            position = chunk_length
        else:
            print( f"Regular size of slide : position = {position}" )
            data_chunk_one = dna_chunk[ position : position + slide_length ]
            position += slide_length
            data_chunk_two = dna_chunk[ position : position + slide_length ]
            position += slide_length

        measurement = count_bases( data_chunk_one )
        records_list.append( measurement )

        measurement = count_bases( data_chunk_two )
        records_list.append( measurement )

dna_file.close()

print( f"Number of calculated slides: {len( records_list ) }" )

print( records_list )

