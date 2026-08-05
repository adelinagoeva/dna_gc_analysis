################################################################################
# slide_window_of_file : 21.JUL.2026 : Read DNA.fasta using Slide Window algo  #
# ====================                                                         #
# Implement Slide Window algorithm for fasta DNA files.                        #
#                                                                              #
#                                                                              #
################################################################################

import os
import sys
import matplotlib.pyplot as plt

from typing import NamedTuple
from datetime import datetime

# Counter for bases in calculated slice
class BasesCount( NamedTuple ):
    A_cnt: int
    T_cnt: int
    G_cnt: int
    C_cnt: int
    GC_ratio: float

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

    total_bases = A_count + T_count + G_count + C_count
    if total_bases > 0:
        GC_percent = 100 * ( ( G_count + C_count ) / ( A_count + T_count + G_count + C_count ) )
    else:
        print( f"count_bases : A = {A_count}, T = {T_count}, G = {G_count}, C = {C_count} : GC_percent will be set as 0.0" )
        GC_percent = 0.0

    return BasesCount( A_cnt = A_count, T_cnt = T_count, G_cnt = G_count, C_cnt = C_count, GC_ratio = GC_percent )


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


start_time = datetime.now()
print( f"{start_time:%Y-%m-%d %H:%M:%S} : begin of {sys.argv[ 0 ]}" )
print( f"DNA filename : {dna_filename}" )
print( f"Window size  : {window_size}" )
print( f"Slide length : {slide_length}\n" )

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

print( "\nCalculate GC ratio for each slide." )
not_last_chunk_flag = True
fragment = ""
while not_last_chunk_flag:
    # read twice as much data as the window size,
    # then clean up from new lines.
    dna_chunk = dna_file.read( 2 * window_size )
    if len( dna_chunk ) < 2 * window_size:
        not_last_chunk_flag = False
#        print( "Last chunk of DNA file was read." )

    dna_chunk = dna_chunk.replace("\r", "").replace("\n", "")
    dna_chunk = fragment + dna_chunk
    chunk_length = len( dna_chunk )
#    print( f"Cleaned chunk has length = {chunk_length}" )
    if chunk_length == 0:
        break

    position = 0
    fragment = ""
    while position < chunk_length:
        if position + window_size > chunk_length:
#            print( f"Leftover is smaller : position = {position}/{chunk_length}" )
            if not_last_chunk_flag:
                fragment = dna_chunk[ position : ] 
                break

            middle = int( ( position + chunk_length ) / 2 )
            data_chunk_one = dna_chunk[ position : middle ]
            data_chunk_two = dna_chunk[ middle : ] # get bytes from position to the end of the dna_chunk.
            position = chunk_length # no need this assignment.
        else:
#            print( f"Regular size of slide : position = {position}" )
            data_chunk_one = dna_chunk[ position : position + slide_length ]
            position += slide_length
            data_chunk_two = dna_chunk[ position : position + slide_length ]
            position += slide_length

        measurement = count_bases( data_chunk_one )
        records_list.append( measurement )

        measurement = count_bases( data_chunk_two )
        records_list.append( measurement )

dna_file.close()

print( f"{len( records_list ) } slides was calculated from provided fasta file." )
# print( records_list )

print( "\nPrepare a list of non-overlaping windows of GC ratio." )
# each window is made up of two adjacent slides.
windows_avg_gc = []
list_size = len( records_list )
for index in range( 0, list_size, 2 ):
    if index + 1 < list_size:
        slide_one = records_list[ index ]
        slide_two = records_list[ index + 1 ]
        windows_gc_ratio = ( slide_one.GC_ratio + slide_two.GC_ratio ) / 2
        windows_avg_gc.append( windows_gc_ratio )
    else:
        # if our records are odd, the last value cannot be paired. Use it as is.
        windows_avg_gc.append( records_list[ index ].GC_ratio )

print( f"{len( windows_avg_gc )} non-overlaping windows was prepared from provided fasta file." )

print( "\nPrepare MatPlotLib to create an graph from the data in the list." )
x_positions = [ i * window_size for i in range ( len( windows_avg_gc ) ) ]

plt.figure( figsize = ( 10, 5 ) )

# plt.plot( x_positions, windows_avg_gc, color = "teal", linestyle = '-', marker = 'o', linewidth = 2, label = 'GC % per window' )
plt.plot( x_positions, windows_avg_gc, color = "teal", linestyle = '-', marker = None, linewidth = 1, label = 'GC % per window' )

total_avg = sum( windows_avg_gc) / len( windows_avg_gc )
plt.axhline( y = total_avg, color = 'crimson', linestyle = '--', linewidth = 1.5, label =f'Total average ({total_avg:.2f}%)' )

file_name = os.path.basename( dna_filename )
organism_name = file_name.split( '.')[0]

plt.title( f"{organism_name} GC-content changes", fontsize = 14, fontweight = 'bold' )
plt.xlabel( "Position in genome (number of bases)", fontsize = 12 )
plt.ylabel( "Average GC, %", fontsize = 12 )

plt.legend( loc = 'upper right' )

plt.grid( True, linestyle='--', alpha = 0.7 )

plt.tight_layout()

# plt.show()
# get the name of the bacteria from fasta filename.
file_name = os.path.basename( dna_filename )
organism_name = file_name.split( '.')[0]
current_timestamp= datetime.now().strftime( "%Y%m%d_%H%M%S" )
image_name = f"{organism_name}.gc_content.{current_timestamp}.png"

plt.savefig( image_name, dpi = 300, bbox_inches = 'tight' )
print( f"The graph was successfully saved in file '{image_name}' on disk\n" )

print( "Print the average GC ratio for each window:" )
print( "num : data (%)\n==============" )
for i, gc in enumerate( windows_avg_gc, start=1 ):
    print( f"{i:04d}: {gc:.3f}" )

stop_time = datetime.now()
print( f"\n{stop_time:%Y-%m-%d %H:%M:%S} : end of {sys.argv[ 0 ]}" )

