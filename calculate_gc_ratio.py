import sys

# add below your own function(s) and global variables.
A_count = 0
T_count = 0
G_count = 0
C_count = 0

def count_bases( dna ):
    global A_count
    global T_count
    global G_count
    global C_count

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

    
# read program arguments.
dna_filename = sys.argv[ 1 ] # the very first argument of our program is DNA filename.

# open file.
dna_file = open( dna_filename, "rt" ) # rt to work with text
print( f"DNA file: {dna_filename}" )

# read huge DNA file on chunks.
dna_chunk = dna_file.readline()
while len( dna_chunk ) > 0:
    
    if dna_chunk[ 0 ] == ">":
        print( f"This is the very first line of the fasta file." )
        print( f"DNA chunk {len( dna_chunk )} bytes : {dna_chunk}" )
        dna_chunk = dna_file.readline()
        continue

    # call here your function that operates over DNA data.
    count_bases( dna_chunk )

    dna_chunk = dna_file.readline()


dna_file.close()

# Put here your formula.
GC_percent = 100 * ( ( G_count + C_count ) / ( A_count + T_count + G_count + C_count ) )

print( f"GC% in {dna_filename} is {GC_percent} %" )

