import matplotlib.pyplot as plt
import numpy as np
import sys

if len( sys.argv ) != 3:
	print( "Error: expected 2 arguments" )
	print( f"{sys.argv[ 0 ]} input_data_file.txt output_graph_file.png" )
	sys.exit( 1 )

input_data_file = sys.argv[ 1 ]
output_graph_file = sys.argv[ 2 ]

group_list = []
bacteria_list = []
gc_list = []

# Read the results from the txt file
with open( input_data_file, "r" ) as file:
    for line in file:
        line = line.strip()

        if not line:
            continue

        group_name, bacteria_name, gc_value = line.split( "," )

        group_list.append( group_name.strip() )
        bacteria_list.append( bacteria_name.strip().replace( "_", " " ) )
        gc_list.append( float( gc_value.strip() ) )
        print( f"name = {bacteria_name} : GC = {gc_value}" )

x = np.arange( 1, len( bacteria_list ) + 1 )
median = np.median( gc_list )

plt.figure( figsize = ( 14, 7 ) )

group_colors = {
    "Thermophiles": "red",
    "Mesophiles": "green",
    "Psychrophiles": "blue"
}

for group, color in group_colors.items():

    indices = [
        i for i, g in enumerate( group_list )
        if g == group
    ]
    # Points
    plt.scatter( 
        [ x[ i ] for i in indices ], 
        [ gc_list[ i ] for i in indices ],
        color = color,
        s = 80,
        label = group
    )

# Median
plt.axhline(
    median,
    color = "black",
    linestyle = "--" ,
    linewidth = 3,
    label = f"Median = {median:.2f}%"
)

plt.xticks(
    x,
    bacteria_list,
    rotation = 45,
    ha = "right"
)

plt.xlabel(
	"Bacteria name",
	fontsize = 15,
	fontweight = "bold"
)

plt.ylabel(
	"GC bases ratio (%)",
	fontsize = 15,
    fontweight = "bold"
)

plt.title(
	"GC content of all thermal groups",
	fontsize = 20,
    fontweight = "bold"
)

plt.legend()
plt.tight_layout()

plt.savefig( output_graph_file, dpi=300, bbox_inches="tight" )
plt.close()

