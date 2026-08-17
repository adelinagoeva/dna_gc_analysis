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
gc3_list = []

# Read the results from the txt file
with open( input_data_file, "r" ) as file:
    for line in file:
        line = line.strip()

        if not line:
            continue

        group_name, bacteria_name, gc_value, gc3_value = line.split( "," )

        group_list.append( group_name.strip() )
        bacteria_list.append( bacteria_name.strip().replace( "_", " " ) )
        gc_list.append( float( gc_value.strip() ) )
        gc3_list.append( float( gc3_value.strip() ) )
        print( f"name = {bacteria_name} : GC = {gc_value} : GC3 = {gc3_value}" )


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
        [ gc_list[i] for i in indices ],
        [ gc3_list[i] for i in indices ],
        color = color,
        s = 80,
        label = group
    )

    plt.plot(
    	[0, 100],
    	[0, 100],
    	color="black",
    	linestyle="--",
    	linewidth=1
    )

plt.xlabel(
	"GC (%)",
	fontsize = 15,
	fontweight = "bold"
)

plt.ylabel(
	"GC3 (%)",
	fontsize = 15,
    fontweight = "bold"
)

plt.title(
	"GC vs GC3 content of all thermal groups",
	fontsize = 20,
    fontweight = "bold"
)

plt.legend()
plt.tight_layout()

plt.savefig( output_graph_file, dpi=300, bbox_inches="tight" )
plt.close()

