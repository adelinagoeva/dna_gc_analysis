import pandas as pd
import glob
import os

files = glob.glob("GC_data/*.fna")

results = []

for file in files:

    # Read the file
    data = pd.read_csv(file)

    # Calculate average GC%
    avg_gc = data["gc_percent"].mean()

    # Get bacterium name from filename
    bacterium = os.path.basename(file).replace(".fna", "")

    # Save result
    results.append({
        "bacterium": bacterium,
        "avg_gc": avg_gc
    })

print(results)
