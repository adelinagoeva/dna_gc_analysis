import glob
import os
import readline
import subprocess
import sys
from pathlib import Path

# --- Directory for results ---
RESULTS_DIR = Path( "results" )

# --- Path to scripts ---
# Option 1: GC from genome
SCRIPTS_DIR_1 = Path( "gc_from_genome" )
CALC_SCRIPT_1 = SCRIPTS_DIR_1 / "calculate_gc_from_genome_of_bacterial_groups.py"
PLOT_SCRIPT_1 = SCRIPTS_DIR_1 / "scatter_gc_from_genome_of_bacterial_groups.py"

# Option 2: GC / GC3 from CDS
SCRIPTS_DIR_2 = Path( "gc__gc3_from_cds" )
CALC_SCRIPT_2 = (
    SCRIPTS_DIR_2 / "calculate_gc__gc3_from_cds_of_bacterial_groups.py"
)
PLOT_SCRIPTS_2 = {
    "1": (
        "Scatter plot for GC3",
        SCRIPTS_DIR_2 / "scatter_gc3_from_cds_of_bacterial_groups.py",
        "scatter_gc3.png",
    ),
    "2": (
        "Scatter plot for both GC and GC3 (at the same place)",
        SCRIPTS_DIR_2 / "scatter_gc__gc3_from_cds_of_bacterial_groups.py",
        "scatter_gc__gc3.png",
    ),
    "3": (
        "Scatter plot for GC versus GC3 (GC to GC3)",
        SCRIPTS_DIR_2 / "scatter_gc_to_gc3_from_cds_of_bacterial_groups.py",
        "scatter_gc_to_gc3.png",
    ),
}


def ensure_results_dir():
    """Create 'results' directory if not exists."""
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)


def path_completer(text, state):
    """File and directory autocompletion (TAB)."""
    expanded_text = os.path.expanduser(text)
    matches = glob.glob(expanded_text + "*")

    results = []
    for match in matches:
        if os.path.isdir(match):
            results.append(match + "/")
        else:
            results.append(match)

    try:
        return results[state]
    except IndexError:
        return None


def setup_readline():
    """Activate TAB autocompletion."""
    readline.set_completer(path_completer)
    readline.parse_and_bind("tab: complete")
    readline.set_completer_delims(" \t\n\"'")


def print_header():
    print( "=" * 63 )
    print( "    BIOINFORMATICS TOOL FOR GC AND GC3 ANALYSIS OF GENOMES    " )
    print( "=" * 63 )


def check_scripts():
    """Check if all needed scripts exist."""
    """Check for programming scripts lability."""
    missing = []

    # Check for Option 1
    if not CALC_SCRIPT_1.exists():
        missing.append( str( CALC_SCRIPT_1 ) )
    if not PLOT_SCRIPT_1.exists():
        missing.append( str( PLOT_SCRIPT_1 ) )

    # Check for Option 2
    if not CALC_SCRIPT_2.exists():
        missing.append( str( CALC_SCRIPT_2 ) )
    for _, script_path, _ in PLOT_SCRIPTS_2.values():
        if not script_path.exists():
            missing.append( str( script_path ) )

    if missing:
        print(
            "\n[ERROR] Missing following python scripts and subdirectories:",
            file=sys.stderr,
        )
        for m in missing:
            print( f"  - {m}", file=sys.stderr )
        sys.exit( 1 )


def get_genome_files():
    """Interactive collecting information about .cds.fna and .genomic.fna files."""
    print( "\n--- Choose genomic/cds FNA files ---" )
    print( "💡 Files and directories autocompletion with TAB is available." )
    print( "Possible way to include data files:" )
    print( "  1. Path to directory and file pattern ('./Mesophiles/*.cds.fna')" )
    print( "  2. Point single file ('./Mesophiles/Escherichia_coli.cds.fna')" )
    print( "  3. Press ENTER (empty line), once ready with data files." )
    print( "-" * 65 )

    selected_files = []

    while True:
        user_input = input(
            f"To add more files type file/pattern or press ENTER to continue (currently added: {len(selected_files)}): "
        ).strip()

        if not user_input:
            if selected_files:
                break
            else:
                print(
                    "[!] At least one file or pattern must be added."
                )
                continue

        matched_files = glob.glob( os.path.expanduser( user_input ) )

        if not matched_files:
            print(f"[!] No file matched the pattern: '{user_input}'")
            continue

        added_count = 0
        for f in matched_files:
            if os.path.isfile( f ) and f not in selected_files:
                selected_files.append( f )
                added_count += 1

        print( f" -> Added {added_count} new files." )

    return selected_files


# ==============================================================================
# OPTION 1: GC FROM GENOME
# ==============================================================================
def run_option_1():
    ensure_results_dir()
    print( "\n" + "=" * 65 )
    print( "OPTION 1: GC ANALYSIS FROM WHOLE DNA GENOME (GC from genome)" )
    print( "=" * 65 )

    genome_files = get_genome_files()
    print( f"\n[DONE] {len(genome_files)} genomic files choosen." )

    print( "\n--- Configure sorting ---" )
    sort_flag = input( "Sorting flag (sort/nosort) [nosort]: " ).strip() or "nosort"

    print( "\n--- Calculations results file (under dir 'results/') ---" )
    raw_store_file = (
        input( "Input result_name.txt file [gc_results.txt]: ").strip()
        or "gc_results.txt"
    )
    store_file = RESULTS_DIR / Path( raw_store_file ).name

    print( "\n--- Output graph (in dir 'results/') ---" )
    raw_output_graph = (
        input( "Input graph_name.png file [gc_scatter.png]: " ).strip()
        or "gc_scatter.png"
    )
    output_graph = RESULTS_DIR / Path( raw_output_graph ).name

    print( "\n" + "-" * 65 )
    print( "STARTING CALCULATIONS..." )
    print( "-" * 65 )

    # 1. Calculations
    print( f"\n[1/2] Calculate GC content ({CALC_SCRIPT_1.name})..." )
    cmd_calc = [
        sys.executable,
        str( CALC_SCRIPT_1 ),
        sort_flag,
        str( store_file ),
    ] + genome_files
    if subprocess.run( cmd_calc ).returncode != 0:
        print( "\n[ERROR] Calculation failed.", file=sys.stderr )
        return

    # 2. Plot
    print( f"\n[2/2] Generate the graph ({PLOT_SCRIPT_1.name})..." )
    cmd_plot = [
        sys.executable,
        str( PLOT_SCRIPT_1 ),
        str( store_file ),
        str( output_graph ),
    ]
    if subprocess.run(cmd_plot).returncode != 0:
        print( "\n[ERROR] Graph generation failed.", file=sys.stderr )
        return

    print( f"\n[DONE] Task is ready." )
    print( f"  -> Calculation file:  '{store_file}'" )
    print( f"  -> Graph file:        '{output_graph}'" )


# ==============================================================================
# OPTION 2: GC / GC3 FROM CDS
# ==============================================================================
def run_option_2():
    ensure_results_dir()
    print( "\n" + "=" * 65 )
    print( "OPTION 2: GC AND GC3 ANALYSIS FROM CDS (GC and GC3 from CDS)" )
    print( "=" * 65 )

    genome_files = get_genome_files()
    print( f"\n[DONE] Selected {len(genome_files)} CDS files." )

    print( "\n--- Configure sorting ---" )
    sort_flag = input( "Sorting flag (sort/nosort) [nosort]: " ).strip() or "nosort"

    print( "\n--- Calculation results file (under dir 'results/') ---" )
    raw_store_file = (
        input( "Input result_name.txt file [gc_gc3_results.txt]: ").strip()
        or "gc_gc3_results.txt"
    )
    store_file = RESULTS_DIR / Path( raw_store_file ).name

    # Step 1: Calculations
    print( "\n" + "-" * 65 )
    print( f"[1/2] Calculate GC and GC3 ({CALC_SCRIPT_2.name})..." )
    print( "-" * 65 )

    cmd_calc = [
        sys.executable,
        str( CALC_SCRIPT_2 ),
        sort_flag,
        str( store_file ),
    ] + genome_files
    if subprocess.run(cmd_calc).returncode != 0:
        print( "\n[ERROR] Calculation failed.", file=sys.stderr )
        return

    if not store_file.exists():
        print( f"\n[ERROR] The file '{store_file}' missing.", file=sys.stderr )
        return

    # Step 2: Choose which graph to generate
    print( "\n--- Choose generation of graph ---" )
    print( "1. Scatter plot just GC3" )
    print( "2. Scatter plot both GC and C3 (on the same place)" )
    print( "3. Scatter plot GC in relation to GC3 (GC to GC3)" )
    print( "4. Generate all 3 graphs" )
    print( "-" * 65)

    plot_choice = input(" Choose option [1-4]: ").strip()

    chosen_plots = []
    if plot_choice in PLOT_SCRIPTS_2:
        chosen_plots.append( PLOT_SCRIPTS_2[ plot_choice ] )
    elif plot_choice == "4":
        chosen_plots = list( PLOT_SCRIPTS_2.values() )
    else:
        print( "\n[!] Invalid choice. Omit graph generation." )
        return

    print( "\n" + "-" * 65 )
    print( "[2/2] Generate choosen graphs under dir 'results/'..." )
    print( "-" * 65 )

    for desc, script_path, default_img in chosen_plots:
        raw_img = (
            input(
                f"Filename for '{desc}' [{default_img}]: "
            ).strip()
            or default_img
        )
        custom_img = RESULTS_DIR / Path( raw_img ).name

        print(f" -> Execute {script_path.name}...")
        cmd_plot = [
            sys.executable,
            str( script_path ),
            str( store_file ),
            str( custom_img ),
        ]

        if subprocess.run( cmd_plot ).returncode != 0:
            print(
                f"[ERROR] Could not generate {custom_img}",
                file=sys.stderr,
            )
        else:
            print( f"    └─ Generated graph: '{custom_img}'" )

    print( f"\n[DONE] Analysis by Option 2 successfully finished." )


# ==============================================================================
# MAIN MENU
# ==============================================================================
def main_menu():
    check_scripts()
    setup_readline()

    while True:
        print_header()
        print( "1. Generate GC content from whole genome (GC from genome)" )
        print( "2. Generate GC & GC3 content from CDS (GC / GC3 from CDS)" )
        print( "0. Exit" )
        print( "-" * 65 )

        choice = input( "Choose option [1/2/0]: " ).strip()

        if choice == "1":
            run_option_1()
            input( "\nPress ENTER to return to main menu...")
        elif choice == "2":
            run_option_2()
            input( "\nPress ENTER to return to main menu...")
        elif choice == "0":
            print( "\nGoodbye..." )
            sys.exit(0)
        else:
            print( "\n[!] Invalid option. Try again ]1/2/0].\n" )


if __name__ == "__main__":
    main_menu()

