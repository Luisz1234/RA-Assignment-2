# Balanced Allocation Strategies Simulation

This repository contains the C++ source code for Assignment 2: Balanced Allocation. The program empirically studies various "balls and bins" allocation strategies, analyzes their performance by measuring the load gap, and generates plots to visualize the results.

## Program Overview

The C++ program performs a comprehensive set of simulations to compare different load balancing schemes. Its core functions include:

1.  **Strategy Simulation:** Implements and simulates a wide range of allocation strategies:
    *   **Standard d-Choice:** One-choice, two-choice, and three-choice strategies.
    *   **Hybrid (1+β)-Choice:** A probabilistic mix between one-choice and two-choice strategies.
    *   **B-Batched Setting:** Simulates the effect of outdated information by processing allocations in batches of size `b`.
    *   **Partial Information:** Models scenarios where decisions are made with limited binary information about bin loads (e.g., "is the load above the median?").

2.  **Metric Calculation:** For each strategy, the program tracks the evolution of key performance metrics as the number of balls (`n`) increases up to the heavy-load scenario of `n = m^2`:
    *   **Average Gap:** The primary metric, `G_n = max(load) - average(load)`, averaged over `T` independent runs.
    *   **Standard Deviation:** The standard deviation of the gap across the `T` runs to measure strategy consistency and predictability.

3.  **Automated Plotting:** Directly interfaces with Gnuplot to automatically generate `.png` figures that compare the performance of the various strategies.

## Prerequisites

To compile and run the program, the following software must be installed and accessible in your system's PATH:

1.  **C++ Compiler:** A C++11 (or newer) compatible compiler (e.g., `g++` or Clang).
2.  **Gnuplot:** The plotting utility, used for generating all figures.

## Instructions

1.  **Compilation:**
    Save the code as `balanced_allocation.cpp`. Navigate to the directory containing the file and compile it using a C++ compiler.

    ```bash
    g++ balancedalloc.cpp
    ```

2.  **Execution:**
    Run the compiled executable from the command line. The program is non-interactive and has the simulation parameters (`m=100`, `T=20`) hardcoded.

    ```bash
    ./a.out
    ```

    Upon execution, the program will first run all simulations, printing progress to the console. It will then automatically generate the data and plot files in the same directory. The entire process may take a few minutes to complete.

## Generated Output

Running the program produces the following files in the current directory:

*   **Data Files (`.dat`):** A set of text files containing the raw data for each experiment, used by Gnuplot for plotting. Examples include:
    *   `one_choice_stats.dat`
    *   `two_choice_stats.dat`
    *   `three_choice_stats.dat`
    *   `beta_0.2_stats.dat`
    *   `batched_b10.dat`
    *   `batched_b100.dat`
    *   `batched_b200.dat`
    *   `batched_b1000.dat`
    *   `partial_k1_stats.dat`
    *   `partial_k2_stats.dat`

*   **Plot Images (`.png`):** Image files containing the final plots for the report.
    *   `plot_standard_strategies.png`
    *   `plot_variance.png`
    *   `plot_batched_strategies.png`
    *   `plot_partial_info.png`
