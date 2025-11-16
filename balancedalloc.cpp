#include <iostream>
#include <vector>
#include <random>
#include <numeric>
#include <algorithm>
#include <fstream>
#include <string>
#include <functional>
#include <cmath>
#include <stdexcept>

using namespace std;

// --- Gnuplot Pipe Class ---
#ifdef _WIN32
    #define POPEN _popen
    #define PCLOSE _pclose
#else
    #define POPEN popen
    #define PCLOSE pclose
#endif
class Gnuplot {
private:
    FILE* pipe;
public:
    Gnuplot() {
        pipe = POPEN("gnuplot", "w"); 
        if (!pipe) { throw runtime_error("Gnuplot not found or failed to open pipe."); }
    }
    ~Gnuplot() {
        if (pipe) { PCLOSE(pipe); }
    }
    void command(const string& cmd) {
        fprintf(pipe, "%s\n", cmd.c_str());
        fflush(pipe);
    }
};

// --- Random Number Generation Setup ---
mt19937 create_random_engine() {
    random_device rd;
    return mt19937(rd());
}
auto rng = create_random_engine();

// --- Helper Functions ---
double get_median(const vector<int>& bins) {
    vector<int> temp = bins;
    size_t n = temp.size() / 2;
    nth_element(temp.begin(), temp.begin() + n, temp.end());
    if (temp.size() % 2 == 1) {
        return static_cast<double>(temp[n]);
    } else {
        nth_element(temp.begin(), temp.begin() + n - 1, temp.end());
        return 0.5 * (temp[n - 1] + temp[n]);
    }
}

int get_percentile_val(const vector<int>& bins, double percentile) {
    vector<int> temp = bins;
    size_t index = static_cast<size_t>(temp.size() * percentile);
    if (index >= temp.size()) index = temp.size() - 1;
    nth_element(temp.begin(), temp.begin() + index, temp.end());
    return temp[index];
}


// --- Allocation Strategies ---

void one_choice_strategy(vector<int>& bins) {
    uniform_int_distribution<int> dist(0, bins.size() - 1);
    bins[dist(rng)]++;
}

void two_choice_strategy(vector<int>& bins) {
    uniform_int_distribution<int> dist(0, bins.size() - 1);
    int idx1 = dist(rng);
    int idx2 = dist(rng);
    if (bins[idx1] < bins[idx2]) {
        bins[idx1]++;
    } else if (bins[idx2] < bins[idx1]) {
        bins[idx2]++;
    } else {
        uniform_int_distribution<int> coin(0, 1);
        if (coin(rng) == 0) bins[idx1]++; else bins[idx2]++;
    }
}

// Picks three bins and places the ball in the one with the minimum load.
void three_choice_strategy(vector<int>& bins) {
    uniform_int_distribution<int> dist(0, bins.size() - 1);
    int idx1 = dist(rng);
    int idx2 = dist(rng);
    int idx3 = dist(rng);

    int min_load = bins[idx1];
    int chosen_idx = idx1;

    if (bins[idx2] < min_load) {
        min_load = bins[idx2];
        chosen_idx = idx2;
    }
    if (bins[idx3] < min_load) {
        chosen_idx = idx3;
    }
    // Note: A simple tie-break rule is to favor the first index found (idx1, then idx2, then idx3).
    bins[chosen_idx]++;
}

void beta_choice_strategy(vector<int>& bins, double beta) {
    uniform_real_distribution<double> dist(0.0, 1.0);
    if (dist(rng) < beta) two_choice_strategy(bins); else one_choice_strategy(bins);
}

void partial_info_strategy(vector<int>& bins, int k) {
    uniform_int_distribution<int> dist(0, bins.size() - 1);
    int idx1 = dist(rng);
    int idx2 = dist(rng);
    if (idx1 == idx2) idx2 = (idx1 + 1) % bins.size();

    // --- First Question: "Is the load above the median?" ---
    double median_load = get_median(bins);
    bool above_median1 = bins[idx1] > median_load;
    bool above_median2 = bins[idx2] > median_load;

    // Case 1: Answers differ
    if (above_median1 != above_median2) {
        if (!above_median1) bins[idx1]++; else bins[idx2]++;
        return;
    }

    // Case 2: We can only ask one question (k=1), so we must break the tie randomly.
    if (k == 1) {
        uniform_int_distribution<int> coin(0, 1);
        if (coin(rng) == 0) bins[idx1]++; else bins[idx2]++;
        return;
    }

    // Case 3: We can ask a second question (k=2).
    if (k == 2) {
        auto random_tie_break = [&]() {
            uniform_int_distribution<int> coin(0, 1);
            if (coin(rng) == 0) bins[idx1]++; else bins[idx2]++;
        };

        // Scenario A: Both bins are ABOVE the median load.
        if (above_median1 && above_median2) {
            // Second Question: "Is this bin among the 25% most loaded?" (i.e., > 75th percentile)
            int p75_load = get_percentile_val(bins, 0.75);
            bool top_25_percent1 = bins[idx1] > p75_load;
            bool top_25_percent2 = bins[idx2] > p75_load;

            if (top_25_percent1 != top_25_percent2) {
                // Answers differ, pick the one that is NOT in the top 25%.
                if (!top_25_percent1) bins[idx1]++; else bins[idx2]++;
            } else {
                // Answers are the same, so we pick at random.
                random_tie_break();
            }
        } 
        // Scenario B: Both bins are AT OR BELOW the median load.
        else {
            // Second Question: "Are they among the 75% most loaded?" (i.e., > 25th percentile)
            int p25_load = get_percentile_val(bins, 0.25);
            bool top_75_percent1 = bins[idx1] > p25_load;
            bool top_75_percent2 = bins[idx2] > p25_load;

            if (top_75_percent1 != top_75_percent2) {
                // Answers differ, pick the one that is NOT in the top 75%.
                if (!top_75_percent1) bins[idx1]++; else bins[idx2]++;
            } else {
                // Answers are the same, so we pick at random.
                random_tie_break();
            }
        }
    }
}


// --- Experiment Runners ---
void run_full_stat_experiment(int m, int T, function<void(vector<int>&)> strategy, const string& filename) {
    const int n_max = m * m;
    vector<double> total_gaps(n_max + 1, 0.0);
    vector<double> total_gaps_squared(n_max + 1, 0.0); // For variance calculation
    
    cout << "Running full-stat experiment for " << filename << "..." << endl;
    for (int t = 0; t < T; ++t) {
        vector<int> bins(m, 0);
        for (int n = 1; n <= n_max; ++n) {
            strategy(bins);
            int max_load = *max_element(bins.begin(), bins.end());
            double avg_load = static_cast<double>(n) / m;
            double gap = max_load - avg_load;
            
            total_gaps[n] += gap;
            total_gaps_squared[n] += gap * gap; 
        }
        cout << ".";
    }
    cout << "\nDone." << endl;

    ofstream outfile(filename);
    outfile << "# n\tavg_gap\tstd_dev\n";
    for (int n = 1; n <= n_max; ++n) {
        double avg_gap = total_gaps[n] / T;
        double avg_gap_sq = total_gaps_squared[n] / T;
        double variance = avg_gap_sq - (avg_gap * avg_gap);
        double std_dev = sqrt(max(0.0, variance));
        
        outfile << n << "\t" << avg_gap << "\t" << std_dev << "\n";
    }
}

void run_batched_experiment(int m, int T, int b, const string& filename) {
    const int n_max = m * m;
    vector<double> total_gaps(n_max + 1, 0.0);
    
    cout << "Running b-batched experiment (b=" << b << ") for " << filename << "..." << endl;
    uniform_int_distribution<int> dist(0, m - 1);

    for (int t = 0; t < T; ++t) {
        vector<int> bins(m, 0);
        int n = 0;
        while (n < n_max) {
            vector<int> bins_snapshot = bins;
            int current_batch_size = min(b, n_max - n);
            for (int i = 0; i < current_batch_size; ++i) {
                int idx1 = dist(rng);
                int idx2 = dist(rng);
                int chosen_idx = (bins_snapshot[idx1] <= bins_snapshot[idx2]) ? idx1 : idx2;
                bins[chosen_idx]++;
            }
            n += current_batch_size;
            int max_load = *max_element(bins.begin(), bins.end());
            double avg_load = static_cast<double>(n) / m;
            total_gaps[n] += (max_load - avg_load);
        }
        cout << ".";
    }
    cout << "\nDone." << endl;

    ofstream outfile(filename);
    outfile << "# n\tavg_gap\n";
    double last_gap = 0.0;
    for (int n = 1; n <= n_max; ++n) {
        if(total_gaps[n] != 0.0) {
            last_gap = total_gaps[n] / T;
        }
        if (last_gap > 0.0) {
            outfile << n << "\t" << last_gap << "\n";
        }
    }
}


// --- Main Program ---
int main() {
    const int m = 100; // Number of bins
    const int T = 20;  // Number of runs to average over

    // --- Part 1: Run All Simulations ---
    cout << "--- Starting Simulations ---" << endl;
    run_full_stat_experiment(m, T, one_choice_strategy, "one_choice_stats.dat");
    run_full_stat_experiment(m, T, two_choice_strategy, "two_choice_stats.dat");
    run_full_stat_experiment(m, T, three_choice_strategy, "three_choice_stats.dat"); // New experiment
    
    run_full_stat_experiment(m, T, [](auto& b){ beta_choice_strategy(b, 0.2); }, "beta_0.2_stats.dat");
    run_full_stat_experiment(m, T, [](auto& b){ beta_choice_strategy(b, 0.5); }, "beta_0.5_stats.dat");
    run_full_stat_experiment(m, T, [](auto& b){ beta_choice_strategy(b, 0.8); }, "beta_0.8_stats.dat");
    
    run_batched_experiment(m, T, 10, "batched_b10.dat");
    run_batched_experiment(m, T, m, "batched_b100.dat");
    run_batched_experiment(m, T, 2 * m, "batched_b200.dat");
    run_batched_experiment(m, T, 10 * m, "batched_b1000.dat");

    run_full_stat_experiment(m, T, [](auto& b){ partial_info_strategy(b, 1); }, "partial_k1_stats.dat");
    run_full_stat_experiment(m, T, [](auto& b){ partial_info_strategy(b, 2); }, "partial_k2_stats.dat");
    cout << "--- All simulations complete. ---" << endl;

    // --- Part 2: Generate Plots Directly ---
    cout << "\n--- Generating Plots via Gnuplot ---" << endl;
    try {
        Gnuplot gp;
        gp.command("set terminal pngcairo size 1024,768 enhanced font 'Verdana,10'");

        // Plot 1: Standard Strategies
        gp.command("set output 'plot_standard_strategies.png'");
        gp.command("set title 'Standard Allocation Strategies (m=" + to_string(m) + ")'");
        gp.command("set xlabel 'Number of Balls (n)'");
        gp.command("set ylabel 'Average Gap G_n'");
        gp.command("set grid");
        gp.command("set logscale y");
        gp.command("set key top left");
        gp.command("set arrow 1 from " + to_string(m) + ",graph 0 to " + to_string(m) + ",graph 1 nohead lc 'gray' dt 2");
        gp.command("plot 'one_choice_stats.dat' using 1:2 with lines title 'One-Choice (d=1)', \\\n"
                   "     'two_choice_stats.dat' using 1:2 with lines title 'Two-Choice (d=2)', \\\n"
                   "     'three_choice_stats.dat' using 1:2 with lines title 'Three-Choice (d=3)', \\\n"
                   "     'beta_0.2_stats.dat' using 1:2 with lines title '(1+0.2)-Choice', \\\n"
                   "     'beta_0.5_stats.dat' using 1:2 with lines title '(1+0.5)-Choice', \\\n"
                   "     'beta_0.8_stats.dat' using 1:2 with lines title '(1+0.8)-Choice'");
        gp.command("unset arrow 1");
        cout << "Generated plot_standard_strategies.png" << endl;
        

        // Plot 2: Batched Strategies 
        gp.command("set output 'plot_batched_strategies.png'");
        gp.command("set title 'Two-Choice in Batched Setting (m=" + to_string(m) + ")'");
        gp.command("set arrow 1 from " + to_string(m) + ",graph 0 to " + to_string(m) + ",graph 1 nohead lc 'gray' dt 2");
        gp.command("plot 'two_choice_stats.dat' using 1:2 with lines title 'Standard Two-Choice (b=1)', \\\n"
                   "     'batched_b10.dat' using 1:2 with lines title 'Batched (b=10)', \\\n"
                   "     'batched_b100.dat' using 1:2 with lines title 'Batched (b=m=" + to_string(m) + ")', \\\n"
                   "     'batched_b200.dat' using 1:2 with lines title 'Batched (b=2m=" + to_string(2*m) + ")', \\\n"
                   "     'batched_b1000.dat' using 1:2 with lines title 'Batched (b=10m=" + to_string(10*m) + ")'");
        gp.command("unset arrow 1");
        cout << "Generated plot_batched_strategies.png" << endl;

        // Plot 3: Partial Information 
        gp.command("set output 'plot_partial_info.png'");
        gp.command("set title 'Partial Information Strategies (m=" + to_string(m) + ")'");
        gp.command("set arrow 1 from " + to_string(m) + ",graph 0 to " + to_string(m) + ",graph 1 nohead lc 'gray' dt 2");
        gp.command("plot 'one_choice_stats.dat' using 1:2 with lines dashtype 2 title 'One-Choice (baseline)', \\\n"
                   "     'two_choice_stats.dat' using 1:2 with lines dashtype 3 title 'Two-Choice (perfect info)', \\\n"
                   "     'partial_k1_stats.dat' using 1:2 with lines title 'Partial Info (k=1 question)', \\\n"
                   "     'partial_k2_stats.dat' using 1:2 with lines title 'Partial Info (k=2 questions)'");
        gp.command("unset arrow 1");
        cout << "Generated plot_partial_info.png" << endl;

        // Plot 4: Visualizing Standard Deviation 
        gp.command("set output 'plot_variance.png'");
        gp.command("set title 'Average Gap with Standard Deviation (m=" + to_string(m) + ")'");
        gp.command("set xlabel 'Number of Balls (n)'"); // Added x-label for clarity
        gp.command("set ylabel 'Average Gap G_n'");
        gp.command("set yrange [0:30]"); // Set a fixed y-range for better comparison
        gp.command("set key top left");
        gp.command("unset logscale y"); // Variance is better viewed on a linear scale
        gp.command("set arrow 1 from " + to_string(m) + ",graph 0 to " + to_string(m) + ",graph 1 nohead lc 'gray' dt 2");
        gp.command("plot 'one_choice_stats.dat' u 1:($2-$3):($2+$3) with filledcurves fs transparent solid 0.2 lc 'purple' title 'One-Choice +/- 1 Std Dev', \\\n"
                   "     '' u 1:2 with lines lc 'purple' notitle, \\\n"
                   "     'two_choice_stats.dat' u 1:($2-$3):($2+$3) with filledcurves fs transparent solid 0.2 lc 'green' title 'Two-Choice +/- 1 Std Dev', \\\n"
                   "     '' u 1:2 with lines lc 'green' notitle");
        gp.command("unset arrow 1");
        cout << "Generated plot_variance.png" << endl;

    } catch (const exception& e) {
        cerr << "Plotting failed: " << e.what() << endl;
        cerr << "Please ensure 'gnuplot' is installed and in your system's PATH." << endl;
        return 1;
    }

    cout << "--- All plots generated successfully. ---" << endl;
    return 0;
}