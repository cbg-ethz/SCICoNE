#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <map>
#include <string>

#include <iterator>


#include "MathOp.h"
#include "Tree.h"
#include "Inference.h"


using namespace std;

vector<vector<double>> read_counts(const string path)
{
    /*
     * Parses the input data into a double vector.
     * */
    vector<vector<double>> mat;

    std::ifstream filein(path);

    int i = 0, j=0;
    for (std::string line; std::getline(filein, line); )
    {

        // push an empty vector
        mat.push_back(vector<double>());

        std::istringstream fline(line);
        j = 0;
        for(;;) {
            double val;
            fline >> val;
            if (!fline) break;
            mat[i].push_back(val);
            j++;
        }
        i++;

    }
    return mat;
}

void disp_vec(vector<vector<long double>>& vec) {
/*
 * displays the long double vector
 * */

    for (auto const &v1: vec) {
        for (auto const &v2: v1)
            cout << v2 << ' ';
        cout << endl;
    }
}
void disp_vec(vector<vector<double>>& vec) {
/*
 * displays the double vector
 * */

    for (auto const &v1: vec) {
        for (auto const &v2: v1)
            cout << v2 << ' ';
        cout << endl;
    }
}


int main() {

    // set a seed number for reproducibility
    //SingletonRandomGenerator::get_generator(42);


    // counts per region per cell
    vector<vector<int>> D = {{39,37,45,49,30},{31,28,34,46,11},{69,58,68,34,21},{72,30,31,46,21},{50,32,20,35,13}};

    // region sizes
    vector<int> r = {4,2,3,5,2};

    // move probabilities
    vector<float> move_probs = {1.0f,1.0f,1.0f,1.0f, 1.0f, 1.0f, 1.0f};

    Inference mcmc(size(r));

    mcmc.initialize_worked_example();
    // mcmc.random_initialize();
    mcmc.compute_t_table(D,r);

    mcmc.infer_mcmc(D, r, move_probs, 50000);
    mcmc.write_best_tree();

    mcmc.destroy();



    // parse input
    vector<vector<double>> mat;
    mat = read_counts("/Users/mtuncel/git_repos/sc-dna/input_data/norm_counts.tsv");

    // compute the AIC scores
    u_int window_size = 1;
    vector<vector<double>> aic_vec = MathOp::likelihood_ratio(mat,window_size);

    vector<vector<double>> sigma;

    int n_breakpoints = aic_vec.size();
    int n_cells = aic_vec[0].size();
    cout <<"n_breakpoints: " << n_breakpoints << " n_cells: " << n_cells <<endl;

    for (auto &vec: aic_vec) // compute sigma matrix
    {
        auto res = MathOp::combine_scores(vec);
        sigma.push_back(res);
    }

    vector<double> log_priors;
    for (int j = 0; j < n_cells; ++j) {
        log_priors.push_back(MathOp::breakpoint_log_prior(j, n_cells,0.5));
    }


    vector<vector<long double>> log_posterior;

    for (int k = 0; k < n_breakpoints; ++k) {
        log_posterior.push_back(vector<long double>());
        for (int j = 0; j < n_cells; ++j) {
            long double val = log_priors[j] + sigma[k][j];
            log_posterior[k].push_back(val);
        }
    }

    cout<<"log posterior: \n";
    disp_vec(log_posterior);


    vector<vector<long double>> posterior;
    int k_star = 4;

    vector<long double> s_p;
    for (int l = 0; l < n_breakpoints; ++l)
    {
        posterior.push_back(vector<long double>());

        long double max_num = *max_element(log_posterior[l].begin(), log_posterior[l].begin()+k_star-1);
        long double max_denom = *max_element(log_posterior[l].begin(), log_posterior[l].end());

        for (int j = 0; j < k_star - 1; ++j) {
            long double val =exp(log_posterior[l][j] - max_num);
            posterior[l].push_back(val);
        }
        for (int k = k_star -1 ; k < log_posterior[l].size(); ++k) {
            long double val =exp(log_posterior[l][k] - max_denom);
            posterior[l].push_back(val);
        }


        long double sp_num = std::accumulate(posterior[l].begin(), posterior[l].begin()+k_star-1, 0.0);
        long double sp_denom = std::accumulate(posterior[l].begin(), posterior[l].end(), 0.0);
        long double fraction = sp_num / sp_denom;

        s_p.push_back(-1*log( 1- (1-fraction)));
    }

    std::ofstream output_file("./s_p.txt");
    for (const auto &e : s_p) output_file << e << "\n";
    return EXIT_SUCCESS;
}
