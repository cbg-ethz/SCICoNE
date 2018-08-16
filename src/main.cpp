#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <unordered_map>
#include <string>

#include <iterator>


#include "MathOp.h"
#include "Tree.h"
#include "Inference.h"
#include "../tests/validation.h"

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

void disp_vec(vector<vector<double>> vec) {
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


    test_swap_label();
    test_weighted_sample();
    test_prune_reattach();
    test_weighted_prune_reattach();
    test_add_remove_event();

    // counts per region per cell
    vector<vector<int>> D = {{39,37,45,49,30},{31,28,34,46,11},{69,58,68,34,21},{72,30,31,46,21},{50,32,20,35,13}};

    // region sizes
    vector<int> r = {4,2,3,5,2};

    // move probabilities
    vector<float> move_probs = {1.0f,1.0f,1.0f,1.0f, 0.0f};

    Inference mcmc(size(r));

    mcmc.random_initialize();
    mcmc.compute_t_table(D,r);
    disp_vec(mcmc.t_scores);

    mcmc.infer_mcmc(D, r, move_probs);
    mcmc.write_best_tree();
    cout<< *mcmc.get_t();

    mcmc.destroy();



//
//    // build tree
//    u_int ploidy = 2; // diploid
//    Tree t(ploidy); // root is created
//
//
//    t.random_insert({{0, 1}, {1, 1}});
//    t.insert_at(1,{{1, 1}, {2, 1}});
//    t.insert_at(2,{{0, -1}});
//    t.insert_at(2,{{3, -1}});
//    t.insert_at(1,{{1, 1}});
//
//    vector<double> sum_scores;
//
//    int n = std::size(D);
//    for (int i = 0; i < n; ++i) {
//        t.compute_tree(D[i], r);
//        sum_scores.push_back(t.sum_score());
//    }
//
//    // log likelihood of data given tree
//    int m = 6; //number of nodes (compute it automatically TODO)
//    double ll_dgt = MathOp::vec_sum(sum_scores) ;
//    ll_dgt -= (n-1 + m) * log(n+1);
//
//    cout<<ll_dgt<<endl;
//
//    Tree t_prime(t, ploidy);
//
//    //t.prune(t.all_nodes[2]);
//
//    t.traverse_tree();
//
//    auto attached_node = t_prime.prune_reattach();
//
//
//    auto v1 = t.get_scores();
//    auto v2 = t_prime.get_scores();
//
//    t = t_prime;
//
//
//    if (attached_node != nullptr)
//        t_prime.compute_tree( D[4], r);
//    //t_prime.compute_stack(attached_node, D[4], n, r);
//
//
//    t_prime.traverse_tree();


//    t.destroy();
//    t_prime.destroy();






    //TODO log normalize the log values before taking the exp values!
    // Afterwards make sure that the results before & after are the same
/*
    // parse input
    vector<vector<double>> mat;
    mat = read_counts("/Users/mtuncel/git_repos/sc-dna/input_data/norm_counts.tsv");

    // compute the AIC scores
    u_int window_size = 3;
    vector<vector<double>> aic_vec = MathOp::likelihood_ratio(mat,window_size);


    // dynamic programming
//    for (auto const &v2: aic_vec[11])
//        cout << v2 << ' ';
//    cout << endl;

    //auto n_aic_vec = log_normalize(aic_vec[0]);

    cout <<aic_vec.size()<<endl;
    int i = 0;
    for (auto vec: aic_vec)
    {
        cout << i++ <<" --> ";
        auto res = MathOp::combine_scores(vec);
        for (auto const &v2: vec)
            cout << v2 << ' ';
        cout <<endl;

    }
*/
    return 0;
}
