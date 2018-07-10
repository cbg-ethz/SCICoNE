#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <unordered_map>
#include <string>

#include <iterator>


#include "MathOp.h"
#include "Tree.h"

using namespace std;

vector<vector<double>> read_counts(const string path)
{
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

    for (auto const &v1: vec) {
        for (auto const &v2: v1)
            cout << v2 << ' ';
        cout << endl;
    }
}

vector<vector<double>> likelihood_ratio(vector<vector<double>> mat, double window_size) {
    MathOp mo = MathOp();
    // the last breakpoint
    u_int bp_size = mat[0].size();

    //cout << bp_size << endl;
    u_int n_cells = mat.size();

    // u_int cell_no = 0;
    vector<vector<double>> aic_vec;
    for (int i = 0; i < bp_size; ++i) {

        int start = i - window_size;
        int end = i + window_size;

        // end cannot be equal to bp_size because cellranger DNA is putting 0 to the end
        if (start < 0 || end >= bp_size) {
            //cout << "continue" << endl;
            continue;
        }

        vector<double> aic_cell = vector<double>();

        for (int j = 0; j < n_cells; ++j) {


            vector<double> vect = mat[j];

            vector<double> lbins = vector<double>(vect.begin() + start, vect.begin() + i);
            vector<double> rbins = vector<double>(vect.begin() + i, vect.begin() + end);
            vector<double> all_bins = vector<double>(vect.begin() + start, vect.begin() + end);
            /*cout << lbins.size() <<' ' << rbins.size() << ' ' << all_bins.size() << ' ' << i;
            cout << endl;
            cout << avg(lbins) << ' ' << avg(rbins) << ' ' << avg(all_bins) <<endl;
            */
            // k is the degrees of freedom of the segment model
            u_int k_segment = 1;
            double ll_segment = mo.log_likelihood(all_bins);
            double aic_segment = 2 * k_segment - 2 * ll_segment;

            // degrees of freedom is 2, lambda_r and lambda_l
            u_int k_break = 2;
            double ll_break = mo.log_likelihood(lbins) + mo.log_likelihood(rbins);
            double aic_break = 2 * k_break - 2 * ll_break;

            //cout << ll_break << ' ' << ll_segment << ' ' << aic_break << ' ' << aic_segment << endl;
            double aic_p = aic_segment - aic_break;
            //cout << aic_p << endl;



            aic_cell.push_back(aic_p);
        }
        aic_vec.push_back(aic_cell);

    }
    return aic_vec;
}

long double log_add(long double val1, long double val2)
{
    /*
     * uses the following identity:
     * log(a + b) = log(a * (1 + b/a)) = log a + log(1 + b/a)
     * make sure a > b
     * to avoid undefined values
     */

    // the input are in the log space
    long double a = val1;
    long double b = val2;

    if (b>a)
    {
        long double temp = a;
        a = b;
        b = temp;
    }

    long double res = a + log(1+ exp(b-a));

//    if (isinf(res))
//        cout << "inf result detected";

    return res;
}



vector<double> combine_scores(vector<double> aic_vec)
{

    vector<double> row1(aic_vec.size(), 0.0);
    vector<double> row2;
    vector<double> res(1,0.0);

    // iterate over n_cells
    for (int i = 0; i < aic_vec.size(); ++i) {
        row2.clear();
        // inner computation of row2
        for (int j = 0; j < aic_vec.size()-i; ++j) {
            double value;
            if (j==0)
            {
                value = row1[j] + aic_vec[j+i];
            }
            else
            {
                value = log_add(row2[j-1] , row1[j] + aic_vec[j+i]);

            }

            if (isinf(value))
                cout << "inf value detected";

            row2.push_back(value);

        }
        //cout << " row1 size: " << row1.size() << " row2 size: " << row2.size() << " max j: " << aic_vec.size()-i<<endl;
        row1.clear();
        row1 = row2;
        double last = row2.back();


        res.push_back(last);

    }



    return res;
}


int main() {


    // counts per region per cell
    int D[5][5] = {{39,37,45,49,30},{31,28,34,46,11},{69,58,68,34,21},{72,30,31,46,21},{50,32,20,35,13}};

    // region sizes
    int r[5] = {4,2,3,5,2};

    // build tree
    u_int ploidy = 2; // diploid
    Tree t(ploidy); // root is created

    t.random_insert({{0, 1}, {1, 1}});
    t.insert_at(1,{{1, 1}, {2, 1}});
    t.insert_at(1,{{1, 1}});
    t.insert_at(2,{{0, -1}});
    t.insert_at(2,{{3, -1}});


    vector<double> avg_scores;

    int n = std::size(D);
    for (int i = 0; i < n; ++i) {
        t.compute_tree(D[i], r);
        avg_scores.push_back(t.average_score());
    }

    // log likelihood of data given tree
    double sum = std::accumulate(avg_scores.begin(), avg_scores.end(), 0.0);


    t.traverse_tree();
    t.destroy();


/*

    //TODO log normalize the log values before taking the exp values!
    // Afterwards make sure that the results before & after are the same

    // parse input
    vector<vector<double>> mat;
    mat = read_counts("/Users/mtuncel/git_repos/sc-dna/input_data/norm_counts.tsv");

    // compute the AIC scores
    u_int window_size = 3;
    vector<vector<double>> aic_vec = likelihood_ratio(mat,window_size);


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
        auto res = combine_scores(vec);
        for (auto const &v2: vec)
            cout << v2 << ' ';
        cout <<endl;

    }
*/
    return 0;
}