#include <iostream>
#include <fstream>
#include <vector>
#include <numeric>

#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix.hpp>

using namespace std;

vector<vector<double>> read_counts(const string path)
{
    vector<vector<double>> mat;

    boost::numeric::ublas::matrix<double> m;
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

template <class T>
double avg(vector<T> v)
{
    float average = accumulate( v.begin(), v.end(), 0.0)/v.size();
    return average;
}

int main() {

    vector<vector<double>> mat;

    mat = read_counts("/Users/mtuncel/git_repos/sc-dna/input_data/norm_counts.tsv");
    // disp_vec(mat);

    // the last breakpoint
    u_int bp_size = mat[0].size();
    u_int window_size = 3;
    cout << bp_size << endl;

    // TODO:  iterate over cell_no values
    u_int cell_no = 0;

    for (int i = 0; i < bp_size; ++i) {

        vector<double> vect = mat[cell_no];

        int start = i - window_size;
        int end = i + window_size;

        // end cannot be equal to bp_size because cellranger DNA is putting 0 to the end
        if (start < 0 || end >= bp_size)
        {
            cout << "continue" <<endl;
            continue;
        }

        vector<double> lbins = vector<double>(vect.begin() + start, vect.begin() + i);
        vector<double> rbins = vector<double>(vect.begin() + i, vect.begin() + end);
        vector<double> all_bins = vector<double>(vect.begin() + start, vect.begin() + end);
        cout << lbins.size() <<' ' << rbins.size() << ' ' << all_bins.size() << ' ' << i;
        cout << endl;
        cout << avg(lbins) << ' ' << avg(rbins) << ' ' << avg(all_bins) <<endl;

        s

    }

    return 0;
}