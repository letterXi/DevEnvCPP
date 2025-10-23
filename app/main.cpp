#include "AmgclSolver.hpp"
#include "CsrMatrix.hpp"
#include <iostream>
#include <vector>

void Poisson(size_t n, std::vector<size_t>& addr, std::vector<size_t>& cols, std::vector<double>& vals,
             std::vector<double>& rhs) {
    size_t n2 = n * n;
    double h = 1.0 / (double)(n - 1);

    addr.clear();
    addr.reserve(n2 + 1);
    addr.push_back(0);

    cols.clear();
    cols.reserve(n2 * 5);

    vals.clear();
    vals.reserve(n2 * 5);

    rhs.resize(n2);

    for (size_t j = 0, k = 0; j < n; ++j) {
        for (size_t i = 0; i < n; ++i, ++k) {
            if (i == 0 || j == 0 || i == n - 1 || j == n - 1) {
                cols.push_back(k);
                vals.push_back(1);
                rhs[k] = 0.0;
            } else {
                cols.push_back(k - n);
                vals.push_back(-1.0 / (h * h));

                cols.push_back(k - 1);
                vals.push_back(-1.0 / (h * h));

                cols.push_back(k);
                vals.push_back(4.0 / (h * h));

                cols.push_back(k + 1);
                vals.push_back(-1.0 / (h * h));

                cols.push_back(k + n);
                vals.push_back(-1.0 / (h * h));

                rhs[k] = 1.0;
            }
            addr.push_back(cols.size());
        }
    }
}

int main() {

    std::vector<size_t> addr;
    std::vector<size_t> cols;
    std::vector<double> vals;
    std::vector<double> rhs;
    size_t n = 1;
    Poisson(n, addr, cols, vals, rhs);

    std::vector<double> x;

    CsrMatrix mat(std::move(addr), std::move(cols), std::move(vals));
    AmgclSolver solver({{"precond.relax.type", "gauss_seidel"}});
    solver.set_matrix(mat);
    solver.solve(rhs, x);
    return 0;
}