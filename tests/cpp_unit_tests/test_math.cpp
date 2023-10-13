#include "openmc/math_functions.h"
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_vector.hpp>


TEST_CASE("Test evaluate_legendre")
{

    int max_order = 10;
    // for (int i = 0; i < test_coeffs.size(); i++) test_coeffs[i] = 0.5 * (2. * i + 1.);
    std::vector<double> test_x_vals = {-1.0, -0.5, 0.0, 0.5, 1.0};
    std::vector<double> ref_vals = {5.5, -0.455976486, -1.353515625, -2.773099899, 60.5};

    std::vector<double> test_coeffs(max_order + 1, 1.0);
    std::vector<double> test_vals(test_x_vals.size());

    for (int i = 0; i < ref_vals.size(); i ++) {
        test_vals[i] = openmc::evaluate_legendre(max_order, test_coeffs.data(), test_x_vals[i]);
    }

    CHECK_THAT(test_vals, Catch::Matchers::Approx(ref_vals).epsilon(1e-09));
}

TEST_CASE("Test calc_zn")
{
    int n = 10;
    double rho = 0.5;
    double phi = 0.5;

    // Reference solution
    std::vector<double> ref_vals =
    {1.00000000e+00, 2.39712769e-01, 4.38791281e-01,
     2.10367746e-01, -5.00000000e-01, 1.35075576e-01,
     1.24686873e-01, -2.99640962e-01, -5.48489101e-01,
     8.84215021e-03, 5.68310892e-02, -4.20735492e-01,
     -1.25000000e-01, -2.70151153e-01, -2.60091773e-02,
     1.87022545e-02, -3.42888902e-01, 1.49820481e-01,
     2.74244551e-01, -2.43159131e-02, -2.50357380e-02,
     2.20500013e-03, -1.98908812e-01, 4.07587508e-01,
     4.37500000e-01, 2.61708929e-01, 9.10321205e-02,
     -1.54686328e-02, -2.74049397e-03, -7.94845816e-02,
     4.75368705e-01, 7.11647284e-02, 1.30266162e-01,
     3.37106977e-02, 1.06401886e-01, -7.31606787e-03,
     -2.95625975e-03, -1.10250006e-02, 3.55194307e-01,
     -1.44627826e-01, -2.89062500e-01, -9.28644588e-02,
     -1.62557358e-01, 7.73431638e-02, -2.55329539e-03,
     -1.90923851e-03, 1.57578403e-02, 1.72995854e-01,
     -3.66267690e-01, -1.81657333e-01, -3.32521518e-01,
     -2.59738162e-02, -2.31580576e-01, 4.20673902e-02,
     -4.11710546e-04, -9.36449487e-04, 1.92156884e-02,
     2.82515641e-02, -3.90713738e-01, -1.69280296e-01,
     -8.98437500e-02, -1.08693628e-01, 1.78813094e-01,
     -1.98191857e-01, 1.65964201e-02, 2.77013853e-04};

    std::vector<double> test_vals((n+1) * (n+2) / 2.);

    openmc::calc_zn(n, rho, phi, test_vals.data());

    CHECK_THAT(test_vals, Catch::Matchers::Approx(ref_vals));
}