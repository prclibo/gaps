#include <iostream>
#include <vector>
#include <numeric>
#include <Eigen/Dense>
#include "mex.h"
#include "matrix.h"

#if $(use_sparse_template)
#include <Eigen/Sparse>
#endif
#if $(use_sturm_eigensolver)
#define MAX_DEG  $(num_basis)
#include "sturm.h"
#include "charpoly.h"
#endif
#if $(use_sturm_dani_eigensolver)
#define MAX_DEG  $(num_basis)
// #include "sturm.h"
#define DEG  $(num_basis)
#include "sturm_mart.h"
#include "charpoly.h"
#endif

using namespace Eigen;

#if $(use_reduced_eigenvector_solver)
void fast_eigenvector_solver(double * eigv, int neig, Eigen::Matrix<double,$(num_basis),$(num_basis)> &AM, MatrixXcd &sols);
#endif

void solver_$(solv_name)($(function_param_declaration))
{
    // Compute coefficients
    
$(code_compute_coefficients)

    // Setup elimination template
    static const int coeffs0_ind[] = { $(coeffs0_ind) };
    static const int coeffs1_ind[] = { $(coeffs1_ind) };
        

$(code_setup_template)


    // Setup action matrix
    Matrix<double,$(num_available), $(num_basis)> RR;
    RR << -C12.bottomRows($(num_reducible)), Matrix<double,$(num_basis),$(num_basis)>::Identity($(num_basis), $(num_basis));

    static const int AM_ind[] = { $(AM_ind) };
    Matrix<double, $(num_basis), $(num_basis)> AM;
    for (int i = 0; i < $(num_basis); i++) {
        AM.row(i) = RR.row(AM_ind[i]);
    }

    MatrixXcd sols($(num_vars), $(num_basis));
    sols.setZero();

    // Solve eigenvalue problem
#if $(use_standard_eigensolver)
    EigenSolver<Matrix<double, $(num_basis), $(num_basis)> > es(AM);
    ArrayXcd D = es.eigenvalues();    
    ArrayXXcd V = es.eigenvectors();

    $(code_normalize_eigenvectors)
    $(code_extract_solutions)
#endif
#if $(use_eigsonly_eigensolver)

    EigenSolver<MatrixXd> es(AM, false);
    ArrayXcd D = es.eigenvalues();

    int nroots = 0;
    double eigv[$(num_basis)];
    for (int i = 0; i < $(num_basis); i++) {
        if (std::abs(D(i).imag()) < 1e-6)
            eigv[nroots++] = D(i).real();
    }

    fast_eigenvector_solver(eigv, nroots, AM, sols);
#endif
#if $(use_sturm_eigensolver)
    double p[1+$(num_basis)];
    Matrix<double, $(num_basis), $(num_basis)> AMp = AM;
    charpoly_$(charpoly_method)(AMp, p);    
    double roots[$(num_basis)];
    int nroots;
    find_real_roots_sturm(p, $(num_basis), roots, &nroots, 8, 0);
    fast_eigenvector_solver(roots, nroots, AM, sols);
#endif

#if $(use_sturm_dani_eigensolver)

    double p[1 + $(num_basis)];
    Matrix<double, $(num_basis), $(num_basis)> T;
    charpoly_danilevsky_piv_T(AM, p, T);
    double roots[$(num_basis)];
    int nroots;
    // find_real_roots_sturm(p, $(num_basis), roots, &nroots, 8, 0);
    nroots = realRoots(p, roots);
    sols.resize($(num_vars), nroots);

    Eigen::MatrixXd V($(num_basis), nroots);
    Eigen::Map<Eigen::MatrixXd> D(roots, 1, nroots);
    V.bottomRows(1).setConstant(1);
    for (int j = $(num_basis) - 2; j >= 0; j--) {
        V.row(j) = V.row(j + 1).array() * D.array();
    }
    V = T * V;
    Eigen::RowVectorXd row = V.row(0);
    V.array().rowwise() /= row.array();

    D.transposeInPlace();
    $(code_extract_solutions)

#endif
    $(code_pack_outputs)

}
$(debug_comments)

#if $(use_reduced_eigenvector_solver)
    void fast_eigenvector_solver(double * eigv, int neig, Eigen::Matrix<double,$(num_basis),$(num_basis)> &AM, MatrixXcd &sols) {
    static const int ind[] = { $(ind_non_trivial) };    
    // Truncated action matrix containing non-trivial rows
    Matrix<double, $(length_ind_non_trivial), $(num_basis)> AMs;
    double zi[$(max_power)];
    
    for (int i = 0; i < $(length_ind_non_trivial); i++)    {
        AMs.row(i) = AM.row(ind[i]);
    }
    for (int i = 0; i < neig; i++) {
        zi[0] = eigv[i];
        for (int j = 1; j < $(max_power); j++)
        {
            zi[j] = zi[j - 1] * eigv[i];
        }
        Matrix<double, $(AA_sz)> AA;
$(code_setup_reduced_eigenvalue_eq)

        Matrix<double, $(ind_unit), 1>  s = AA.leftCols($(ind_unit)).colPivHouseholderQr().solve(-AA.col($(ind_unit)));
$(code_extract_solutions)
    }
}
#endif

mxArray* convertToMatlabCell(std::vector<Eigen::MatrixXcd> const& sols) {
    mxArray* cell = mxCreateCellMatrix(1, sols.size());
    for (int i = 0; i < sols.size(); ++i) {
        Eigen::MatrixXcd sol = sols.at(i);
        mxArray* mat = mxCreateDoubleMatrix(sol.rows(), sol.cols(), mxCOMPLEX);
        memcpy(mxGetComplexDoubles(mat), sol.data(), sizeof(mxComplexDouble) * sol.size());
        // double* mat_r = mxGetPr(mat);
        // double* mat_i = mxGetPi(mat);
        // for (int r = 0; r < sol.rows(); ++r) {
        //     for (int c = 0; c < sol.cols(); ++c) {
        //         mat_r[c * sol.rows() + r] = sol(r, c).real();
        //         mat_i[c * sol.rows() + r] = sol(r, c).imag();
        //     }
        // }
        mxSetCell(cell, i, mat);
    }
    return cell;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    std::ios_base::sync_with_stdio(false);
    if (nrhs != $(input_size)) {
        mexErrMsgIdAndTxt("automatic_generator_cvpr:$(solv_name):nrhs", "One input required.");
    }
    if (nlhs != $(output_size)) {
        mexErrMsgIdAndTxt("automatic_generator_cvpr:$(solv_name):nlhs", "One output required.");
    }    
    if (!mxIsDouble(prhs[0]) || mxIsComplex(prhs[0])) {
        mexErrMsgIdAndTxt("automatic_generator_cvpr:$(solv_name):notDouble", "Input data must be type double.");
    }
    $(code_mex_function_body)
}

