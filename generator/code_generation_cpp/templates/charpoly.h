#pragma once
#include <Eigen/Eigen>

// Note: These methods work inplace on the matrix. I.e. if you plan on using A afterwards, make a copy.

/* Computes minimal poly using Kyrlov's method. */
template<typename Derived>
void charpoly_krylov(Eigen::MatrixBase<Derived> &A, double *p) {
	int n = A.rows();
	Eigen::MatrixXd V(n, n);

	// TODO: Figure out what a good choice is here
	V.col(0) = Eigen::VectorXd::Ones(n);
	V(n - 1, 0) = 0;

	for (int i = 1; i < n; i++)
		V.col(i) = A * V.col(i - 1);

	Eigen::VectorXd c(n);
	c = V.partialPivLu().solve(A*V.col(n - 1));

	p[n] = 1;
	for (int i = 0; i < n; i++)
		p[i] = -c(i);
}

/* Computes characteristic poly using Danilevsky's method. */
template<typename Derived>
void charpoly_danilevsky(Eigen::MatrixBase<Derived> &A, double *p, double pivtol = 1e-12) {
	int n = A.rows();

	for (int i = n - 1; i > 0; i--) {
		double piv = A(i, i - 1);

		if (std::abs(piv) < pivtol) {
			int piv_ind = 0;
			double piv_new = std::abs(piv);
			for (int j = 0; j < i - 1; j++) {
				if (std::abs(A(i, j)) > piv_new) {
					piv_new = std::abs(A(i, j));
					piv_ind = j;
				}
			}
			// Perform permutation
			A.row(i - 1).swap(A.row(piv_ind));
			A.col(i - 1).swap(A.col(piv_ind));
			piv = A(i, i - 1);
		}

		Eigen::VectorXd v = A.row(i);
		A.row(i - 1) = v.transpose()*A;

		Eigen::VectorXd vinv = (-1.0)*v;
		vinv(i - 1) = 1;
		vinv /= piv;
		vinv(i - 1) -= 1;
		Eigen::VectorXd Acol = A.col(i - 1);
		for (int j = 0; j <= i; j++)
			A.row(j) = A.row(j) + Acol(j)*vinv.transpose();


		A.row(i) = Eigen::VectorXd::Zero(n);
		A(i, i - 1) = 1;
	}
	p[n] = 1;
	for (int i = 0; i < n; i++)
		p[i] = -A(0, n - i - 1);
}


/* Computes characteristic poly using Danilevsky's method with full pivoting */
template<typename Derived>
void charpoly_danilevsky_piv(Eigen::MatrixBase<Derived> &A, double *p) {
	int n = A.rows();

	for (int i = n - 1; i > 0; i--) {

		int piv_ind = i - 1;
		double piv = std::abs(A(i, i - 1));

		// Find largest pivot
		for (int j = 0; j < i - 1; j++) {
			if (std::abs(A(i, j)) > piv) {
				piv = std::abs(A(i, j));
				piv_ind = j;
			}
		}
		if (piv_ind != i - 1) {
			// Perform permutation
			A.row(i - 1).swap(A.row(piv_ind));
			A.col(i - 1).swap(A.col(piv_ind));
		}
		piv = A(i, i - 1);

		Eigen::VectorXd v = A.row(i);
		A.row(i - 1) = v.transpose()*A;

		Eigen::VectorXd vinv = (-1.0)*v;
		vinv(i - 1) = 1;
		vinv /= piv;
		vinv(i - 1) -= 1;
		Eigen::VectorXd Acol = A.col(i - 1);
		for (int j = 0; j <= i; j++)
			A.row(j) = A.row(j) + Acol(j)*vinv.transpose();


		A.row(i) = Eigen::VectorXd::Zero(n);
		A(i, i - 1) = 1;
	}
	p[n] = 1;
	for (int i = 0; i < n; i++)
		p[i] = -A(0, n - i - 1);
}


/* Computes characteristic poly using Danilevsky's method.
   Also returns the transform T for the eigenvectors */
template<typename Derived>
void charpoly_danilevsky_T(Eigen::MatrixBase<Derived> &A, double *p, Eigen::MatrixBase<Derived> &T, double pivtol = 1e-12) {
	int n = A.rows();
	T.setIdentity();

	for (int i = n - 1; i > 0; i--) {
		double piv = A(i, i - 1);

		if (std::abs(piv) < pivtol) {
			int piv_ind = 0;
			double piv_new = std::abs(piv);
			for (int j = 0; j < i - 1; j++) {
				if (std::abs(A(i, j)) > piv_new) {
					piv_new = std::abs(A(i, j));
					piv_ind = j;
				}
			}
			// Perform permutation
			A.row(i - 1).swap(A.row(piv_ind));
			A.col(i - 1).swap(A.col(piv_ind));
			T.col(i - 1).swap(T.col(piv_ind));

			piv = A(i, i - 1);
		}

		Eigen::VectorXd v = A.row(i);
		A.row(i - 1) = v.transpose()*A;

		Eigen::VectorXd vinv = (-1.0)*v;
		vinv(i - 1) = 1;
		vinv /= piv;
		vinv(i - 1) -= 1;
		Eigen::VectorXd Acol = A.col(i - 1);
		for (int j = 0; j <= i; j++)
			A.row(j) = A.row(j) + Acol(j)*vinv.transpose();

		Eigen::VectorXd Tcol = T.col(i - 1);
		for (int j = 0; j < n - 1; j++)
			T.row(j) = T.row(j) + Tcol(j)*vinv.transpose();

		A.row(i) = Eigen::VectorXd::Zero(n);
		A(i, i - 1) = 1;
	}
	p[n] = 1;
	for (int i = 0; i < n; i++)
		p[i] = -A(0, n - i - 1);
}

/* Computes characteristic poly using Danilevsky's method with full pivoting.
   Also returns the transform T for the eigenvectors */
template<typename Derived>
void charpoly_danilevsky_piv_T(Eigen::MatrixBase<Derived> &A, double *p, Eigen::MatrixBase<Derived> &T) {	

	int n = A.rows();
	T.setIdentity();

	for (int i = n - 1; i > 0; i--) {

		int piv_ind = i - 1;
		double piv = std::abs(A(i, i - 1));

		// Find largest pivot
		for (int j = 0; j < i - 1; j++) {
			if (std::abs(A(i, j)) > piv) {
				piv = std::abs(A(i, j));
				piv_ind = j;
			}
		}
		if (piv_ind != i - 1) {
			// Perform permutation
			A.row(i - 1).swap(A.row(piv_ind));
			A.col(i - 1).swap(A.col(piv_ind));
			T.col(i - 1).swap(T.col(piv_ind));
		}
		piv = A(i, i - 1);

		Eigen::VectorXd v = A.row(i);
		A.row(i - 1) = v.transpose()*A;

		Eigen::VectorXd vinv = (-1.0)*v;
		vinv(i - 1) = 1;
		vinv /= piv;
		vinv(i - 1) -= 1;
		Eigen::VectorXd Acol = A.col(i - 1);
		for (int j = 0; j <= i; j++)
			A.row(j) = A.row(j) + Acol(j)*vinv.transpose();

		Eigen::VectorXd Tcol = T.col(i - 1);
		for (int j = 0; j < n - 1; j++)
			T.row(j) = T.row(j) + Tcol(j)*vinv.transpose();

		A.row(i) = Eigen::VectorXd::Zero(n);
		A(i, i - 1) = 1;
	}
	p[n] = 1;
	for (int i = 0; i < n; i++)
		p[i] = -A(0, n - i - 1);
}


/* Computes characteristic poly using La Budde's method.
	https://arxiv.org/abs/1104.3769
 */
template<typename Derived>
void charpoly_la_budde(Eigen::MatrixBase<Derived> &A, double *p) {
	int n = A.rows();

	// Compute hessenberg form
	Eigen::HessenbergDecomposition<Eigen::MatrixXd> hess(A);
	Eigen::MatrixXd H = hess.matrixH();

	Eigen::VectorXd beta = H.diagonal(-1);
	Eigen::MatrixXd c(n, n);
	c.setZero();

	// Precompute the beta products needed
	Eigen::MatrixXd beta_prod(n - 1, n - 1);
	beta_prod.col(0) = beta;
	for (int i = 0; i < n - 1; i++) {
		for (int j = 1; j < n - i - 1; j++) {
			beta_prod(i, j) = beta_prod(i, j - 1) * beta(i + j);
		}
	}

	c(0, 0) = -H(0, 0);
	c(0, 1) = c(0, 0) - H(1, 1);
	c(1, 1) = H(0, 0)*H(1, 1) - H(0, 1)*beta(0);

	for (int i = 2; i < n; i++) {
		c(0, i) = c(0, i - 1) - H(i, i);
		for (int j = 1; j <= i - 1; j++) {
			c(j, i) = c(j, i - 1) - H(i, i)*c(j - 1, i - 1);
			for (int m = 1; m <= j - 1; m++) {
				c(j, i) -= beta_prod(i - m, m - 1) * c(j - m - 1, i - m - 1) * H(i - m, i);
			}
			c(j, i) -= beta_prod(i - j, j - 1) * H(i - j, i);
		}
		c(i, i) = -H(i, i)*c(i - 1, i - 1);
		for (int m = 1; m <= i - 1; m++) {
			c(i, i) -= beta_prod(i - m, m - 1)*c(i - m - 1, i - m - 1)*H(i - m, i);
		}
		c(i, i) -= H(0, i)*beta_prod(0, i - 1);
	}

	p[n] = 1;
	for (int i = 0; i < n; i++)
		p[i] = c(n - i - 1, n - 1);
}











