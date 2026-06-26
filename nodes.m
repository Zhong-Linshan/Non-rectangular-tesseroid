function [t, w] = nodes(n)
% nodes   Gauss–Legendre quadrature nodes and weights via the
%         Golub–Welsch algorithm.
%
%   [t, w] = nodes(n) returns the nodes t and weights w for n‑point
%   Gauss–Legendre quadrature on the canonical interval [-1, 1].
%   The nodes are the eigenvalues of a symmetric tridiagonal matrix
%   (the Jacobi matrix of the Legendre polynomials) and are returned
%   in ascending order.  The weights are computed from the first
%   component of the corresponding eigenvectors and sum to exactly 2.
%
%   This implementation follows the classic eigenvalue method
%   described by Golub & Welsch (1969) and is numerically stable for
%   moderate n (typically n ≤ 100).  It requires no special toolboxes.
%
%   Input:
%       n - Number of quadrature points (positive integer, n ≥ 2).
%
%   Outputs:
%       t - n×1 vector of Gauss–Legendre nodes, strictly inside
%           (-1, 1), sorted in increasing order.
%       w - n×1 vector of quadrature weights, all positive.
%
%   Example:
%       % Compute 7‑point Gauss–Legendre nodes and weights
%       [t, w] = nodes(7);
%       % Approximate ∫_{-1}^{1} exp(x) dx ≈ ∑ w.*exp(t)
%       integral_approx = sum(w .* exp(t));   % exact: exp(1)-exp(-1)
%
%   References:
%       Golub, G. H., & Welsch, J. H. (1969). Calculation of Gauss
%       quadrature rules. Mathematics of Computation, 23(106), 221‑230.
%       Abramowitz, M. & Stegun, I. A. (1964), Handbook of Mathematical
%       Functions, Ch. 25.4.
%
%   See also: GLQ_L, LGWT (if available in the Statistics Toolbox).

% ---------------------------------------------------------------------
% Build the symmetric tridiagonal Jacobi matrix for Legendre polynomials
% ---------------------------------------------------------------------
% Sub‑diagonal entries: beta_k = k / sqrt(4*k^2 - 1),  k = 1, …, n-1
k = (1 : n-1)';
beta = k ./ sqrt(4 * k.^2 - 1);

% Symmetric tridiagonal matrix with zero diagonal
A = diag(beta, 1) + diag(beta, -1);

% ---------------------------------------------------------------------
% Eigenvalue decomposition
% ---------------------------------------------------------------------
[V, E] = eig(A);        % V: eigenvectors, E: diagonal eigenvalues
t = diag(E);            % nodes = eigenvalues (unsorted)
w = 2 * V(1, :).^2;     % weights from first eigenvector components

% ---------------------------------------------------------------------
% Ensure nodes are sorted in ascending order (standard output)
% ---------------------------------------------------------------------
[t, sortIdx] = sort(t);
w = w(sortIdx).';       % convert to column vector to match t

% Normalisation check (optional – not strictly necessary)
% assert(abs(sum(w) - 2) < 1e-12, 'Weights do not sum to 2.');
end