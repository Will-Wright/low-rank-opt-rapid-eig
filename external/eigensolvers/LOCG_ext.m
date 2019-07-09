function [lam, v, res, info] = LOCG_ext(A, n, u0, nx, A_norm, Udeflate, shift, tol, maxitn, opts)
%
%   Locally optimal extended conjugate gradient method
%
% Input
%   A             n-by-n Hermitian matrix (numerical or function representation)
%   n             Dimension of eigenvalue problem
%   u0            Initial approximate eigenvector
%   nx            Size of extended CG subspace
%                 (i.e., nx = 0 is CG, nx = 1 is LOCG, nx = 2 is extended LOCG)
%   nrmA          Estimate of 2-norm of A
%   Udeflate      n-by-k0 matrix with columns spanning approximate invariant subspace 
%                 corresponding to 1st to k0th largest eigenvalues.  Assume U'*U = I.
%   shift         Real shift for Udeflate
%   tol           Relative error tolerance.
%                 Terminate on ||A*u - lam*u|| / (||A|| + |lam|) <= tol
%   maxitn        Maximum number of iterations
%   opts          TODO
%  

if isnumeric(A)
  Afun = @(x) A*x;
else
  Afun = A;
end



end
