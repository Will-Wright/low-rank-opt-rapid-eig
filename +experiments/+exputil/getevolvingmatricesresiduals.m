function [xDiffNorm, yDiffNorm, yRelErr, ENorm, ANorm, ARelErr] = getevolvingmatricesresiduals(results)
% EXPERIMENT.exputil.getevolvingmatricesresiduals
%
% Computes the norm of x_k, y_k and E_k = A_(k+1) - A_k in the sequence of 
% saga_sd iterates.

import util.*;

n = results.gendatapl.A.n;
maxIt = results.saga_sd.solInfo.nit;
xDiffNorm = zeros(maxIt, 1);
yDiffNorm = zeros(maxIt, 1);
yRelErr = zeros(maxIt, 1);
ENorm = zeros(maxIt, 1);
ANorm = zeros(maxIt, 1);
ARelErr = zeros(maxIt, 1);
x = results.saga_sd.solInfo.iterData.x{1};
y = results.saga_sd.solInfo.iterData.y{1};

k = 1;
sigma = 'LM';     % makes eigs find largest magnitude eigval for ||A||_2
eigsopts.tol = eps;
eigsopts.maxit = 1000;
eigsopts.p = min(n, 20);
eigsopts.issym = true;
eigsopts.isreal = false;

Afun  = @(x)util.vec(results.gendatapl.A.adjoint(y,x));
[~,dA] = eigs(Afun, n, k, sigma, eigsopts);
ANorm(1, 1) = abs(dA);

for i = 2:maxIt
   xPrev = x;
   x = results.saga_sd.solInfo.iterData.x{1, i};   
   xDiffNorm(i, 1) = norm(xPrev(:) - x(:));
   
   yPrev = y;
   y = results.saga_sd.solInfo.iterData.y{1, i};   
   yDiffNorm(i, 1) = norm(yPrev(:) - y(:));   
   yRelErr(i, 1) = yDiffNorm(i, 1) / norm(y(:));
   
   yDiff = y - yPrev;
   % Note: using 'vec' instead of 'rvec' to have lam1 match ||A||_2
   Efun  = @(x)util.vec(results.gendatapl.A.adjoint(yDiff,x));
   [~,dE] = eigs(Efun, n, k, sigma, eigsopts);
   % The two lines below verify lam1 matches ||A||_2 for full matrix A
   % E = experiments.exputil.getfulldualmatrix(results.gendatapl.A, yDiff);
   % ENorm(i-1, 1) = norm(E);
   ENorm(i-1, 1) = abs(dE);
   
   Afun  = @(x)util.vec(results.gendatapl.A.adjoint(y,x));
   [~,dA] = eigs(Afun, n, k, sigma, eigsopts);
   ANorm(i, 1) = abs(dA);
   ARelErr(i, 1) = ENorm(i-1, 1) / ANorm(i, 1);
end

   
end