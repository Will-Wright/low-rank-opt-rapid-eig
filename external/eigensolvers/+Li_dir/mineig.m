function [V,D]=mineig(A,k)
%
% This function return k smallest eigenvalues of a Hermitian matrix A
%
% Copyright by Ren-Cang Li, 2/6/2013
%
%---------------------------------------------------------------
%
%  Input
%
%        A     n-by-n; Hermitian
%        k     integer: k<=n
%              number of smallest eigenvalues asked
%
%  Output
%
%        V     n-by-k, eigenvector matrix with associated eigenvalues in D
%        D     k-by-1, eigenvalues D(1) <= ... <= D(k)
%
%---------------------------------------------------------------

n=size(A,1);
if k>n
   disp('mineig: too many eigenvalues asked');
   return
end

%[V,D]=eig(A);

[V,D]=schur(A);

D=real(diag(D)); [D,idx]=sort(D,'ascend');
D=D(1:k); V=V(:,idx(1:k));

err=norm(V'*V-eye(k),1);
if err>=1000*n*eps
   disp(strcat('mineig: returned eigenvectors are not orthonormal; err=',num2str(err)));
   %[n norm(A'-A,1)]
end
