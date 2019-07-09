function [lam, Y, res, hsty, info, iter_cache]=LOBCGgS(A_in, n, B, X, BU, nrmAB, shift, tol, maxitn, opts)
%
%  Locally Optimal Block Conjugate Gradient (LOBCG) method to compute 
%  the (k0+1)st to (k0+nb)th smallest eigenvalues of A-lamda B, 
%  where A and B are Hermitian, and B is positive definite, k0 is the 
%  number of columns in BU=B*U and U (not passing in) contains k0 
%  approximate eigenvectors associated with the first k0 smallest eigenvalues of A-lamda B.
%
%  Deflation is done by shifting away known eigenvalues. Basically LOBCG on
%
%      A+shift*B*U*(B*U)'=A+shift*BU*BU', B
%
% Copyright by Ren-Cang Li, 6/5/2013
%---------------------------------------------------------------
%
%  Input
%
%        A     n-by-n Hermitian matrix
%        B     n-by-n Hermitian matrix and positive definite
%        X     n-by-nb whose columns span an approximate invariant subspace 
%              associated with (k0+1)st to (k0+nb)th smallest eigenvalues
%        BU    n-by-k0, B*U, where U's columns span an approximate invariant subspace 
%              (accurate enough, however) associated with 1st to k0th smallest 
%              eigenvalues. Assume U'*B*U=I.
%              Possibly k0=0, i.e., U is empty array.
%        nrmAB 2-by-1, nrmAB(1): estimate of ||A||
%                      nrmAB(2): estimate of ||B||
%        shift real for shifting away known eigenvalues
%
%  Output
%
%        lam   computed eigenvalues
%        Y     corresponding eigenvectors
%        res   residual error  ||A y - lam B y||
%        hsty  struct 
%              hsty.eig -- row i contains the eigenvalue approximations of ith iteration
%              hsty.res -- row i contains normalized residuals 
%                                ||A yi - lam_i B yi||/(||A||*||y_i||+|lam_i|*||B||*||y_i||)
%                          for ith iteration
%        info  struct
%              info.itn number of iterations
%
% 
%---------------------------------------------------------------
%


if isnumeric(A_in)
  A = @(x) A_in*x;
else
  A = A_in;
end
if isempty(B)
  B_is_identity = true;
  B = sparse(eye(n,n));
else
  B_is_identity = false;
end
if nargout > 5
  iter_cache.V = [];
  iter_cache.lam = [];
  do_iter_cache = true;
else
  do_iter_cache = false;
end

if isempty(opts)
  opts = [];
end
if ~isfield(opts, 'display') || opts.display ~=1
	opts.display = 0;
end
if opts.display
  PrintBanner();
end

info.num_mat_vec_mults = 0;
num_mat_vec_mults_prev = 0;

%n=max(size(A)); 
%maxitn=min(round(0.05*n),400);


nb=size(X,2); 
k0=size(BU,2); 
if isfield(opts, 'nb_desired')
  nb_desired = opts.nb_desired;
else
  nb_desired = nb;
end

which_imp=2;

% Commented out from original code
%X0=MGSg(B,X);

% If B = I, then X0 = X
%BX=B*X;
%XBX=X'*BX; cholXBX=chol(XBX); invcholXBX=inv(cholXBX);
%X0=X*invcholXBX;

%[1 norm(X0'*B*X0-eye(nb),1)]

%AX=A(X0);
AX = A(X);
info.num_mat_vec_mults = info.num_mat_vec_mults + size(X, 2);
if k0 > 0 
   AX=AX+shift*(BU*(BU'*X));
end
Rho=X'*AX; 
%BX=B*X0; 
%R=AX-BX*Rho;
R = AX - X*Rho;

hsty.eig=[]; hsty.res=[];

Q=Li_dir.MGSg(B,R,X);
%[2 norm(X0'*B*X0-eye(k0),1) norm(Q'*B*Q-eye(k0),1) norm([X0 Q]'*B*[X0 Q]-eye(k0*2),1) norm(X0'*AX-diag(D),1)]

% [X0,Q]'*A*[X0,Q]
tmp=Q'*AX; AQ=A(Q);
info.num_mat_vec_mults = info.num_mat_vec_mults + size(Q, 2);
if k0 > 0 
   AQ=AQ+shift*(BU*(BU'*Q));
end

Rho1=[Rho  tmp';
     tmp     Q'*AQ];
[V,D]=Li_dir.mineig(Rho1,nb);   % V  should have orthonormal columns
 
X1=[X,Q]*V;           % X1 should have B-orthonormal columns
AX=[AX, AQ]*V;         % shift is already included.
%BX=B*X1;
R=AX-X1*diag(D);
if which_imp==2
   Y=Q;
end

nrmR=sqrt(sum(R.*conj(R)))./( nrmAB(1)+abs(D')*nrmAB(2) );
hsty.eig=[hsty.eig;       
           D.']; 
hsty.res=[hsty.res;     
           nrmR];

% Only considers error for leading eigenvalue
%err = norm(R(:, 1))/(nrmAB(1)+abs(Rho(1,1)*nrmAB(2)));

for i = 1:nb
  err(i) = norm(R(:, i))/(nrmAB(1)+abs(D(i, 1)*nrmAB(2)));
end
err = max(nrmR(1:nb_desired));

if opts.display
  PrintIter(0, nrmR, info.num_mat_vec_mults - num_mat_vec_mults_prev);
end
num_mat_vec_mults_prev = info.num_mat_vec_mults;


itn=1;

while err > tol && itn < maxitn,
 
%    XBX=X1'*BX; cholXBX=chol(XBX); invcholXBX=inv(cholXBX);
%    X1=X1*invcholXBX; AX=AX*invcholXBX; 
    Rho=X1'*AX;
    
    if which_imp==1
       Q=Li_dir.MGSg(B,[R X0],X1);
    else
       Q=Li_dir.MGSg(B,[R Y],X1);
    end
    
    % [X1,Q]'*A*[X1,Q]
    tmp=Q'*AX; AQ=A(Q);
    info.num_mat_vec_mults = info.num_mat_vec_mults + size(Q, 2);
    if k0 > 0 
       AQ=AQ+shift*(BU*(BU'*Q));
    end
    Rho1=[Rho  tmp';
          tmp Q'*AQ];
       
    [V,D]=Li_dir.mineig(Rho1,nb);  % V  should have orthonormal columns
    X0=X1;
    X1=[X1,Q]*V;          % X1 should have B-orthonormal columns
    AX=[AX, AQ]*V;        % shift is already included.
    %
    % TODO: This line throws dimension errors when algorithm has converged
    %
    if which_imp==2
       Y=Q*V(nb+1:3*nb,:);
    end
    
%    BX=B*X1;
    R=AX-X1*diag(D);
    
    nrmR=sqrt(sum(R.*conj(R)))./( nrmAB(1)+abs(D')*nrmAB(2) );
    hsty.eig=[hsty.eig;       
               D.']; 
    hsty.res=[hsty.res;     
               nrmR];

    for i = 1:nb
      err(i) = norm(R(:, i))/(nrmAB(1)+abs(D(i, 1)*nrmAB(2)));
    end
		err = max(nrmR(1:nb_desired));

		if opts.display
       PrintIter(itn, nrmR, info.num_mat_vec_mults - num_mat_vec_mults_prev);
    end
    num_mat_vec_mults_prev = info.num_mat_vec_mults;

    itn = itn+1;
end
info.itn = itn; 
lam=D; res=nrmR; Y=X1;

end % end of function LOBCGgS



function [] = PrintBanner()
  fprintf('\n              LOBCG Prototype Solver\n');
  fprintf('             rel          rel       num  \n')
  fprintf(' Iter |     err 1        err 2      A*x  \n');
  fprintf('------|-------------------------------------------------\n');
end

function [] = PrintIter(itn, nrmR, num_mat_vec_mults)

    print_array = [itn, nrmR, num_mat_vec_mults];
    print_str = ' %3i  |  ';
    for i = 1:length(nrmR)
      print_str = strcat(print_str, '  %1.5e');
    end
    print_str = strcat(print_str, '  %4i \n');
    fprintf(print_str, print_array);
end


