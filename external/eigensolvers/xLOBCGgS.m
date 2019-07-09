function [lam, Y, res, hsty, info, iter_cache]=xLOBCGgS(A_in, n, B, X, BU, nrmAB, Nd, shift, tol, maxitn)
%
%  Expertly Locally Optimal Block Conjugate Gradient (LOBCG) method to compute 
%  the Nd smallest eigenvalues of A-lamda B on the orthogonal complement of span(U), 
%  i.e., the kth to (k+nb-1)st smallest eigenvalues,
%  where A and B are Hermitian, and B is positive definite, k=k0+1 and k0 is the 
%  number of columns in U, size(X)=[n,nb].
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
%              eigenvalues Assume U'*B*U=I.
%              Possibly k0=0, i.e., U is empty array.
%        nrmAB 2-by-1, nrmAB(1): estimate of ||A||
%                      nrmAB(2): estimate of ||B||
%        Nd    integer, # of eigenpairs needed to be computed by this call
%        shift real for shifting away known eigenvalues
%        tol   tolerance for testing convergence 
%
%  Output
%
%        lam  converged eigenvalue
%        Y     corresponding converged eigenvector in this call
%        res   residual error  ||A y - lam B y||
%        hsty  struct 
%              hsty.eig -- row i contains the eigenvalue approximations of ith iteration
%              hsty.res -- row i contains normalized residuals 
%                                ||A yi - lam_i B yi||/(||A||*||y_i||+|lam_i|*||B||*||y_i||)
%                          for ith iteration
%        info  struct
%              info.itn number of iterations
%              info.kc  number of converged eigenpairs by this call
%
% Note: The code will always run block nb LOBCG. Naturally if one eigenpair is determined converged,
%       we append a working block by the next approximation in the process.
% 
%---------------------------------------------------------------
%

% n=max(size(A)); 
% maxitn=min(round(0.2*n),2000);

if isnumeric(A_in)
  A = @(x) A_in*x;
else
  A = A_in;
end
if isempty(B)
  B = sparse(eye(n,n));
end
iter_cache.V = [];
iter_cache.lam = [];

info.num_Ax_mults = 0;

nb=size(X,2); 
k0=size(BU,2); 

kc=k0;  % kc = converged current
        % # of converged eigenpairs (including pass-ins)
Y=[];

X0=Li_dir.MGSg(B,X);
                
%[1 norm(X0'*B*X0-eye(nb),1)]

AX=A(X0);
info.num_Ax_mults = info.num_Ax_mults + size(X0, 2);
if kc > 0 
   AX=AX+shift*(BU*(BU'*X0));
end
Rho=X0'*AX; 
BX=B*X0;
R=AX-BX*Rho; 

hsty.eig=[]; hsty.res=[];
lam=[]; res=[];    % converged eigenvalues and corresponding residuals

Q=Li_dir.MGSg(B,R,X0);

% [X0,Q]'*A*[X0,Q]
tmp=Q'*AX; AQ=A(Q);
info.num_Ax_mults = info.num_Ax_mults + size(Q, 2);
if kc > 0 
   AQ=AQ+shift*(BU*(BU'*Q));
end
%   Rho = [X'*A*X   X'*A*Q;
%          Q'*A*X   Q'*A*Q];
Rho=[Rho      tmp';
     tmp     Q'*AQ];
[V,D]=Li_dir.mineig(Rho,2*nb);        % V  should have orthonormal columns
 
X1=[X0,Q]*V(:,1:nb);           % X1 should have B-orthonormal columns
AX1=[AX, AQ]*V(:,1:nb);        % shift is already included.
BX1=B*X1;
R=AX1-BX1*diag(D(1:nb));

nrmR=sqrt(sum(R.*conj(R)))./( nrmAB(1)+abs(D(1:nb)')*nrmAB(2) );
hsty.eig=[hsty.eig;       
           D(1:nb).']; 
hsty.res=[hsty.res;       
           nrmR];

% Convergence test
% Suppose that the smaller the eigenvalues, the faster the convergence
kk=1;
while kk<=nb & nrmR(kk)<=tol
   kk=kk+1;
end
kk=kk-1;
% # of converged eigenpairs is kk.  Theorectically, all nb pairs could coverge at once. 
% But that's unlikely, we assume they don't, i.e., kk<nb.
lam=[lam; D(1:kk)]; res=[res nrmR(1:kk)]; 
if kk>0
   Y=[Y, X1(:,1:kk)];
   BU=[BU, BX1(:,1:kk)]; kc=kc+kk; 
   
   % make working block to nb again
   X1a=[X0,Q]*V(:,nb+1:nb+kk);           % X1a should have B-orthonormal columns
   AXa=[AX, AQ]*V(:,nb+1:nb+kk);         % shift is already included.
   BXa=B*X1a;
   Ra=AXa-BXa*diag(D(nb+1:nb+kk));
   
   R=[R(:,kk+1:nb) Ra]; D=D(kk+1:kk+nb); X1=[X1(:,kk+1:nb) X1a];  
   AX=[AX1(:,kk+1:nb) AXa]; BX=[BX1(:,kk+1:nb) BXa];
   % for X0, we take what remains   
   X0=X0(:,kk+1:nb); 
else
   AX=AX1; BX=BX1;
end

err = norm(R(:, 1))/(nrmAB(1)+abs(Rho(1,1)*nrmAB(2)));

itn=1;

err_orth=[];
while err > tol && itn < maxitn && kc<k0+Nd,

%    XBX=X1'*BX; RX=chol(XBX); invRX=inv(RX);
%    X1=X1*invRX; AX=AX*invRX; 
    Rho0=X1'*AX;
    
    Q=Li_dir.MGSg(B,[R X0],X1);

    % [X1,Q]'*A*[X1,Q]
    tmp=Q'*AX; AQ=A(Q);
    info.num_Ax_mults = info.num_Ax_mults + size(Q, 2);
    if kc > 0 
       AQ=AQ+shift*(BU*(BU'*Q));
    end
    Rho=[Rho0 tmp';
          tmp Q'*AQ];
       
    [V,D]=Li_dir.mineig(Rho,2*nb);   % V  should have orthonormal columns
    X0=X1;
    X1=[X0,Q]*V(:,1:nb);      % X1 should have B-orthonormal columns
    AX1=[AX, AQ]*V(:,1:nb);   % shift is already included.
    BX1=B*X1;
    
    R=AX1-BX1*diag(D(1:nb));
    
    nrmR=sqrt(sum(R.*conj(R)))./( nrmAB(1)+abs(D(1:nb)')*nrmAB(2) );
    hsty.eig=[hsty.eig;       
               D(1:nb).']; 
    hsty.res=[hsty.res;       
               nrmR];
    
    % Convergence test
    % Suppose that the smaller the eigenvalues, the faster the convergence
    kk=1;
    while kk<=nb & nrmR(kk)<=tol
       kk=kk+1;
kk
    end
    kk=kk-1;
    % # of converged eigenpairs is kk.
    lam=[lam; D(1:kk)]; res=[res nrmR(1:kk)]; 
    if kk>0    
       Y=[Y, X1(:,1:kk)];
       BU=[BU, BX1(:,1:kk)]; kc=kc+kk; 
   
       % make working block to nb again
       X1a=[X0,Q]*V(:,nb+1:nb+kk);           % X1a should have B-orthonormal columns
       AXa=[AX, AQ]*V(:,nb+1:nb+kk);         % shift is already included.
       BXa=B*X1a;
       Ra=AXa-BXa*diag(D(nb+1:nb+kk));
   
       R=[R(:,kk+1:nb) Ra]; D=D(kk+1:kk+nb); X1=[X1(:,kk+1:nb) X1a];  
       AX=[AX1(:,kk+1:nb) AXa]; BX=[BX1(:,kk+1:nb) BXa];
       % for X0, we take what remains   
       X0=X0(:,kk+1:nb);  
    else
       AX=AX1; BX=BX1;       
    end

    if isempty(Y)
      iter_cache.V = [iter_cache.V, X1(:, 1)];
    else
      iter_cache.V = [iter_cache.V, Y(:, 1)];
    end
    if ~isempty(lam)
      iter_cache.lam = [iter_cache.lam; lam(1,1)];
    end
    err = norm(R(:, 1))/(nrmAB(1)+abs(Rho(1,1)*nrmAB(2)));
if kk > 0
fprintf('kk > 0, new eig locked\n');
end
if mod(itn,5) == 0
fprintf('iter = %4i, err %1.5e\n', itn, err);
end
    itn = itn+1;
end
info.itn = itn; info.kc=kc-k0; Y=[Y X1];
