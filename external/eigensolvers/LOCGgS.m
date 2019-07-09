function [lam, y, res, hsty, info, iter_cache]=LOCGgS(A_in, n, B, x, BU, nrmAB, shift, tol, maxitn)
%
%  Locally Optimal Conjugate Gradient (LOCG) method to compute the kth smallest eigenvalue of A-lamda B,
%  where A and B are Hermitian, and B is positive definite, k=k0+1 and k0 is the number 
%  of columns in U.
%
%  Deflation is done by shifting away known eigenvalues. Basically LOBCG on
%
%      A+shift*B*U*(B*U)'=A+shift*BU*BU', B
%
% RCL 7/21/2013
%---------------------------------------------------------------
%
%  Input
%
%        A     n-by-n Hermitian matrix
%        B     n-by-n Hermitian matrix and positive definite
%        x     n-vector, initial guess to the kth eigenvector
%              associated with (k0+1)st to (k0+k)th smallest eigenvalues
%        BU    n-by-k0, B*U, where U's columns span an approximate invariant subspace 
%              (accurate enough, however) associated with 1st to k0th smallest 
%              eigenvalues. Assume U'*B*U=I.
%              Possibly k0=0, i.e., U is empty array.
%        nrmAB 2-by-1, nrmAB(1): estimate of ||A||
%                      nrmAB(2): estimate of ||B||
%        shift real for shifting away known eigenvalues
%        tol   tolerance for testing convergence
%
%  Output
%
%        lam   computed eigenvalues
%        Y     corresponding eigenvectors
%        res   residual error  ||A y - lam B y||
%        hsty  *-by-2 
%              hsty(:,1) -- eigenvalue approximations
%              hsty(:,2) -- normalized residuals ||A yi - lam_i B*yi||/(||A||*||y_i||+|lam_i|*||B||*||y_i||)
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


%maxitn=min(round(0.1*n),400);
%n=max(size(A)); 



k0=size(BU,2);


Bx=B*x; xBx=real(x'*Bx); nrmx=sqrt(xBx);
x0=(1/nrmx)*x;  Bx=(1/nrmx)*Bx;

if k0==0,
   
   Ax=A(x0); rho=real(x0'*Ax);
   r0=Ax-rho*Bx; 
   nrmr0=norm(r0);
   err=nrmr0/(nrmAB(1)+abs(rho)*nrmAB(2)); 
   
   hsty=[rho err];
   
   itn=0;

%   Q=r0-x0*(Bx'*r0);  % now Q'*B*x0=0;
Q=-r0+x0*(Bx'*r0);
   BQ=B*Q; pBp=real(Q'*(BQ)); nrmp=sqrt(pBp);  Q=(1/nrmp)*Q;
   
%   disp('Here 1');
%   disp(norm([x0 Q]'*B*[x0 Q]-eye(2),1));
   AQ = A(Q); tmp = x0'*AQ;
   Rho1=[rho     tmp;
       conj(tmp) real(Q'*AQ)];
       
   [v,rho]=Li_dir.mineig(Rho1,1);
   x1=[x0,Q]*v; % x1'*B*x1 approx 1   
x1=-x1;
   y=Q; 
  
   Bx=B*x1; xBx=real(x1'*Bx); nrmx=sqrt(xBx);
   x1=(1/nrmx)*x1;  Bx=(1/nrmx)*Bx;
   
   Ax=[Ax,AQ]*(v/nrmx); % A*x1;
   %Ax=A*x1;
   
   r=Ax-rho*Bx; 
   
   while err > tol & itn < maxitn,
   
      Q=Li_dir.MGSg(B,[y, r], x1);
%      disp('Here 2')
%      disp(norm([x1 Q]'*B*[x1 Q]-eye(3),1));
      AQ = A(Q);
      
      tmp=x1'*AQ;
      Rho1=[rho     tmp;
           tmp'   Q'*AQ];
       
      [v,rho]=Li_dir.mineig(Rho1,1);
      x0=x1;
      x1=[x0,Q]*v; % x1'*B*x1 approx 1
      try
        y=Q*v(2:3); 
      catch ME
        ME
        itn
        break;
      end

      Bx=B*x1; xBx=real(x1'*Bx); nrmx=sqrt(xBx);
      x1=(1/nrmx)*x1;  Bx=(1/nrmx)*Bx;
   
      Ax=[Ax,AQ]*(v/nrmx); % A*x1;
      %Ax=A*x1;
   
      r=Ax-rho*Bx; 
      nrmr=norm(r);
      err=nrmr/(nrmAB(1)+abs(rho)*nrmAB(2));

      hsty=[hsty; rho err];
 
      itn = itn+1;
   end
   
else  % k0>0
   
   Ax=A(x0)+shift*(BU*(BU'*x0)); 
   rho=real(x0'*Ax);
   r0=Ax-rho*Bx; 
   nrmr0=norm(r0);
   err=nrmr0/(nrmAB(1)+abs(rho)*nrmAB(2)); 
   
   hsty=[rho err];
   
   itn=0;

   Q=r0-x0*(Bx'*r0);  % now Q'*B*x0=0;
   BQ=B*Q; pBp=real(Q'*(BQ)); nrmp=sqrt(pBp);  Q=(1/nrmp)*Q;
   
%   disp('Here 1');
%   disp(norm([x0 Q]'*B*[x0 Q]-eye(2),1));
   
   AQ=A(Q)+shift*(BU*(BU'*Q)); tmp=x0'*AQ;
   Rho1=[rho     tmp;
       conj(tmp) real(Q'*AQ)];
       
   [v,rho]=Li_dir.mineig(Rho1,1);
   x1=[x0,Q]*v; % x1'*B*x1 approx 1   
   y=Q; 
   
   Bx=B*x1; xBx=real(x1'*Bx); nrmx=sqrt(xBx);
   x1=(1/nrmx)*x1;  Bx=(1/nrmx)*Bx;
   
   Ax=[Ax,AQ]*(v/nrmx); % A*x1+shift*(BU*(BU'*x1));
   %Ax=A*x1+shift*(BU*(BU'*x1));
   
   r=Ax-rho*Bx; 
   
   while err > tol & itn < maxitn,
   
      Q=Li_dir.MGSg(B,[y, r], x1);
%      disp('Here 2')
%      disp(norm([x1 Q]'*B*[x1 Q]-eye(3),1));

      AQ = A(Q) + shift*(BU*(BU'*Q)); tmp=x1'*AQ;
      Rho1=[rho     tmp;
           tmp'   Q'*AQ];
       
      [v,rho]=Li_dir.mineig(Rho1,1);
      x0=x1;
      x1=[x0,Q]*v; % x1'*B*x1 approx 1   
      y=Q*v(2:3); 
   
      Bx=B*x1; xBx=real(x1'*Bx); nrmx=sqrt(xBx);
      x1=(1/nrmx)*x1;  Bx=(1/nrmx)*Bx;
   
      Ax=[Ax,AQ]*(v/nrmx); % A*x1+shift*(BU*(BU'*x1));
      %Ax=A*x1+shift*(BU*(BU'*x1));
   
      r=Ax-rho*Bx; 
      nrmr=norm(r);
      err=nrmr/(nrmAB(1)+abs(rho)*nrmAB(2));

      hsty=[hsty; rho err];
 
      itn = itn+1;
   end

end

info = itn; y=x1; lam=rho; res=err;




