function [lam, y, res, hsty, info, iter_cache]=CGgS(A_in, n, B, x, BU, nrmAB, shift, tol, maxitn, opts)
%
%  Conjugate Gradient method to compute the kth smallest eigenvalue of A-lamda B,
%  where A and B are Hermitian, and B is positive definite, k=k0+1 and k0 is the number 
%  of columns in U.
%
%  Deflation is done by shifting away known eigenvalues. Basically CG on
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
%        x     n-vector, initial guess to the kth eigenvector
%        BU    n-by-k0, B*U, where U's columns span an approximate invariant subspace 
%              (accurate enough, however) associated with 1st to k0th smallest 
%              eigenvalues. Assume U'*B*U=I.
%              Possibly k0=0, i.e., BU is empty array.
%        nrmAB 2-by-1, nrmAB(1): estimate of ||A||
%                      nrmAB(2): estimate of ||B||
%        shift real for shifting away known eigenvalues
%        tol   tolerance for testing convergence
%
%  Output
%
%        lam  converged eigenvalue
%        y     corresponding eigenvector
%        res   residual error  ||A y - lam B*y||
%        hsty  *-by-2 
%              hsty(:,1) -- eigenvalue approximations
%              hsty(:,2) -- normalized residuals ||A yi - lam_i B*yi||/(||A||*||y_i||+|lam_i|*||B||*||y_i||)
%        info  the number of CG steps
%
%----------------------------------------------------------------------
%
if isnumeric(A_in)
  A = @(x) A_in*x;
else
  A = A_in;
end
if isempty(B)
  B_is_identity = true;
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

info.num_mat_vec_mults = 0;
num_mat_vec_mults_prev = 0;

if ~isfield(opts, 'display') || opts.display ~=1
   opts.display = 0;
end
if opts.display
   PrintBanner();
end


%n=max(size(A));
%maxitn=round(0.1*n); 
whichbeta = 1; 
k0=size(BU,2); 
LSopts.Line=1;   % Line Search opts
howLS=0;         % =0 solve 2-by-2 eigenvalue problem: CAN FAIL if 'eig' returns x1 = -x1_true      
                 % =1 use Line Search function LNSRCHg(...)
if B_is_identity
  Bx = x;
else
  Bx=B*x; 
end
xBx=real(x'*Bx); nrmx=sqrt(xBx);
x0=(1/nrmx)*x;  Bx=(1/nrmx)*Bx;

if k0==0,
   
   Ax=A(x0); rho=x0'*Ax;
   info.num_mat_vec_mults = info.num_mat_vec_mults + size(x0, 2);
   r0=Ax-rho*Bx; 
   nrmr0=norm(r0);
   err=nrmr0/(nrmAB(1)+abs(rho)*nrmAB(2)); itn=0;

   p0=-r0; 
   
   hsty=[rho err];

   if opts.display
      PrintIter(itn, err, info.num_mat_vec_mults - num_mat_vec_mults_prev);
   end
   num_mat_vec_mults_prev = info.num_mat_vec_mults;

% NOTE: The following left-justified modifications in this while loop
%       were inserted to determine why LS significantly outperformed
%       eigmin, and how to fix this.
   while err > tol & itn < maxitn,
   
       if howLS==1,
          if B_is_identity
            Bp0 = p0;
          else
            Bp0=B*p0; 
          end
          pBp=real(p0'*(Bp0)); nrmp=sqrt(pBp);  p0n=(1/nrmp)*p0;
          Ap0n=A(p0n);
				  info.num_mat_vec_mults = info.num_mat_vec_mults + size(p0n, 2);
          rhos=[rho p0n'*(Ap0n)];
%q=p0-x0*(Bx'*p0);
%q = (1/norm(q))*q;
%p0n = q;
%Ap0n=A(p0n);
%rhos = [rho p0n'*(Ap0n)];

          LSopts.rhos=rhos;
          
          [x1, t, lam]=Li_dir.LNSRCHg(Ax,Bx,x0,p0n,LSopts);

%format long
%x1temp = x1;
%[t, lam]

       else
          q=p0-x0*(Bx'*p0);  % now q'*B*x0=0;
          if B_is_identity
            Bq = q;
          else
            Bq=B*q; 
          end
          pBp=real(q'*(Bq)); nrmp=sqrt(pBp);  q=(1/nrmp)*q;
   
          Aq=A(q); tmp=x0'*Aq;
					info.num_mat_vec_mults = info.num_mat_vec_mults + size(q, 2);
          Rho1=[rho     tmp;
              conj(tmp) real(q'*Aq)];
 
        % Same results for both eigmin methods
        % [v1,lam1]=Li_dir.mineig(Rho1,1);

          [v, lam] = eig(Rho1);
%lam
          [lam, lam_idx] = min(diag(lam));
          v = v(:, lam_idx);

          x1=[x0,q]*v; % x1'*B*x1 approx 1 
% TODO: Replace this, maybe with rigorous heuristic on signs of x0, q?
x1 = -x1;

% TESTING IF THERE'S A CORRECT SIGN FOR x1
%if sign(x1temp(1,1)) ~= sign(x1(1,1))
%  x1 = -x1;
%end

%Ax1temp = A(x1temp); Ax1 = A(x1);
%[x1temp(1:5, 1), x1(1:5, 1); ...
% norm(Ax1temp - x1temp'*Ax1temp*x1temp), norm(Ax1 - x1'*Ax1*x1); ...
% x1temp'*Ax1temp, x1'*Ax1]

       end

       if B_is_identity
         Bx = x1;
       else
         Bx=B*x1; 
       end
       
       xBx=real(x1'*Bx); nrmx=sqrt(xBx);
       x1=(1/nrmx)*x1;  Bx=(1/nrmx)*Bx;
       
       Ax=A(x1); rho=lam; %rho=x1'*Ax;
			 info.num_mat_vec_mults = info.num_mat_vec_mults + size(x1, 2);
       
       r1=Ax-rho*Bx;
       nrmr1=norm(r1);
       err=nrmr1/(nrmAB(1)+abs(rho)*nrmAB(2));

       hsty=[hsty; rho err];
       
       if whichbeta == 1,
          beta=nrmr1*nrmr1/(nrmr0*nrmr0); 
       else
          tmp=nrmx*r1;
          beta=r1'*(tmp-r0)/(nrmr0*nrmr0); 
       end
       p1=-r1+beta*p0;

			 if opts.display && mod(itn, 10) == 0
		 		  PrintIter(itn, err, info.num_mat_vec_mults - num_mat_vec_mults_prev);
		 	 end
			 num_mat_vec_mults_prev = info.num_mat_vec_mults;
       
       x0=x1; r0=r1; p0=p1; nrmr0=nrmr1;
       if do_iter_cache
         iter_cache.V = [iter_cache.V, x0];
         iter_cache.lam = [iter_cache.lam, lam];
       end
       itn = itn+1;
   end
   
else  % k0>0

   
   Ax=A(x0)+shift*(BU*(BU'*x0)); rho=x0'*Ax;
	 info.num_mat_vec_mults = info.num_mat_vec_mults + size(x0, 2);
   r0=Ax-rho*Bx; 
   nrmr0=norm(r0);
   err=nrmr0/(nrmAB(1)+abs(rho)*nrmAB(2)); itn=0;

   p0=-r0;
   
   hsty=[rho err];

   if opts.display
      PrintIter(itn, err, info.num_mat_vec_mults - num_mat_vec_mults_prev);
   end
   num_mat_vec_mults_prev = info.num_mat_vec_mults;

   
   while err > tol & itn < maxitn,
   
       if howLS==1,
          if B_is_identity
            Bp0 = p0;
          else
            Bp0=B*p0; 
          end
          pBp=real(p0'*(Bp0)); nrmp=sqrt(pBp);  p0n=(1/nrmp)*p0;
          Ap0n=A(p0n)+shift*(BU*(BU'*p0n));
					info.num_mat_vec_mults = info.num_mat_vec_mults + size(p0n, 2);
          rhos=[rho p0n'*(Ap0n)];    
          LSopts.rhos=rhos;
          
          [x1, t, lam]=Li_dir.LNSRCHg(Ax,Bx,x0,p0n,LSopts);
       else
          q=p0-x0*(Bx'*p0);  % now q'*B*x0=0;
          if B_is_identity
            Bq = q;
          else
            Bq=B*q; 
          end
          pBp=real(q'*(Bq)); nrmp=sqrt(pBp);  q=(1/nrmp)*q;
   
          Aq=A(q)+shift*(BU*(BU'*q)); tmp=x0'*Aq;
					info.num_mat_vec_mults = info.num_mat_vec_mults + size(q, 2);
          Rho1=[rho     tmp;
              conj(tmp) real(q'*Aq)];
       
         [v,lam]=Li_dir.mineig(Rho1,1);
         x1=[x0,q]*v; % x1'*B*x1 approx 1 
       end
   
       if B_is_identity
         Bx = x1;
       else
         Bx=B*x1; 
       end
       xBx=real(x1'*Bx); nrmx=sqrt(xBx);
       x1=(1/nrmx)*x1; Bx=(1/nrmx)*Bx;
       
       Ax=A(x1)+shift*(BU*(BU'*x1)); rho=lam; % rho=x1'*Ax;
	     info.num_mat_vec_mults = info.num_mat_vec_mults + size(x1, 2); 
       
       r1=Ax-rho*Bx;
       nrmr1=norm(r1);
       err=nrmr1/(nrmAB(1)+abs(rho)*nrmAB(2));

       hsty=[hsty; rho err];

       if opts.display && mod(itn, 10) == 0
          PrintIter(itn, err, info.num_mat_vec_mults - num_mat_vec_mults_prev);
       end
       num_mat_vec_mults_prev = info.num_mat_vec_mults;
       
       if whichbeta == 1,
          beta=nrmr1*nrmr1/(nrmr0*nrmr0); 
       else
          tmp=nrmx*r1;
          beta=r1'*(tmp-r0)/(nrmr0*nrmr0); 
       end
       p1=-r1+beta*p0;
       
       x0=x1;  r0=r1; p0=p1; nrmr0=nrmr1;
       if do_iter_cache
         iter_cache.V = [iter_cache.V, x0];
         iter_cache.lam = [iter_cache.lam, lam];
       end
       itn = itn+1;
   end

end

info.itn = itn; y=x0; lam=rho; res=err;

end % end function CGgS


function [] = PrintBanner()
  fprintf('\n        CG Prototype Solver\n');
  fprintf('             rel        num \n')
  fprintf(' Iter |     err 1       A*x \n');
  fprintf('------|-------------------------------------------------\n');
end


function [] = PrintIter(itn, err, num_mat_vec_mults)
   fprintf(' %3i  |  %1.5e   %4i     \n', itn, err, num_mat_vec_mults);
end

