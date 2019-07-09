function [lam, y, res, hsty, info, iter_cache]=SDgS(A_in, n, B, x, BU, nrmAB, shift, tol, maxitn, opts)
%
%  Steepest decent method to compute the kth smallest eigenvalue of A-lamda B,
%  where A and B are Hermitian, and B is positive definite, k=k0+1 and k0 is the number 
%  of columns in BU=B*U.
%
%  Deflation is done by shifting away known eigenvalues. Basically SD on
%
%      A+shift*B*U*(B*U)'=A+shift*BU*BU', B
%
% Copyright by Ren-Cang Li, 7/14/2013
%
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
%        info  number of iterations
% 
%---------------------------------------------------------------
%


%n=max(size(A));
%maxitn=round(0.2*n); 

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

iter_cache.V = [];
iter_cache.lam = [];


if ~isfield(opts, 'display') || opts.display ~=1
	opts.display = 0;
end
if opts.display
	PrintBanner();
end
info.num_mat_vec_mults = 0;
num_mat_vec_mults_prev = 0;

k0=size(BU,2); 
opts.Line=1;

if B_is_identity
  Bx = x;
else
  Bx=B*x; 
end
xBx=real(x'*Bx); nrmx=sqrt(xBx);
x=(1/nrmx)*x; Bx=(1/nrmx)*Bx;

if k0==0,

   Ax=A(x); rho=x'*Ax; 
   info.num_mat_vec_mults = info.num_mat_vec_mults + size(x, 2);   

   r=Ax-rho*Bx;   
   err=norm(r)/( nrmAB(1)+abs(rho)*nrmAB(2) ); itn=0;
   
   if opts.display
   		PrintIter(itn, err, info.num_mat_vec_mults - num_mat_vec_mults_prev);
 	 end
 	 num_mat_vec_mults_prev = info.num_mat_vec_mults;

   hsty=[rho err];

   while err > tol & itn < maxitn,    
       if B_is_identity
         Br = r;
       else
         Br=B*r; 
       end
       rBr=real(r'*Br); nrmr=sqrt(rBr); r=(1/nrmr)*r;
       rhos=[rho r'*(A(r))];
       opts.rhos=rhos;
       
       [y, t, lam]=Li_dir.LNSRCHg(Ax,Bx,x,r,opts);
       x=y; 
       if B_is_identity
         Bx = x;
       else
         Bx=B*x; 
       end
       
       xBx=real(x'*Bx); nrmx=sqrt(xBx);
       x=(1/nrmx)*x; Bx=(1/nrmx)*Bx;
       
       Ax=A(x); rho=x'*Ax;
       info.num_mat_vec_mults = info.num_mat_vec_mults + size(x, 2); 
       
       r=Ax-rho*Bx; 
       err=norm(r)/( nrmAB(1)+abs(rho)*nrmAB(2) );
       
			 if opts.display
				 PrintIter(itn, err, info.num_mat_vec_mults - num_mat_vec_mults_prev);
			 end
			 num_mat_vec_mults_prev = info.num_mat_vec_mults;

       hsty=[hsty; rho err];
       itn = itn+1;
   end

else  % k0>0  

   Ax=A(x)+shift*(BU*(BU'*x)); rho=x'*Ax; 
   info.num_mat_vec_mults = info.num_mat_vec_mults + size(x, 2);   

   r=Ax-rho*Bx;   
   err=norm(r)/( nrmAB(1)+abs(rho)*nrmAB(2) ); itn=0;
		 
	 if opts.display
		 PrintIter(0, err, info.num_mat_vec_mults - num_mat_vec_mults_prev);
	 end
	 num_mat_vec_mults_prev = info.num_mat_vec_mults;

   hsty=[rho err];
   while err > tol & itn < maxitn,    
       if B_is_identity
         Br = r;
       else
         Br=B*r; 
       end
       
       rBr=real(r'*Br); nrmr=sqrt(rBr); r=(1/nrmr)*r;
       Ar=A(r)+shift*(BU*(BU'*r));
 
       rhos=[rho r'*(Ar)];
       opts.rhos=rhos;
       
       [y, t, lam]=Li_dir.LNSRCHg(Ax,Bx,x,r,opts);
       x=y; 
       
       if B_is_identity
         Bx = x;
       else
         Bx=B*x; 
       end
       
       xBx=real(x'*Bx); nrmx=sqrt(xBx);
       x=(1/nrmx)*x; Bx=(1/nrmx)*Bx;
       Ax=A(x)+shift*(BU*(BU'*x)); rho=x'*Ax; 
       info.num_mat_vec_mults = info.num_mat_vec_mults + size(x, 2);
 
       r=Ax-rho*Bx; 
       err=norm(r)/( nrmAB(1)+abs(rho)*nrmAB(2) );
       
			 if opts.display
				 PrintIter(itn, err, info.num_mat_vec_mults - num_mat_vec_mults_prev);
			 end
			 num_mat_vec_mults_prev = info.num_mat_vec_mults;

       hsty=[hsty; rho err];
       itn = itn+1;
   end

end


info.itn = itn; y=x; lam=rho; res=nrmr;

end % end function SDgS

function [] = PrintBanner()
	fprintf('\n        PSD Prototype Solver\n');
	fprintf('             rel        num \n')
	fprintf(' Iter |     err 1       A*x  \n');
	fprintf('------|-------------------------------------------------\n');
end


function [] = PrintIter(itn, err, num_mat_vec_mults)
	if length(err) == 1
		fprintf(' %3i  |  %1.5e   %4i     \n', itn, err, num_mat_vec_mults);
	else

	end
end


