function [lam, y, res, hsty, info]=PSDgS(A_in, n, B, x, BU, nrmAB, shift, tol, maxitn, opts)
%
%  Preconditioned Steepest Decent Method to compute the kth smallest eigenvalue of A-lamda B,
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
%        opts  options ...
%              opts.precond=1,  use 1st type preconditioner K=(A-sgm B)^{-1}
%                   1) opts.precond_one=1: direct solver
%                   2) opts.precond_one=2: A-sgm B is (effectively) SPD. K is implemented by CG   
%                   3) opts.precond_one=3: A-sgm B is indefinite. K is implemented by MINRES 
%              opts.precond=2,  use 2nd type preconditioner; suppose A-sgm B approx LDL'.
%                   1) opts.precond_two=1: K=(LL')^{-1}
%                   2) opts.precond_two=2: K=[LL'+shift*B*U*(B*U)]^{-1}
%
%              opts.precond_sgm  the shift parameter to build a preconditioner
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

sgm=opts.precond_sgm;

if isnumeric(A_in)
 A = @(x) A_in*x;
 C = @(x) A_in*x - (sgm)*(B*x);
else
 A = A_in;
 C = @(x) A_in(x) - (sgm)*(B*x);
end

if ~isfield(opts, 'display') || opts.display ~=1
 opts.display = 0;
end
if opts.display
 PrintBanner();
end
info.num_mat_vec_mults = 0;
num_mat_vec_mults_prev = 0;

k0=size(BU,2);  
LSopts.Line=1;   % Line Search opts 
howLS=0;         % =0 solve 2-by-2 eigenvalue problem
                 % =1 use Line Search function LNSRCHg(...)
           
if opts.precond==1
   if opts.precond_one==1
      [L,U,pp]=lu(C,'vector');
   elseif opts.precond_one==2
      CGopts.nitn=10;
      CGopts.tol=1e-2;
%      CGopts.met=2;
			CGopts.met = 0;
      CG_M=0;
   elseif opts.precond_one==3
      % TBI
   end
elseif opts.precond==2
   iLUopts.thresh = 0;
   iLUopts.udiag = 1;
   iLUopts.milu = 1;
   [L, U] = luinc(C, iLUopts);   
   %  Scale diagonals to get L.
   d = diag (U);
   d = sqrt(abs(d));
   for i = 1:n
      if (d(i) < 1e-2)
         d(i) = 1.0;
      end
   end
   L = L * spdiags(d, 0, n, n);
end  

Bx=B*x; xBx=real(x'*Bx); nrmx=sqrt(xBx);
x=(1/nrmx)*x; Bx=(1/nrmx)*Bx;

if k0==0,

   if opts.precond==1
      if opts.precond_one==1
         % no action
      elseif opts.precond_one==2
         CGopts.update = 0;
      elseif opts.precond_one==3
             % TBI
      end
   elseif opts.precond==2
      % no action
   end    

   Ax=A(x); rho=real(x'*Ax);
   info.num_mat_vec_mults = info.num_mat_vec_mults + size(x, 2); 
   
   r=Ax-rho*Bx;   
   err=norm(r)/( nrmAB(1)+abs(rho)*nrmAB(2) ); itn=0;
   
   hsty=[rho err];

   if opts.display
			PrintIter(0, err, info.num_mat_vec_mults - num_mat_vec_mults_prev);
   end
	 num_mat_vec_mults_prev = info.num_mat_vec_mults;

   while err > tol & itn < maxitn, 
   
       if opts.precond==1
          if opts.precond_one==1
             % C*p=r
             p=U\(L\r(pp));
          elseif opts.precond_one==2
             CGx0=zeros(n,1);
             [p, error, iter, flag, info_cg] = Li_dir.LinCG(C, n, CGx0, r, CG_M, CGopts);
             info.num_mat_vec_mults = info.num_mat_vec_mults + info_cg.num_mat_vec_mults;
          elseif opts.precond_one==3
             % TBI
          end
       elseif opts.precond==2
          % (L*L')p=r
          p=(L')\(L\r);
       end   
       
       if howLS==1,
          Bp=B*p; pBp=real(p'*Bp); nrmr=sqrt(pBp); p=(1/nrmr)*p;
          rhos=[rho p'*(A(p))];
				  info.num_mat_vec_mults = info.num_mat_vec_mults + size(p, 2);
          LSopts.rhos=rhos;
       
          [y, t, lam]=Li_dir.LNSRCHg(Ax,Bx,x,p,LSopts);
       else
          q=p-x*(Bx'*p);  % now q'*B*x=0;
          Bq=B*q; pBp=real(q'*(Bq)); nrmp=sqrt(pBp);  q=(1/nrmp)*q;
   
          Aq=A(q); tmp=x'*Aq;
					info.num_mat_vec_mults = info.num_mat_vec_mults + size(q, 2);
          Rho1=[rho     tmp;
              conj(tmp) real(q'*Aq)];
       
         [v,lam]=Li_dir.mineig(Rho1,1);
         y=[x,q]*v; % y'*B*y approx 1 
       end
       
       x=y;
       Bx=B*x; xBx=real(x'*Bx); nrmx=sqrt(xBx);
       x=(1/nrmx)*x; Bx=(1/nrmx)*Bx;
       Ax=A(x); rho=lam; % rho=x'*Ax; 
		   info.num_mat_vec_mults = info.num_mat_vec_mults + size(x, 2);
       
       r=Ax-rho*Bx; 
       %disp([rho, lam, rho-lam])
       err=norm(r)/( nrmAB(1)+abs(rho)*nrmAB(2) );
       
       hsty=[hsty; rho err];

       if opts.display && mod(itn, 10) == 0
         PrintIter(itn, err, info.num_mat_vec_mults - num_mat_vec_mults_prev);
       end
       num_mat_vec_mults_prev = info.num_mat_vec_mults;

       itn = itn+1;
   end

else  % k0>0

   if opts.precond==1
      if opts.precond_one==1
         % [C+shift*V*V']^{-1} = C^{-1}-shift*C^{-1}*V*[I+shift*V'*C^{-1}*V]^{-1}*V'*C^{-1}, where V=BU
         CiV=U\(L\BU(pp,:));
         T=eye(k0)+shift*BU'*CiV; invT=inv(T);
         CiViT=CiV*(shift*invT);
      elseif opts.precond_one==2
         CGopts.update = 1;
         CGopts.shift =shift;
         CGopts.V = BU;
      elseif opts.precond_one==3
             % TBI
      end
   elseif opts.precond==2
      if opts.precond_two==2
         % [C+shift*V*V']^{-1} = C^{-1}-shift*C^{-1}*V*[I+shift*V'*C^{-1}*V]^{-1}*V'*C^{-1}, where V=BU, C approx LL'
         CiV=(L')\(L\BU);
         T=eye(k0)+shift*BU'*CiV; invT=inv(T);
         CiViT=CiV*(shift*invT);
      end
   end    

   Ax=A(x)+shift*(BU*(BU'*x)); rho=real(x'*Ax);
	 info.num_mat_vec_mults = info.num_mat_vec_mults + size(x, 2); 
   
   r=Ax-rho*Bx;   
   err=norm(r)/( nrmAB(1)+abs(rho)*nrmAB(2) ); itn=0;
   
   hsty=[rho err];

 	 if opts.display
			PrintIter(0, err, info.num_mat_vec_mults - num_mat_vec_mults_prev);
	 end
	 num_mat_vec_mults_prev = info.num_mat_vec_mults;


   while err > tol & itn < maxitn, 
   
       if opts.precond==1
          if opts.precond_one==1
             % [C+shift*BU*(BU)']*p=r
             Cinvr=U\(L\r(pp));
             p=Cinvr-CiViT*(CiV'*r);
          elseif opts.precond_one==2
             CGx0=zeros(n,1);
             [p, error, iter, flag, info_cg] = Li_dir.LinCG(C, n, CGx0, r, CG_M, CGopts);
             info.num_mat_vec_mults = info.num_mat_vec_mults + info_cg.num_mat_vec_mults;
          elseif opts.precond_one==3
             % TBI
          end
       elseif opts.precond==2
          if opts.precond_two==1,
             % (L*L')p=r
             p=(L')\(L\r);
          elseif opts.precond_two==2
             % [LL'+shift*B*U*(B*U)]p=r
             Cinvr=(L')\(L\r);
             p=Cinvr-CiViT*(CiV'*r);
          end
       end    
              
       if howLS==1,
          % This option has difficulty reducing residual below abount 1e-12 on an example,
          % But the other option doesn't. So I think the function LNSRCHg(...) is likely
          % the reason.  
          Bp=B*p; pBp=real(p'*Bp); nrmr=sqrt(pBp); p=(1/nrmr)*p;
          Ap=A(p)+shift*(BU*(BU'*p));
					info.num_mat_vec_mults = info.num_mat_vec_mults + size(p, 2);
          rhos=[rho p'*(Ap)];
          LSopts.rhos=rhos;
       
          [y, t, lam]=Li_dir.LNSRCHg(Ax,Bx,x,p,LSopts);
       else
          q=p-x*(Bx'*p);  % now q'*B*x=0;
          Bq=B*q; pBp=real(q'*(Bq)); nrmp=sqrt(pBp);  q=(1/nrmp)*q;
   
          Aq=A(q)+shift*(BU*(BU'*q)); tmp=x'*Aq;
			    info.num_mat_vec_mults = info.num_mat_vec_mults + size(q, 2);
          Rho1=[rho     tmp;
              conj(tmp) real(q'*Aq)];
       
         [v,lam]=Li_dir.mineig(Rho1,1);
         y=[x,q]*v; % y'*B*y approx 1 
       end
       
%       Bp=B*p; pBp=real(p'*Bp); nrmr=sqrt(pBp); p=(1/nrmr)*p;
%       Ap=A*p+shift*(BU*(BU'*p));
%       rhos=[rho p'*(Ap)];
%       LSopts.rhos=rhos;
%       
%       [y, t, lam]=LNSRCHg(Ax,Bx,x,p,LSopts);
       
       x=y;        
       Bx=B*x; xBx=real(x'*Bx); nrmx=sqrt(xBx);
       x=(1/nrmx)*x; Bx=(1/nrmx)*Bx;
       Ax=A(x)+shift*(BU*(BU'*x)); rho=lam; %rho=x'*Ax;
			 info.num_mat_vec_mults = info.num_mat_vec_mults + size(x, 2); 
       
       r=Ax-rho*Bx; 
       err=norm(r)/( nrmAB(1)+abs(rho)*nrmAB(2) );
       
       hsty=[hsty; rho err];

			 if opts.display && mod(itn, 10) == 0
					PrintIter(itn, err, info.num_mat_vec_mults - num_mat_vec_mults_prev);
			 end
			 num_mat_vec_mults_prev = info.num_mat_vec_mults;

       itn = itn+1;
   end

end

info.itn = itn; y=x; lam=rho; res=err;

end % end function PSDgS


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







