function [lam, y, res, hsty, info]=PCGgS(A_in, n, B, x, BU, nrmAB, shift, tol, maxitn, opts)
%
%  Preconditioned Conjugate Gradient Method to compute the kth smallest eigenvalue of A-lamda B,
%  where A and B are Hermitian, and B is positive definite, k=k0+1 and k0 is the number 
%  of columns in U.
%
%  Deflation is done by shifting away known eigenvalues. Basically CG on
%
%      A+shift*B*U*(B*U)'=A+shift*BU*BU', B
%
% Copyright by Ren-Cang Li, 7/21/2013
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
%        info  the number of CG steps
%
%----------------------------------------------------------------------
%

%n=max(size(A));
%maxitn=round(0.1*n); 

whichbeta = 1; 
k0=size(BU,2); k=k0+1; 
LSopts.Line=1;
           
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
x0=(1/nrmx)*x; Bx=(1/nrmx)*Bx;

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

   Ax=A(x0); rho=x0'*Ax;
   info.num_mat_vec_mults = info.num_mat_vec_mults + size(x0, 2); 
   r0=Ax-rho*Bx; 
   
   nrmr0=norm(r0);
   err=nrmr0/(nrmAB(1)+abs(rho)*nrmAB(2)); itn=0;     
   
   hsty=[rho err];

   %p0=-r0;    
   
   if opts.precond==1
      if opts.precond_one==1
         % C*Kr=r0
         Kr=U\(L\r0(pp));
      elseif opts.precond_one==2
         CGx0=zeros(n,1);
         [Kr, error, iter, flag, info_cg] = Li_dir.LinCG(C, n, CGx0, r0, CG_M, CGopts);
				 info.num_mat_vec_mults = info.num_mat_vec_mults + info_cg.num_mat_vec_mults; 
      elseif opts.precond_one==3
         % TBI
      end
   elseif opts.precond==2
      % (L*L')*Kr=r0
      Kr=(L')\(L\r0);
   end 
   p0=-Kr;
   rKr0=real(r0'*Kr);

	 if opts.display
     PrintIter(0, err, info.num_mat_vec_mults - num_mat_vec_mults_prev);
   end
   num_mat_vec_mults_prev = info.num_mat_vec_mults;   

   while err > tol & itn < maxitn,
   
      Bp0=B*p0; pBp=real(p0'*(Bp0)); nrmp=sqrt(pBp);  p0n=(1/nrmp)*p0;
      Ap0n=A(p0n);
      info.num_mat_vec_mults = info.num_mat_vec_mults + size(p0n, 2);
      rhos=[rho p0n'*(Ap0n)];    
      LSopts.rhos=rhos;
      
      [x1, t, lam]=Li_dir.LNSRCHg(Ax,Bx,x0,p0n,LSopts);

      Bx=B*x1; xBx=real(x1'*Bx); nrmx=sqrt(xBx);
      x1=(1/nrmx)*x1; Bx=(1/nrmx)*Bx;
       
      Ax=A(x1); rho=lam; %rho=x1'*Ax;
			info.num_mat_vec_mults = info.num_mat_vec_mults + size(x1, 2); 
       
      r1=Ax-rho*Bx;
      nrmr1=norm(r1);
      err=nrmr1/(nrmAB(1)+abs(rho)*nrmAB(2));

      hsty=[hsty; rho err];
          
      if opts.precond==1
         if opts.precond_one==1
            % C*Kr=r1
            Kr=U\(L\r1(pp));
         elseif opts.precond_one==2
            CGx0=zeros(n,1);
            [Kr, error, iter, flag, info_cg] = Li_dir.LinCG(C, n, CGx0, r1, CG_M, CGopts); 
					  info.num_mat_vec_mults = info.num_mat_vec_mults + info_cg.num_mat_vec_mults;
         elseif opts.precond_one==3
            % TBI
         end
      elseif opts.precond==2
         % (L*L')*Kr=r1
         Kr=(L')\(L\r1);
      end 
      rKr1=real(r1'*Kr);
       
      if whichbeta == 1,
         beta=rKr1/(rKr0); 
      else
         tmp=nrmx*r1;
         beta=Kr'*(tmp-r0)/(rKr0*rKr0); 
      end
      p1=-Kr+beta*p0;
      
      x0=x1; r0=r1; p0=p1; rKr0=rKr1; 
      
			if opts.display
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
   
   Ax=A(x0)+shift*(BU*(BU'*x0)); rho=x0'*Ax;
   info.num_mat_vec_mults = info.num_mat_vec_mults + size(x0, 2); 
   r0=Ax-rho*Bx; 
   
   nrmr0=norm(r0);
   err=nrmr0/(nrmAB(1)+abs(rho)*nrmAB(2)); 
   
   itn=0;
   
   hsty=[rho err];

   %p0=-r0;   
   
    if opts.precond==1
       if opts.precond_one==1
          % [C+shift*BU*(BU)']*Kr=r0
          Cinvr=U\(L\r0(pp));
          Kr=Cinvr-CiViT*(CiV'*r0);
       elseif opts.precond_one==2
          CGx0=zeros(n,1);
          [Kr, error, iter, flag, info_cg] = Li_dir.LinCG(C, n, CGx0, r0, CG_M, CGopts);
				  info.num_mat_vec_mults = info.num_mat_vec_mults + info_cg.num_mat_vec_mults;
       elseif opts.precond_one==3
          % TBI
       end
    elseif opts.precond==2
       if opts.precond_two==1,
          % (L*L')Kr=r0
          Kr=(L')\(L\r0);
       elseif opts.precond_two==2
          % [LL'+shift*B*U*(B*U)]Kr=r0
          Cinvr=(L')\(L\r0);
          Kr=Cinvr-CiViT*(CiV'*r0);
       end
    end 
%   Kr=r0;
   p0=-Kr;
   rKr0=real(r0'*Kr);
  
	 if opts.display
		  PrintIter(0, err, info.num_mat_vec_mults - num_mat_vec_mults_prev);
	 end
	 num_mat_vec_mults_prev = info.num_mat_vec_mults;
 
   while err > tol & itn < maxitn,
   
      Bp0=B*p0; pBp=real(p0'*(B*p0)); nrmp=sqrt(pBp);  p0n=(1/nrmp)*p0;
      Ap0n=A(p0n)+shift*(BU*(BU'*p0n));
			info.num_mat_vec_mults = info.num_mat_vec_mults + size(p0n, 2);
      rhos=[rho p0n'*(Ap0n)];    
      LSopts.rhos=rhos;
       
      [x1, t, lam]=Li_dir.LNSRCHg(Ax,Bx,x0,p0n,LSopts);

      Bx=B*x1; xBx=real(x1'*Bx); nrmx=sqrt(xBx);
      x1=(1/nrmx)*x1; Bx=(1/nrmx)*Bx;
      
      Ax=A(x1)+shift*(BU*(BU'*x1)); rho=lam; % rho=x1'*Ax;
			info.num_mat_vec_mults = info.num_mat_vec_mults + size(x1, 2); 
       
      r1=Ax-rho*Bx;
      nrmr1=norm(r1);
      err=nrmr1/(nrmAB(1)+abs(rho)*nrmAB(2));

      hsty=[hsty; rho err];
   
      if opts.precond==1
         if opts.precond_one==1
            % [C+shift*BU*(BU)']*Kr=r1
            Cinvr=U\(L\r1(pp));
            Kr=Cinvr-CiViT*(CiV'*r1);
         elseif opts.precond_one==2
            CGx0=zeros(n,1);
            [Kr, error, iter, flag, info_cg] = Li_dir.LinCG(C, n, CGx0, r1, CG_M, CGopts);
						info.num_mat_vec_mults = info.num_mat_vec_mults + info_cg.num_mat_vec_mults;
         elseif opts.precond_one==3
            % TBI
         end
      elseif opts.precond==2
         if opts.precond_two==1,
            % (L*L')Kr=r1
            Kr=(L')\(L\r1);
         elseif opts.precond_two==2
            % [LL'+shift*B*U*(B*U)]Kr=r1
            Cinvr=(L')\(L\r1);
            Kr=Cinvr-CiViT*(CiV'*r1);
         end
      end    
%      Kr=r1;
      rKr1=real(r1'*Kr);
       
      if whichbeta == 1,
         beta=rKr1/rKr0; 
      else
         tmp=nrmx*r1;
         beta=Kr'*(tmp-r0)/rKr0; 
      end
      p1=-Kr+beta*p0;
       
      x0=x1; r0=r1; p0=p1; rKr0=rKr1; 

			if opts.display
				 PrintIter(itn, err, info.num_mat_vec_mults - num_mat_vec_mults_prev);
			end
			num_mat_vec_mults_prev = info.num_mat_vec_mults;
  
      itn = itn+1;
   end

end

info.itn = itn; y=x0; lam=rho; res=err;

end % end function PCGgS


function [] = PrintBanner()
	fprintf('\n        PCG Prototype Solver\n');
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






