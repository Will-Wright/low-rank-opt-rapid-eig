function [lam, Y, res, hsty, info]=BPSDgS(A_in, n, B, X, BU, nrmAB, shift, tol, maxitn, opts)
%
%  Block Preconditioned Steepest Decent Method to compute the kth to (k+nb-1)st smallest eigenvalues of A-lamda B,
%  where A and B are Hermitian, and B is positive definite, k=k0+1 and k0 is the number 
%  of columns in BU=B*U, size(X)=[n,nb].
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
%        X     n-by-nb matrix initial guess to the eigenspace asscoated with the kth to (k+nb-1)st smallest eigenvalues
%        BU    n-by-k0, B*U, where U's columns span an approximate invariant subspace 
%              (accurate enough, however) associated with 1st to k0th smallest 
%              eigenvalues. Assume U'*B*U=I.
%              Possibly k0=0, i.e., BU is empty array.
%        nrmAB 2-by-1, nrmAB(1): estimate of ||A||
%                      nrmAB(2): estimate of ||B||
%        shift real for shifting away known eigenvalues
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
%maxitn=round(0.05*n); 
if ~isfield(opts, 'display') || opts.display ~=1
 opts.display = 0;
end
if opts.display
 PrintBanner();
end
info.num_mat_vec_mults = 0;
num_mat_vec_mults_prev = 0;

k0=size(BU,2); nb=size(X,2);
if isfield(opts, 'nb_desired')
  nb_desired = opts.nb_desired;
else
  nb_desired = nb;
end

sgm=opts.precond_sgm;

if isnumeric(A_in)
 A = @(x) A_in*x;
 C = @(x) A_in*x - sgm*(B*x);
else
 A = A_in;
 C = @(x) A_in(x) - sgm*(B*x);
end

if opts.precond==1
   if opts.precond_one==1
      [L,U,pp]=lu(C,'vector');
   elseif opts.precond_one==2
      CGopts.nitn=10;
      CGopts.tol=1e-2;
%      CGopts.met=2;
			CGopts.met = 0;
      CG_M=0;
      P=zeros(n,nb);
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

BX=B*X;
XBX=real(X'*BX); cholXBX=chol(XBX); invcholXBX=inv(cholXBX);
X=X*invcholXBX;   % now X'*B*X=I. can also be implemented by B-MGS
BX=BX*invcholXBX;  

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

   AX=A(X); Rho=X'*AX; 
   info.num_mat_vec_mults = info.num_mat_vec_mults + size(X, 2);
   
   R=AX-BX*Rho;  
   itn=0;

   hsty.eig=[]; hsty.res=[];
   
   for i = 1:nb
     nrmR(i) = norm(R(:, i))/(nrmAB(1)+abs(Rho(i, i)*nrmAB(2)));
   end

   err = max(nrmR(1:nb_desired));

	 if opts.display
		  PrintIter(0, nrmR, info.num_mat_vec_mults - num_mat_vec_mults_prev);
	 end
	 num_mat_vec_mults_prev = info.num_mat_vec_mults;

   while itn < maxitn && err > tol 
   
       if opts.precond==1
          if opts.precond_one==1
             % C*P=R
             P=U\(L\R(pp,:));
          elseif opts.precond_one==2
             CGx0=zeros(n,1);
             for i=1:nb,
                 [P(:,i), error, iter, flag, info_cg] = Li_dir.LinCG(C, n, CGx0, R(:,i), CG_M, CGopts);
             info.num_mat_vec_mults = info.num_mat_vec_mults + info_cg.num_mat_vec_mults;
						 end
          elseif opts.precond_one==3
             % TBI
          end
       elseif opts.precond==2
          % (L*L')*P=R
          P=(L')\(L\R);
       end   
       
       % X already with B-orthonormal columns
       Q=Li_dir.MGSg(B,P,X);        
       
       tmp=Q'*AX; AQ=A(Q);
			 info.num_mat_vec_mults = info.num_mat_vec_mults + size(Q, 2);
       Rho=[Rho  tmp';
             tmp Q'*AQ];
       [V,D]=Li_dir.mineig(Rho,nb);
       
       X=[X Q]*V; % It should be B-orthogonal
       
       BX=B*X;
       XBX=real(X'*BX); cholXBX=chol(XBX); invcholXBX=inv(cholXBX);
       X=X*invcholXBX;   % now X'*B*X=I. can also be implemented by B-MGS
       BX=BX*invcholXBX;

       % Depending on the sparsity of A, it may be cheaper to compute AX as
       % AX=[AX,AQ]*V*invcholXBX       
       AX=A(X); Rho=X'*AX;
			 info.num_mat_vec_mults = info.num_mat_vec_mults + size(X, 2);
   
       R=AX-BX*Rho; 
       nrmR=sqrt(sum(R.*conj(R)))./( nrmAB(1)+abs(D')*nrmAB(2) );
       hsty.eig=[hsty.eig; D.']; hsty.res=[hsty.res; nrmR];
       
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

   AX=A(X)+shift*(BU*(BU'*X)); Rho=X'*AX;
   info.num_mat_vec_mults = info.num_mat_vec_mults + size(X, 2); 
   
   R=AX-BX*Rho;   
   itn=0;
   
   hsty.eig=[]; hsty.res=[];
   
   for i = 1:nb
      nrmR(i) = norm(R(:, i))/(nrmAB(1)+abs(Rho(i, i)*nrmAB(2)));
   end
   err = max(nrmR(1:nb_desired));

	 if opts.display
		  PrintIter(0, nrmR, info.num_mat_vec_mults - num_mat_vec_mults_prev);
	 end
	 num_mat_vec_mults_prev = info.num_mat_vec_mults;

   while  itn < maxitn && err > tol 
   
       if opts.precond==1
          if opts.precond_one==1
             % [C+shift*BU*(BU)']*P=R
             CinvR=U\(L\R(pp,:));
             P=CinvR-CiViT*(CiV'*R);
          elseif opts.precond_one==2
             CGx0=zeros(n,1);
             for i=1:nb,
                 [P(:,i), error, iter, flag, info_cg] = Li_dir.LinCG(C, n, CGx0, R(:,i), CG_M, CGopts);
             info.num_mat_vec_mults = info.num_mat_vec_mults + info_cg.num_mat_vec_mults;
						 end
          elseif opts.precond_one==3
             % TBI
          end
       elseif opts.precond==2
          if opts.precond_two==1,
             % (L*L')*P=R
             P=(L')\(L\R);
          elseif opts.precond_two==2
             % [LL'+shift*B*U*(B*U)]*P=R
             CinvR=(L')\(L\R);
             P=CinvR-CiViT*(CiV'*R);
          end
       end    
       
       % X already with B-orthonormal columns
       Q=Li_dir.MGSg(B,P,X);        
       
       tmp=Q'*AX; AQ=A(Q)+shift*(BU*(BU'*Q));
	     info.num_mat_vec_mults = info.num_mat_vec_mults + size(Q, 2);
       Rho=[Rho  tmp';
             tmp Q'*AQ];
       [V,D]=Li_dir.mineig(Rho,nb);
       
       X=[X Q]*V; % It should be B-orthogonal
       
       BX=B*X;
       XBX=real(X'*BX); cholXBX=chol(XBX); invcholXBX=inv(cholXBX);
       X=X*invcholXBX;   % now X'*B*X=I. can also be implemented by B-MGS
       BX=BX*invcholXBX;
       
       % Depending on the sparsity of A, it may be cheaper to compute AX as
       % AX=[AX,AQ]*V*invcholXBX
       AX=A(X)+shift*(BU*(BU'*X)); Rho=X'*AX;
			 info.num_mat_vec_mults = info.num_mat_vec_mults + size(X, 2);
   
       R=AX-BX*Rho; 
       nrmR=sqrt(sum(R.*conj(R)))./( nrmAB(1)+abs(D')*nrmAB(2) );
       hsty.eig=[hsty.eig; D.']; hsty.res=[hsty.res; nrmR];
       
 			 for i = 1:nb
   				err(i) = norm(R(:, i))/(nrmAB(1)+abs(D(i, 1)*nrmAB(2)));
 			 end
%       if err(1) < tol
%          C = @(x) A_in(x) - (-0.001 + D(2))*(B*x);
%       end

 			 err = max(nrmR(1:nb_desired));

			 if opts.display
				  PrintIter(itn, nrmR, info.num_mat_vec_mults - num_mat_vec_mults_prev);
			 end
			 num_mat_vec_mults_prev = info.num_mat_vec_mults;

       itn = itn+1;
   end

end


info.itn = itn; Y=X; lam=D; res=nrmR;

end % end function BPSDgS


function [] = PrintBanner()
  fprintf('\n                 BPSD Prototype Solver\n');
  fprintf('             rel          rel       num \n')
  fprintf(' Iter |     err 1        err 2      A*x  \n');
  fprintf('------|--------------------------------------\n');
end


function [] = PrintIter(itn, nrmR, num_mat_vec_mults)
 if length(nrmR) >= 2
    print_array = [itn, nrmR, num_mat_vec_mults];
    print_str = ' %3i  |  ';
    for i = 1:length(nrmR)
      print_str = strcat(print_str, '  %1.5e');
    end
    print_str = strcat(print_str, '  %4i \n');
    fprintf(print_str, print_array);
 else

 end
end
