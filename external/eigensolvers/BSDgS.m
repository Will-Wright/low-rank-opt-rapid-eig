function [lam, Y, res, hsty, info, iter_cache]=BSDgS(A_in, n, B, X, BU, nrmAB, shift, tol, maxitn, opts)
%
%  Block Steepest Decent Method to compute the kth to (k+nb-1)st smallest eigenvalues of A-lamda B,
%  where A and B are Hermitian, and B is positive definite, k=k0+1 and k0 is the number 
%  of columns in BU=B*U, size(X)=[n,nb].
%
%  Deflation is done by shifting away known eigenvalues. Basically SD on
%
%      A+shift*B*U*(B*U)'=A+shift*BU*BU', B
%
% Copyright by Ren-Cang Li, 7/20/2013
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
%
%  Output
%
%        lam  converged eigenvalue
%        y     corresponding eigenvector
%        res   residual error  ||A y - lam B*y||
%        hsty  struct 
%              hsty.eig -- row i contains the eigenvalue approximations of ith iteration
%              hsty.res -- row i contains normalized residuals 
%                                ||A yi - lamb_i B yi||/(||A||*||y_i||+|lamb_i|*||B||*||y_i||)
%                          for ith iteration
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


k0=size(BU,2); nb=size(X,2);
if isfield(opts, 'nb_desired')
  nb_desired = opts.nb_desired;
else
  nb_desired = nb;
end

BX=B*X;
XBX=real(X'*BX); cholXBX=chol(XBX); invcholXBX=inv(cholXBX);
X=X*invcholXBX;   % now X'*B*X=I. can also be implemented by B-MGS
BX=BX*invcholXBX;

if k0==0,

   AX=A(X); Rho=X'*AX;
	 info.num_mat_vec_mults = info.num_mat_vec_mults + size(X, 2); 
   
   R=AX-BX*Rho;   
   itn=0;

	 for i = 1:nb
   	 err(i) = norm(R(:, i))/(nrmAB(1)+abs(Rho(i, i)*nrmAB(2)));
 	 end
 	 err = max(err(1:nb_desired));
   
   hsty.eig=[]; hsty.res=[];

   while itn < maxitn && err > tol  
       % X already with B-orthonormal columns
       Q=Li_dir.MGSg(B,R,X);        
       
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
       
       AX=A(X); Rho=X'*AX;
			 info.num_mat_vec_mults = info.num_mat_vec_mults + size(X, 2);
   
       R=AX-BX*Rho; 
       nrmR=sqrt(sum(R.*conj(R)))./( nrmAB(1)+abs(D')*nrmAB(2) );
       hsty.eig=[hsty.eig; D.']; hsty.res=[hsty.res; nrmR];

       err = max(nrmR(1:nb_desired));

			 if opts.display && mod(itn, 10) == 0
         PrintIter(itn, nrmR, info.num_mat_vec_mults - num_mat_vec_mults_prev);
       end
       num_mat_vec_mults_prev = info.num_mat_vec_mults;

       itn = itn+1;
   end

else  % k0>0

   AX=A(X)+shift*(BU*(BU'*X)); Rho=X'*AX;
   info.num_mat_vec_mults = info.num_mat_vec_mults + size(X, 2); 
   
   R=AX-BX*Rho;   
   itn=0;

	 for i = 1:nb
   	 err(i) = norm(R(:, i))/(nrmAB(1)+abs(Rho(i, i)*nrmAB(2)));
 	 end
 	 err = max(err(1:nb_desired));
   
   hsty.eig=[]; hsty.res=[];

   while itn < maxitn && err > tol
       % X already with B-orthonormal columns
       Q=Li_dir.MGSg(B,R,X);        
       
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
       
       AX=A(X)+shift*(BU*(BU'*X)); Rho=X'*AX;
			 info.num_mat_vec_mults = info.num_mat_vec_mults + size(X, 2);
   
       R=AX-BX*Rho; 
       nrmR=sqrt(sum(R.*conj(R)))./( nrmAB(1)+abs(D')*nrmAB(2) );
       hsty.eig=[hsty.eig; D.']; hsty.res=[hsty.res; nrmR];
   
       err = max(nrmR(1:nb_desired));

			 if opts.display && mod(itn, 10) == 0
   		   PrintIter(itn, nrmR, info.num_mat_vec_mults - num_mat_vec_mults_prev);
 			 end
 			 num_mat_vec_mults_prev = info.num_mat_vec_mults;
    
       itn = itn+1;
   end

end


info.itn = itn; Y=X; lam=D; res=nrmR;

end % end function BSDgS


function [] = PrintBanner()
 fprintf('\n                BSD Prototype Solver\n');
 fprintf('             rel          rel       num  \n')
 fprintf(' Iter |     err 1        err 2      A*x  \n');
 fprintf('------|-------------------------------------------------\n');
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
