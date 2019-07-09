function [lam, Y, res, hsty, info]=xLOBPCGgS(A_in, n, B, X, BU, nrmAB, Nd, shift, tol, maxitn, opts)
%
%  Expertly Locally Optimal Block Preconditioned Conjugate Gradient (LOBCG) method to compute 
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
%        Nd    integer, # of eigenpairs needed to be computed by this call
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
%              opts.display      toggles display (on == 1, off == 0)
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
%              info.kc  number of converged eigenpairs by this call
%              info.num_mat_vec_mults  number of matrix-vector multiplications
%
% 
%---------------------------------------------------------------
%

%n=max(size(A)); 

if ~isfield(opts, 'display') || opts.display ~=1
  opts.display = 0;
end
if opts.display
  PrintBanner();
end

nb=size(X,2); 
k0=size(BU,2); 
if isfield(opts, 'nb_desired')
  nb_desired = opts.nb_desired;
else
  nb_desired = nb;
end
info.num_mat_vec_mults = 0;
num_mat_vec_mults_prev = 0;
%
% NOTE: info.num_mat_vec_mults is only valid/correct if using CG,
%       i.e., if user selects opts.precond = 1; opts.precond_one = 2;

kc=k0;  % # of converged eigenpairs (including pass-ins)
Y=[];

which_imp=2;

%maxitn=min(round(0.05*n),100);
          
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
      if isfield(opts, 'cg_nitn')
        CGopts.nitn = opts.cg_nitn;
      else
        CGopts.nitn=10;
      end
      if isfield(opts, 'cg_tol')
        CGopts.tol = opts.cg_tol;
      else
        CGopts.tol=1e-2;
      end
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

if k0==0
   if opts.precond==1
      if opts.precond_one==1
         CiV=[];
      elseif opts.precond_one==2
         CGopts.update = 0;
      elseif opts.precond_one==3
             % TBI
      end
   elseif opts.precond==2
      if opts.precond_two==2
         CiV=[];
      end
   end  
else % k0>0
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
end    

%X0=Li_dir.MGSg(B,X);
BX=B*X;
XBX=X'*BX; cholXBX=chol(XBX); invcholXBX=inv(cholXBX);
X0=X*invcholXBX;

AX=A(X0);
info.num_mat_vec_mults = info.num_mat_vec_mults + size(X0, 2);
if k0 > 0 
   AX=AX+shift*(BU*(BU'*X0));
end
Rho=X0'*AX; 
BX=B*X0; 
R=AX-BX*Rho;

hsty.eig=[]; hsty.res=[];
lam=[]; res=[];    % converged eigenvalues and corresponding residuals

if k0 == 0
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
else % k0 > 0
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
end 
%P=R;

Q=Li_dir.MGSg(B,P,X0);

% [X0,Q]'*A*[X0,Q]
tmp=Q'*AX; AQ=A(Q);
info.num_mat_vec_mults = info.num_mat_vec_mults + size(Q, 2);
if k0 > 0 
   AQ=AQ+shift*(BU*(BU'*Q));
end

Rho1=[Rho  tmp';
     tmp     Q'*AQ];
[V,D]=Li_dir.mineig(Rho1,2*nb);   % V  should have orthonormal columns

X1=[X0,Q]*V(:,1:nb);           % X1 should have B-orthonormal columns           
if which_imp==2
   Z=Q;
end

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
      
   if opts.precond==1
      if opts.precond_one==1
         % [C+shift*V*V']^{-1} = C^{-1}-shift*C^{-1}*V*[I+shift*V'*C^{-1}*V]^{-1}*V'*C^{-1}, where V=BU
         CiV=[CiV, U\(L\BX1(pp,1:kk))];
         T=eye(kc)+shift*BU'*CiV; invT=inv(T);
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
         CiV=[CiV, (L')\(L\BX1(:,1:kk))];
         T=eye(kc)+shift*BU'*CiV; invT=inv(T);
         CiViT=CiV*(shift*invT);
      end
   end 
   
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


if opts.display
  PrintIter(0, nrmR, info.num_mat_vec_mults, kc, kk);
end
num_mat_vec_mults_prev = info.num_mat_vec_mults;


itn=1;

while itn < maxitn && kc<k0+Nd,
    
   XBX=X1'*BX; cholXBX=chol(XBX); invcholXBX=inv(cholXBX);
   X1=X1*invcholXBX; AX=AX*invcholXBX; Rho=X1'*AX;
    
   if kc == 0
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
   else % kc > 0
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
   end 
   %P=R;
    
   if which_imp==1
      Q=Li_dir.MGSg(B,[P X0],X1);
   else
      Q=Li_dir.MGSg(B,[P Z],X1);
   end
    
   % [X1,Q]'*A*[X1,Q]
   tmp=Q'*AX; AQ=A(Q);
   info.num_mat_vec_mults = info.num_mat_vec_mults + size(Q, 2);
   if kc > 0
      AQ=AQ+shift*(BU*(BU'*Q));
   end
   Rho1=[Rho  tmp';
         tmp Q'*AQ];
       
   [V,D]=Li_dir.mineig(Rho1,2*nb);  % V  should have orthonormal columns
   X0=X1;    
    
   %XQ=[X0,Q]; mm=size(XQ,2); disp('Here 1'); disp([nb, norm(XQ'*B*XQ-eye(mm),1)]);
   X1=[X0,Q]*V(:,1:nb);      % X1 should have B-orthonormal columns
   AX1=[AX, AQ]*V(:,1:nb);   % shift is already included.    
   BX1=B*X1;
   R=AX1-BX1*diag(D(1:nb));
    
   %mm=size(V,2); disp('Here 2'); disp(norm(V'*V-eye(mm),1));
   %disp('Here 5'); disp([kk, norm(X1'*BX-eye(nb),1)]);
    
   if which_imp==2
      nV=size(V,1);
      Z=Q*V(nb+1:nV,:);
   end
    
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
   % # of converged eigenpairs is kk.
   lam=[lam; D(1:kk)]; res=[res nrmR(1:kk)]; 
   if kk>0    
      Y=[Y, X1(:,1:kk)];
      BU=[BU, BX1(:,1:kk)]; kc=kc+kk; 
      
      if opts.precond==1
         if opts.precond_one==1
            % [C+shift*V*V']^{-1} = C^{-1}-shift*C^{-1}*V*[I+shift*V'*C^{-1}*V]^{-1}*V'*C^{-1}, where V=BU
            CiV=[CiV, U\(L\BX1(pp,1:kk))];
            T=eye(kc)+shift*BU'*CiV; invT=inv(T);
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
            CiV=[CiV, (L')\(L\BX1(:,1:kk))];
            T=eye(kc)+shift*BU'*CiV; invT=inv(T);
            CiViT=CiV*(shift*invT);
         end
      end 
   
      % make working block to nb again
      X1a=[X0,Q]*V(:,nb+1:nb+kk);           % X1a should have B-orthonormal columns
      AXa=[AX, AQ]*V(:,nb+1:nb+kk);         % shift is already included.
      BXa=B*X1a;
      Ra=AXa-BXa*diag(D(nb+1:nb+kk));
  
      R=[R(:,kk+1:nb) Ra]; D=D(kk+1:kk+nb); X1=[X1(:,kk+1:nb) X1a];  
      AX=[AX1(:,kk+1:nb) AXa]; BX=[BX1(:,kk+1:nb) BXa];
      % for X0, we take what remains   
      X0=X0(:,kk+1:nb);  
      %disp('Here 3'); disp([kk, norm(X1'*B*X1-eye(nb),1)]);
   else
      AX=AX1; BX=BX1;  
      %disp('Here 4'); disp([kk, norm(X1'*BX-eye(nb),1)]);      
   end
   
   if opts.display
      PrintIter(itn, nrmR, info.num_mat_vec_mults - num_mat_vec_mults_prev, kc, kk);
   end
   num_mat_vec_mults_prev = info.num_mat_vec_mults;
   itn = itn+1;

end

info.itn = itn; info.kc=kc-k0; Y=[Y X1];

end % end function xLOBPCGgS


function [] = PrintBanner()
  fprintf('\n              xLOBPCG Prototype Solver\n');
  fprintf('             rel          rel       num   num   current\n')
  fprintf(' Iter |     err 1        err 2      A*x  convgd  convgd\n');
  fprintf('------|-------------------------------------------------\n');
end


function [] = PrintIter(itn, nrmR, num_mat_vec_mults, kc, kk)
 if length(nrmR) >= 2
    print_array = [itn, nrmR, num_mat_vec_mults, kc, kk];
    print_str = ' %3i  |  ';
    for i = 1:length(nrmR)
      print_str = strcat(print_str, '  %1.5e');
    end
    print_str = strcat(print_str, '  %4i      %1i      %1i\n');
    fprintf(print_str, print_array);
 else

 end
end
