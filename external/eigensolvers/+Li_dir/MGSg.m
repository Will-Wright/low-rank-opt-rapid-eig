 function  V=MGSg(B,X,U)
%
%---------------------------------------------------------------
%  Modified Gram-Schmit on X with selective re-orthogonalization
%  with respect to B-inner product
%---------------------------------------------------------------
%
%  Input:
%
%      B      (n-by-n) Hermitian and positive definite
%      X      (n-by-k)
%      U      (n-by-k0) U'*B*U=I already. Optional argument
%             Possibly k0=0, i.e., U is empty array
%
%  Output:
%
%      V      n-by-nV  B-Orthogonalized vectors from columns of X
%
%                                 Copyright by Ren-Cang Li, 02/20/2013
%---------------------------------------------------------------


% NOTE: Removed B for phase retrieval problem


[n,k]=size(X);
rtol_re=1.0e-4; % relative tolerance to perform reorthogonalization
%
V = zeros(n,k); nV=1;
if nargin==2,

   nrm2 = sqrt(real(X(:,1)'*X(:,1))); 
   if nrm2>0,
      V(:,nV)=X(:,1)/nrm2; nV=nV+1;
   end
   for j=2:k,
     vh = X(:,j); nrm2 = sqrt(real(vh'*vh));
     tolorth=rtol_re*nrm2;
     %  by MGS
     for i=1:nV-1,
           vh = vh - V(:,i)*( V(:,i)'*vh );
     end
     nrm2=sqrt(real(vh'*vh));
     %  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
     %  Perform re-orthogonalization once by MGS when deemed necessary.
     if nrm2 <= tolorth
        for i=1:nV-1,
           vh = vh - V(:,i)*( V(:,i)'*vh );
        end
        nrm2=sqrt(real(vh'*vh));
     end
     if nrm2>0
        V(:,nV)=vh/nrm2; %R(j,j)=nrm2;
        nV=nV+1;
     end 
   end

elseif nargin==3,

   k0=size(U,2);
   vh = X(:,1); nrm2 = sqrt(real(vh'*vh));
   tolorth=rtol_re*nrm2;
   % by MGS
   for i=1:k0
       vh = vh - U(:,i)*( ( U(:,i) )'* vh );
   end
   nrm2=sqrt(real(vh'*vh));
   %  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
   %  Perform re-orthogonalization once by MGS when deemed necessary.
   if nrm2 < tolorth    
      for i=1:k0
          vh = vh - U(:,i)*( ( U(:,i) )'*vh );
      end
      nrm2=sqrt(real(vh'*vh)); 
   end
   if nrm2 > 0,
      V(:,nV)=vh/nrm2; nV=nV+1;
   end 
   
   for j=2:k,
     vh = X(:,j); nrm2 = sqrt(real(vh'*vh));
     tolorth=rtol_re*nrm2;
     %  by MGS
     for i=1:k0
         vh = vh - U(:,i)*( ( U(:,i) )'*vh );
     end
     for i=1:nV-1,
           vh = vh - V(:,i)*( V(:,i)'*vh );
     end
     nrm2=sqrt(real(vh'*vh));
     %  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
     %  Perform re-orthogonalization once by MGS when deemed necessary.
     if nrm2 <= tolorth
        for i=1:k0
            vh = vh - U(:,i)*( ( U(:,i) )'*vh );
        end
        for i=1:nV-1,
           vh = vh - V(:,i)*( V(:,i)'*vh );
        end
        nrm2=sqrt(real(vh'*vh));
     end
     if nrm2>0
        V(:,nV)=vh/nrm2; %R(j,j)=nrm2;
        nV=nV+1;
     end
   end

end

V=V(:,1:nV-1);
