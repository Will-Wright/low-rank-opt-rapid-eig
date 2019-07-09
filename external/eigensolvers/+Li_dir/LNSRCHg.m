function [y, t, lam]=LNSRCHg(Ax,Bx,x,p,opts)
%
%  It perform line search to solve
%
%            inf_t rho(x+t*p),
%
%  where
%
%       rho(x)=x' A x/(x' B x), Rayleigh quotient of a Hermitian matrix pencil A-lamda B
%       with B positive definite
%
%
% Copyright by Ren-Cang Li, 1/30/2013
%---------------------------------------------------------------
%
%  Input
%
%        Ax     n-by-1; it is A*x
%        Bx     n-by-1; it is B*x
%        x      n-vector, x'*B*x=1
%        p      n-vector, p'*B*p=1
%        opts   struct
%               opts.rhos = [rho(x)  rho(p)]   
%               opts.Line = 0,  through diff(rho(x+t*p),t)=0
%                         = 1,  through solving [x,p]'*(A-lamda B)*[x,p]
%
%  Output
%
%        y     x+t*p, in case t=infty, y=p instead
%        t     optimal argument
%              = 0    inf_t rho(x+t*p)=rho(x):      no improvement
%        lam  new approximate eigenvalue
%

%AA=[x'; p']*A*[x p];
%if opts.rhos  == 0
%   rhos=[x'*Ax  p'*A*p]; % because x'*B*x=p'*B*p=1 upon entry!
%else
   rhos=opts.rhos;
%end

tmp=Ax'*p;
AA=[rhos(1)  tmp;
    conj(tmp) rhos(2)];
    
%BB=[x'; p']*B*[x p];
tmp=Bx'*p;
BB=[ 1    tmp;
    conj(tmp) 1];
%disp(cond(BB))
    
% solve eig(AA,BB)
tmp=BB(1,2)*BB(2,1); a=1.0-tmp;   
if abs(a) <= 100*eps*(1.0+abs(tmp))  
   % BB is singular => x and p parallel since B is definite => no improvement
   t = 0; y=x;  
   lam=rhos(1);
   return
else
   b=-AA(1,1)-AA(2,2)+AA(1,2)*BB(2,1)+BB(1,2)*AA(2,1);
   c=AA(1,1)*AA(2,2)-AA(1,2)*AA(2,1);
   tmp = b*b-4*a*c;
   % smallest eigebvalue lam is
   if a>0
      if b >= 0,
         lam = (-b-sqrt(abs(tmp)))/(2*a);
      else
         lam = 2*c/(-b+sqrt(abs(tmp)));
      end
   else  % a>0
      if b <= 0,
         lam = (-b+sqrt(abs(tmp)))/(2*a);
      else
         lam = -2*c/(b+sqrt(abs(tmp)));
      end
   end
end

if opts.Line==0 % through diff(rho(x+t*p),t)=0
   % a=AA(2,2)*(BB(1,2)+BB(2,1))-(AA(1,2)+AA(2,1))*BB(2,2);
   tmp1 = AA(2,2)*(BB(1,2)+BB(2,1)); tmp2 = (AA(1,2)+AA(2,1))*BB(2,2);
   a=tmp1-tmp2;
   if abs(a) <= 4*eps*(abs(tmp1)+abs(tmp2))
      a = 0;
   end
   
   % b=2*( AA(2,2)*BB(1,1)-AA(1,1)*BB(2,2) );
   tmp3 = AA(2,2); tmp4 = AA(1,1);
   b = 2*(tmp3-tmp4);
   if abs(b) <= 8*eps*(abs(tmp3)+abs(tmp4))
      b = 0;
   end
   
   % c=(AA(1,2)+AA(2,1))*BB(1,1)-AA(1,1)*(BB(1,2)+BB(2,1));
   tmp5 = (AA(1,2)+AA(2,1)); tmp6 = AA(1,1)*(BB(1,2)+BB(2,1));
   c = tmp5 - tmp6;
   if abs(c) <= 4*eps*(abs(tmp5)+abs(tmp6))
      c = 0;
   end
   
   if a ~= 0,
      tmp = b*b-4*a*c;
      if b <= 0,
         t = (-b+sqrt(abs(tmp)))/(2*a);
      else
         t = -2*c/(b+sqrt(abs(tmp)));
      end
      y=x+t*p;
   else  % a == 0
      if  b == 0 % a = b == 0, too
          y = x; t = 0;      % no improvement
      elseif b < 0   % a = 0, b < 0
          t=1/0; y=p;
      else % a = 0, b>0
         t=-c/b; y= x+t*p;
      end
   end
   
else % opts.Line==1, through solving [x,p]'*(A-lamda B)*[x,p]
      
   y=[-(AA(1,2)-lam*BB(1,2)); AA(1,1)-lam];
   % corresponding eigenvector y=[x,p]*y
   if abs(y(1))>0.0
      t=y(2)/y(1); y=x+t*p;
   else
      t=1/0.0; y=p;
   end
   
end
