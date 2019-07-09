function Py = project(y,b,e,s,varargin)
% function Py = project(y,b,e,s,varargin)
%
% Projects onto b'*y - e*||y||_norm >= s
%
% `norm` options:
%     2        2-norm (default)
%     inf      infinity-norm

ip = inputParser;

% Sets generated data parameters
ip.addParameter('proj_norm', 2); % can be `2` or `inf`
ip.addParameter('verbosity', false);
ip.parse(varargin{:});

proj_norm = ip.Results.proj_norm;
verbosity = ip.Results.verbosity;

%proj_norm = Inf;

if proj_norm == 2
   % Project y onto real(dot(y,b)) - e * norm(y) >= s
   tol = 1e-12;
   if (nargin < 4) || isempty(s)
      s = 1;
   end
   if (nargin < 3) || isempty(e)
      e = 0;
   end
   yb = real(dot(y(:),b(:)));
   yy = dot(y(:),y(:));
   bb = dot(b(:),b(:));
   nb = realsqrt(bb);
   if e <= tol
      Py = y+(max(0,s-yb)/nb)*(b/nb); % project onto halfspace
   %       Py = y+((s-yb)/nb)*(b/nb);      % project onto hyperplane
   elseif yb-e*realsqrt(yy)>=s-tol
      Py = y;
   else
   %       ee  = e*e;
   %       c   =  [0,e*s*conv([e,s-yb],[e,s-yb]),0] ...
   %             +[0,2*s*conv([bb-ee,-e*s],[e,s-yb]),0] ...
   %             -[conv([bb-ee,-e*s],conv([e,s-yb],[e,s-yb])),0] ...
   %             +conv([-1,0,yy],conv([bb-ee,-e*s],[bb-ee,-e*s]));
   %       nPy = roots(c);
   %       nPy = real(nPy((real(nPy)>=-tol)&(abs(imag(nPy))<=tol)));
   %       if isempty(nPy)
   %          error('[project] No real positive roots!');
   %       end
   %       l   = (s-yb+e*nPy)./(bb-ee-e*s./nPy);
   %       nPy = nPy(l>0);
   %       l   = l(l>0);
   %       if isempty(l)
   %          keyboard;
   %          error('[project] No positive multipliers!');
   %       elseif numel(l) > 1
   %          keyboard;
   %          error('[project] Multiple positive multipliers!');
   %       end
   %       Py = y+l*b; Py = nPy*(Py/norm(Py(:)));

      ee = e*e;
      c  =  [0,e*s*conv([e,s-yb],[e,s-yb]),0] ...
            +[0,2*s*conv([bb-ee,-e*s],[e,s-yb]),0] ...
            -[conv([bb-ee,-e*s],conv([e,s-yb],[e,s-yb])),0] ...
            +conv([-1,0,yy],conv([bb-ee,-e*s],[bb-ee,-e*s]));
   %       if any(isnan(c)) || any(~isfinite(c))
   %          keyboard;
   %       end
      r  = roots(c);
      n  = real(r((real(r)>0)&(abs(imag(r))<tol*abs(r))));
      if isempty(n)
         error('[project] No real positive roots!');
      end
      l  = (s-yb+e*n)./(bb-ee-e*s./n);
      n  = n(l>0);
      l  = l(l>0);
      [n,i] = max(n);
      l = l(i(1));
      if isempty(l)
         keyboard;
         error('[project] No positive multipliers!');
      elseif numel(l) > 1
         keyboard;
         error('[project] Multiple positive multipliers!');
      end
      [m1, m2, m3] = size(y);
      Py = y(:)+l*b(:); Py = n*(Py/norm(Py(:)));
      Py = reshape(Py, [m1, m2, m3]);
   end
   
elseif proj_norm == Inf
   
   tol = 1e-12;
   if (nargin < 4) || isempty(s)
      s = 1;
   end
   if (nargin < 3) || isempty(e)
      e = 0;
   end
   yb = real(dot(y(:),b(:)));
   y_inf = norm(y(:), inf);
   bb = dot(b(:),b(:));
   nb = realsqrt(bb);
   
   if e <= tol
      Py = y+(max(0,s-yb)/nb)*(b/nb); % project onto halfspace
   %       Py = y+((s-yb)/nb)*(b/nb);      % project onto hyperplane
   
   elseif yb-e*realsqrt(y_inf)>=s-tol
      Py = y;
      
   else   
      [m1, m2, m3] = size(y);
      if m2 ~= 1 || m3 ~= 1
         y = y(:);
      end
      [m1, m2, m3] = size(b);
      if m2 ~= 1 || m3 ~= 1
         b = b(:);
      end
      [m, ~] = size(b);

      if verbosity
         cvx_begin
      else
         cvx_begin quiet
      end
      variable Proj_y(m)
      minimize norm(Proj_y - y, 2)
      subject to
         b'*Proj_y - e*norm(Proj_y , inf) >= 1;
      cvx_end
      Py = reshape(Proj_y, [m1, m2, m3]);
   end
   
end
   
   
end
