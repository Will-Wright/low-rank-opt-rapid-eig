function [Proj_y] = project_l_infty(y, b, epsilon)
% function [Proj_y] = project_l_infty(y, b, epsilon)
%
% projects onto b'*y - ep*||y||_inf >= 1
%

verbosity = 1;
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
      b'*Proj_y - epsilon*norm(Proj_y , inf) >= 1;
cvx_end

end