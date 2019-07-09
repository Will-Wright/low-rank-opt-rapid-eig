function [num_eigs_update] = adaptive_eigs_params(num_matvecs, num_eigs, iter, opts)
% [num_eigs] = adaptive_eigs_params(num_eigs, num_matvecs)
% Determine new num_eigs (parameter `k` for `eigs`)

num_eigs_min = opts.num_eigs_min;
num_eigs_max = opts.num_eigs_max;
num_lin_interp_pts = opts.num_lin_interp_pts;

num_matvecs = num_matvecs(:);
num_eigs = num_eigs(:);

if length(num_eigs) == 1 || num_eigs(iter) == num_eigs_min
   num_eigs_update = num_eigs(iter,1) + 1;
elseif num_eigs(iter) == num_eigs_max
   num_eigs_update = num_eigs(iter,1) - 1;
elseif iter < num_lin_interp_pts
   delta_2 = sign( num_eigs(iter, 1) - num_eigs(iter-1, 1) )...
      *sign( num_matvecs(iter-1, 1) - num_matvecs(iter, 1) );
   if delta_2 == 0
      delta_2 = sign( num_eigs(iter, 1) - num_eigs(iter-1, 1) );
   end
   num_eigs_update = num_eigs(iter, 1) + delta_2;
else
   delta_2 = sign( num_eigs(iter, 1) - num_eigs(iter-1, 1) )...
      *sign( num_matvecs(iter-1) - num_matvecs(iter) );
   if delta_2 == 0
      delta_2 = sign( num_eigs(iter, 1) - num_eigs(iter-1, 1) );
   end
      
   A = [ones(num_lin_interp_pts, 1), num_eigs(iter - num_lin_interp_pts + 1 : iter, 1)];
   b = num_matvecs(iter - num_lin_interp_pts + 1 : iter, 1);
   x = A\b;
   delta_4 = -sign(x(2,1));
      
   if delta_2 == delta_4
      delta = 2*delta_2;
   else
      delta = delta_2;
   end
   
   num_eigs_update = num_eigs(iter, 1) + delta;
   num_eigs_update = min( max(num_eigs_update, num_eigs_min), num_eigs_max);
end


