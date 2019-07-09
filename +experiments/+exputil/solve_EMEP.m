function [data] = solve_EMEP(varargin)
%(results, yCell, eigsTolCell, itStart, itEnd, maxItsTestSolver)
%
% Solves evolving matrix eigenvalue problem with selected solver.
% Matrix sequence can be passed the following ways:
%     `results` struct from `experiments.evolvingmatricespl`
%     `A_data` struct from `experiments.exputil.convert_results_to_A_matrix_file`
%          

addpath external/eigensolvers/
import util.*;

ip = inputParser;
ip.addParameter('results', []);
ip.addParameter('A_data', []);
ip.addParameter('iter_start', 1);
ip.addParameter('iter_end', []);
ip.addParameter('solver', 'eigs');  % default eigenvalue solver set to `eigs`
                                    % valid solvers: `eigs`, `eigifp`, `bleigifp`, `LOPCG`

ip.addParameter('num_eigs', 2);  % default set to 2 eigenvalues, as required by saga_sd
ip.addParameter('block_size', 20); % used in `bleigifp` and `eigs` with default 20 in original saga_sd, unused in `eigifp`
ip.addParameter('inner_iter', 10); % used in `eigifp` and `bleigifp`, unused in `eigs`
ip.addParameter('inner_iter_max', 40); % used in `eigifp` only
ip.addParameter('eigifp_use_default_block_method', false); 
ip.addParameter('bleigifp_use_default_block_method', false);
ip.addParameter('num_eigs_min', 2); % used in `eigs` only
ip.addParameter('num_eigs_max', 30); % used in `eigs` only
ip.addParameter('eigs_use_adaptive_method', false);
ip.addParameter('max_block_size', 5*20);
ip.addParameter('iter_max_solver', 5000);

ip.addParameter('eigifp_tol_max', 1e-8); % prevents `eigifp` from failing to deflate properly

% matches all solver residuals for near-machine-precision
ip.addParameter('eigs_tol_scale_tight', 3e4); % parameter to get equivalent residuals for all solvers
ip.addParameter('eigifp_tol_scale_tight', 2e3); % parameter to get equivalent residuals for all solvers
ip.addParameter('bleigifp_tol_scale_tight', 1e3); % parameter to get equivalent residuals for all solvers
% matches all solver residuals for weaker precision results
ip.addParameter('eigifp_tol_scale_loose', 1e-1); % parameter to get equivalent residuals for all solvers
ip.addParameter('bleigifp_tol_scale_loose', 1e-1); % parameter to get equivalent residuals for all solvers
ip.addParameter('disp', 0);
ip.parse(varargin{:});

% Creates matrix class object to pass A'y matrix function to eigenvalue solver
if ~isempty(ip.Results.results)
   fprintf('Converting ''results'' struct to ''A_class'' object.\n');
   A_data = experiments.exputil.convert_results_to_A_matrix_file(ip.Results.results);
   AClassObj = experiments.exputil.A_class(A_data);
elseif ~isempty(ip.Results.A_data)
   fprintf('Converting ''A_data'' struct to ''A_class'' object.\n');
   A_data = ip.Results.A_data;
   AClassObj = experiments.exputil.A_class(A_data);
else
   fprintf('No valid input for evolving matrix eigenvalue problem.  Terminating program.\n');
   return;
end
   
% Initializes parameters based on input parser
iter_start = ip.Results.iter_start;
if isempty(ip.Results.iter_end) || (ip.Results.iter_end > size(A_data.y_matrix, 2))
   iter_end = size(A_data.y_matrix, 2);
else
   iter_end = ip.Results.iter_end;
end
solver = ip.Results.solver;
n = A_data.n;
data = [];
num_eigs = ip.Results.num_eigs;


opts = get_opts(ip);
opts.isreal = A_data.isreal;


if strcmp(solver, 'eigifp') && opts.eigifp_use_default_block_method
   opts.inner_iter = 0; % this setting makes `eigifp` run in default mode
end

fprintf('\nSolving evolving matrix eigenvalue problem with %s\n', solver);
fprintf('NOTE: Non-consecutive iters will not initialize correctly on first subiter\n\n');
fprintf('               |     eigs   matvec      num    1st eig   2nd eig\n');
fprintf('y_iter  sub_it |      tol    mults restarts    rel err   rel err\n');


% Sets flag to compute correct initial eigenvector
if iter_start > 1
   is_init_iter = true;
   iter_start = iter_start - 1;
else
   is_init_iter = false;
end
V_prev = [];
d_prev = ones(num_eigs, 1);
num_block_restarts_array = [];
num_mat_vec_array = [];


% Main loop
for y_iter = iter_start:iter_end
   
   opts.tol = A_data.eigs_tols(1, y_iter);
   opts.V_prev = V_prev;
   opts.d_prev = d_prev;
   opts.tol = A_data.eigs_tols(1, y_iter);

   Aty_fun = @(x) AClassObj.A_fun(y_iter, x);
   


   % Passes matrix function to eigenvalue solver
   [V, d, eig_solver_data] = solve_eigval_problem(Aty_fun, n, solver, opts);   
   data.eig_solver_data{y_iter} = eig_solver_data;   
   data.d{y_iter} = d;
   data.num_mat_vec{y_iter} = eig_solver_data.num_mat_vec;
   data.num_block_restarts{y_iter} = eig_solver_data.num_block_restarts;
   
   if opts.eigs_use_adaptive_method
      opts = get_eigs_param_update(opts, data, is_init_iter, y_iter);
   end
   
   %nFFTs = (2*eig_solver_data.num_mat_vec+1)*A_data.m3;      

   if ~isempty(V)
      rel_err1 = norm(Aty_fun(V(:, 1)) - d(1, 1)*V(:, 1)) / (abs(d(1,1)) + abs(d(1,1)));
      %VtV1 = abs(V(:, 1)'*V_eigs(:, 1));
      
      if length(d) >= 2
         rel_err2 = norm(Aty_fun(V(:, 2)) - d(2, 1)*V(:, 2)) / (abs(d(1,1)) + abs(d(2,1)));
      else
         rel_err2 = nan;
      end
      %VtV2 = abs(V(:, 2)'*V_eigs(:, 2));

      data.rel_err1{y_iter} = rel_err1;
      data.rel_err2{y_iter} = rel_err2;
      
      % Skips saving results on initialization iterate
      if ~is_init_iter
         num_block_restarts_array = [num_block_restarts_array; eig_solver_data.num_block_restarts];
         num_mat_vec_array = [num_mat_vec_array; eig_solver_data.num_mat_vec];
      end
      
      % Skips print on initialization iterate
      if ~is_init_iter
         print_array = [y_iter, A_data.y_is_backtrack_step(1, y_iter), ...
                        opts.tol, eig_solver_data.num_mat_vec, ...
                        eig_solver_data.num_block_restarts, rel_err1, rel_err2];
         print_str = '  %4i     %3i | %1.2e    %5i     %4i   %1.2e  %1.2e ';
         fprintf(print_str, print_array);
         fprintf('\n');
      end
   else
      if ~is_init_iter
         print_array = [y_iter, A_data.y_is_backtrack_step(1, y_iter), ...
                        opts.tol];
         print_str = '  %4i     %3i | %1.2e  ';
         fprintf(print_str, print_array);
         fprintf('Error: solver failed internally.')
         fprintf('\n');
      end
   end


   d_prev = d;
   V_prev = V;
   is_init_iter = false;
end


avg_num_block_restarts = mean(num_block_restarts_array);
avg_num_mat_vec = mean(num_mat_vec_array);

fprintf(' The average number of eigs block restarts was %1.3f\n', avg_num_block_restarts);
fprintf(' The average number of eigs A*x mults was %1.3f\n\n\n', avg_num_mat_vec)      

end



function [V, d, eig_solver_data] = solve_eigval_problem(Aty_fun, n, solver, opts)

   % Solves eigenvalue subproblem with user-selected solver
   if strcmp(solver, 'eigs')
      eigs_opts.isreal = opts.isreal;
      eigs_opts.k = opts.num_eigs;
      k = opts.num_eigs;
      if ~isempty(opts.V_prev)
         eigs_opts.v = opts.V_prev(:, 1); 
      end
      eigs_opts.p = opts.block_size;
      
      eigs_opts.tol = max(opts.eigs_tol_scale_tight*eps, opts.tol);      
      %{
      if opts.tol <= 1e3*eps
         eigs_opts.tol = opts.eigs_tol_scale*opts.tol;
      else
         eigs_opts.tol = opts.tol;
      end
      %}
      eigs_opts.maxit = opts.iter_max_solver;
      eigs_opts.issym  = true;
      eigs_opts.returnIterData = true;
      eigs_opts.disp = opts.disp;
      if opts.isreal
         sigma = 'LA';
      else
         sigma = 'LR';
      end
   
      try
         [V, D, ~, eigsIterData] = eigs_mod(Aty_fun, n, k, sigma, eigs_opts);
         
         eig_solver_data.eigsIterData = eigsIterData;
         eig_solver_data.eigs_opts = eigs_opts;
         
         d = diag(D);

         if isnan(d(1,1)) || isinf(d(1,1))
            Ax1 = Aty_fun(V(:,1));
            d(1,1) = reshape(V(:,1)'*Ax1, [1, 1]) / reshape(V(:,1)'*V(:,1), [1, 1]);
            if size(V, 2) >= 2
               Ax2 = Aty_fun(V(:,2));
               d(2,1) = reshape(V(:,2)'*Ax2, [1, 1]) / reshape(V(:,2)'*V(:,2), [1, 1]);
            end
         end
         
         %num_mat_vec = AClassObj.num_mat_vec_mults(y_iter);
         eig_solver_data.num_mat_vec = sum(eigsIterData.iter_vals(:, 2)) - 1;
         eig_solver_data.num_block_restarts = eigsIterData.iter_vals(end, 1);
         % Verifies eigenvalues are in descending order
         [d, d_idx] = sort(d, 'descend'); V = V(:, d_idx);

      catch ME
         d = nan;
         V = [];
         eig_solver_data.num_mat_vec = nan;
         eig_solver_data.num_block_restarts = nan;
      end

      
   elseif strcmp(solver, 'eigifp')
      k = opts.num_eigs;
      eigifp_opts.DISP = opts.disp;
      eigifp_opts.SIZE = n;
      
      % Preconditioning settings (requires A and B dense, not functions)
      %eigifp_opts.USEPRECON = 'NO';
      %eigifp_opts.USEPRECON = 1; % appx eigenvalue
      
      %eigifp_opts.ILUTHRESH = 0.5; % 0 = exact LU; 
                                    % (0,1) uses eigifp builtin function
                                    % 1 = 0-level ILU, using matlab `ilu`
                                    % labeled as `eta` in `eigifp`
      
      if ~isempty(opts.V_prev)
         eigifp_opts.V0 = opts.V_prev;
      end
      if ~isnan(opts.d_prev(1,1))
         eigifp_opts.NORMA = opts.d_prev(1,1);
      else
         eigifp_opts.NORMA = 1;
      end
      
      eigifp_opts.MAXIT = opts.iter_max_solver;
      eigifp_opts.INNERIT = opts.inner_iter;
      eigifp_opts.MAXINNERIT = opts.inner_iter_max;

      Aty_fun_neg = @(x) -Aty_fun(x);

      % scales tol to get appx same residual as other solvers
      eigifp_opts.TOL = max(opts.eigifp_tol_scale_loose*opts.tol, opts.eigifp_tol_scale_tight*eps);

      % guarantees tol is sufficiently small to get convergence
      eigifp_opts.TOL = min(eigifp_opts.TOL, opts.eigifp_tol_max);
      
           
      try
         [D_eigifp, V, eigifp_iter, num_mat_vec, ~] ...
             = eigifp_mod(Aty_fun_neg, k, eigifp_opts);
         d = -D_eigifp;

         if isnan(d(1,1)) || isinf(d(1,1))
            Ax1 = Aty_fun(V(:,1));
            d(1,1) = reshape(V(:,1)'*Ax1, [1, 1]) / reshape(V(:,1)'*V(:,1), [1, 1]);
            if size(V, 2) >= 2
               Ax2 = Aty_fun(V(:,2));
               d(2,1) = reshape(V(:,2)'*Ax2, [1, 1]) / reshape(V(:,2)'*V(:,2), [1, 1]);
            end
         end
         
         eig_solver_data.num_mat_vec = num_mat_vec;
         eig_solver_data.num_block_restarts = sum(eigifp_iter);
         eig_solver_data.eigifp_iter = eigifp_iter;
         % Verifies eigenvalues are in descending order
         [d, d_idx] = sort(d, 'descend'); V = V(:, d_idx);

      catch ME
         d = nan;
         V = [];
         eig_solver_data.num_mat_vec = nan;
         eig_solver_data.num_block_restarts = nan;
         eig_solver_data.eigifp_iter = nan;
      end


   elseif strcmp(solver, 'bleigifp')
      k = opts.num_eigs;
      bleigifp_opt.ADAPTTOL = 0.1; % default 0.1
%      bleigifp_opt.ADAPTTOL = eps; % forces block size to stay fixed
      if ~opts.bleigifp_use_default_block_method
         bleigifp_opt.UPDATEP = 'no'; % forces block size to stay fixed
         bleigifp_opt.UPDATEM = 'no'; % forces inner iter to stay fixed
         bleigifp_opt.INNERIT = opts.inner_iter;
         % BLOCKSIZE should be multiplicity or the cluster size of the eigenvalues desired per user manual.
         bleigifp_opt.BS = opts.block_size; % default 2
         bleigifp_opt.MAXBS = opts.max_block_size;         
      end      
%      bleigifp_opt.ARGSA = [];
      bleigifp_opt.DISP = opts.disp;
      bleigifp_opt.ILUTHRESH = 1e-04; % default:  1e-04
      if ~isempty(opts.V_prev)
         bleigifp_opt.INITIALVEC = opts.V_prev;
      end
      bleigifp_opt.MAXIT = opts.iter_max_solver;
      bleigifp_opt.NORMA = opts.d_prev(1,1);
      bleigifp_opt.NORMB = 1;
      bleigifp_opt.SIGMA = 'LA';
      bleigifp_opt.SIZE = n;
      
      bleigifp_opt.TOLERANCE = max(opts.bleigifp_tol_scale_loose*opts.tol, opts.bleigifp_tol_scale_tight*eps);
      %{
      if opts.tol <= 1e3*eps
         bleigifp_opt.TOLERANCE = opts.bleigifp_tol_scale*eps;
      else
         bleigifp_opt.TOLERANCE = 1e-1*opts.tol;
      end
      %}
      
%      bleigifp_opt.UPDATEM = 'no'; % default is 'yes'
%      bleigifp_opt.UPDATEP = 'no'; % default is 'yes'
%      bleigifp_opt.USEPRECON = opts.d_prev(1,1); % only valid for numeric A
      
      try
         [d, V, resHist, cummatvec] = bleigifp(Aty_fun, k, bleigifp_opt);
         if isnan(d(1,1)) || isinf(d(1,1))
            Ax1 = Aty_fun(V(:,1));
            d(1,1) = reshape(V(:,1)'*Ax1, [1, 1]) / reshape(V(:,1)'*V(:,1), [1, 1]);
            if size(V, 2) >= 2
               Ax2 = Aty_fun(V(:,2));
               d(2,1) = reshape(V(:,2)'*Ax2, [1, 1]) / reshape(V(:,2)'*V(:,2), [1, 1]);
            end
         end
         eig_solver_data.num_mat_vec = cummatvec(end);
         eig_solver_data.num_block_restarts = length(cummatvec);
         % Verifies eigenvalues are in descending order
         [d, d_idx] = sort(d, 'descend'); V = V(:, d_idx);

      catch ME
         d = nan;
         V = [];
         eig_solver_data.num_mat_vec = nan;
         eig_solver_data.num_block_restarts = nan;
      end
      
      
      
   elseif strcmp(solver, 'LOPCG')
      Aty_fun_neg = @(x) -Aty_fun(x);
      LOPCG_opts.display = opts.disp;
      nrmAB = [opts.d_prev(1, 1); 1];
      shift = 0.0;       
      if opts.tol <= eps
         LOPCG_tol1 = 1e2*eps;
         LOPCG_tol2 = LOPCG_tol1;
      else
         LOPCG_tol1 = 1e0*opts.tol;
         LOPCG_tol2 = LOPCG_tol1;
      end

      LOPCG_opts.precond = 1;
      LOPCG_opts.precond_one = 2;
      LOPCG_opts.precond_sgm = -0.001 - opts.d_prev(1,1); %-(0.001 + (D_prev(1,1) - D_prev(2,1))); %-0.1 + D(1,1);
      
      if ~isempty(opts.V_prev)
         V_prev1 = opts.V_prev(:, 1); V_prev2 = opts.V_prev(:, 2);
      else
         V_prev1 = rand(n,1); V_prev1 = V_prev1/norm(V_prev1);
         V_prev2 = rand(n,1); V_prev2 = V_prev2/norm(V_prev2);
      end
      [lam1, V_LOPCG, res, hsty, LOPCG_info1] ...
        = LOPCGgS(Aty_fun, n, sparse(eye(n,n)), V_prev1, [], nrmAB, ...
                  shift, LOPCG_tol1, opts.iter_max_solver, LOPCG_opts);

      % Call solver 2nd time to get 2nd largest eigenvalue
      LOPCG_opts.precond_sgm = -0.001 - opts.d_prev(2,1); %-(0.001 + (D_prev(1,1) - D_prev(2,1)));
      shift = -0.01-lam1;

      [lam2, V_LOPCG2, res, hsty, LOPCG_info2] ...
        = LOPCGgS(Aty_fun, n, sparse(eye(n,n)), V_prev2, V_LOPCG, nrmAB, ...
                  shift, LOPCG_tol2, opts.iter_max_solver, LOPCG_opts);        

      V = [V_LOPCG, V_LOPCG2];
      d = -[lam1; lam2];
      
      if isnan(d(1,1)) || isinf(d(1,1))
         Ax1 = Aty_fun(V(:,1));
         d(1,1) = reshape(V(:,1)'*Ax1, [1, 1]) / reshape(V(:,1)'*V(:,1), [1, 1]);
         if size(V, 2) >= 2
            Ax2 = Aty_fun(V(:,2));
            d(2,1) = reshape(V(:,2)'*Ax2, [1, 1]) / reshape(V(:,2)'*V(:,2), [1, 1]);
         end
      end
      % Verifies eigenvalues are in descending order
      [d, d_idx] = sort(d, 'descend'); V = V(:, d_idx);
    
      eig_solver_data.num_mat_vec = LOPCG_info1.num_mat_vec_mults + LOPCG_info2.num_mat_vec_mults;
      eig_solver_data.num_block_restarts = nan;
   end
   
end



function opts = get_eigs_param_update(opts, data, is_init_iter, y_iter)

% Initial iterate
if is_init_iter || (y_iter==1)
   opts.num_eigs_prev = opts.num_eigs;
   opts.num_eigs = opts.num_eigs + 1;   
else
   shift_dir = opts.num_eigs - opts.num_eigs_prev;
   num_mvs = data.num_mat_vec{y_iter};
   num_mvs_prev = data.num_mat_vec{y_iter-1};
   opts.num_eigs_prev = opts.num_eigs;
   % Shifts away from min and max num eigs
   if opts.num_eigs == opts.num_eigs_min
      opts.num_eigs = opts.num_eigs + 1;
   elseif opts.num_eigs == opts.num_eigs_max
      opts.num_eigs = opts.num_eigs - 1;
   % Shift in direction of fewer mat-vecs
   elseif num_mvs < num_mvs_prev
      opts.num_eigs = opts.num_eigs + shift_dir;
   else
      opts.num_eigs = opts.num_eigs - shift_dir;
   end   
end
   
end



function [opts] = get_opts(ip)
opts.num_eigs = ip.Results.num_eigs;
opts.block_size = ip.Results.block_size;
opts.max_block_size = ip.Results.max_block_size;
opts.inner_iter = ip.Results.inner_iter;
opts.inner_iter_max = ip.Results.inner_iter_max;
opts.iter_max_solver = ip.Results.iter_max_solver;
opts.eigifp_use_default_block_method = ip.Results.eigifp_use_default_block_method;
opts.bleigifp_use_default_block_method = ip.Results.bleigifp_use_default_block_method;
opts.eigs_use_adaptive_method = ip.Results.eigs_use_adaptive_method;
opts.num_eigs_max = ip.Results.num_eigs_max;
opts.num_eigs_min = ip.Results.num_eigs_min;
opts.eigifp_tol_max = ip.Results.eigifp_tol_max;
opts.eigs_tol_scale_tight = ip.Results.eigs_tol_scale_tight;
opts.eigifp_tol_scale_tight = ip.Results.eigifp_tol_scale_tight;
opts.bleigifp_tol_scale_tight = ip.Results.bleigifp_tol_scale_tight;
opts.eigifp_tol_scale_loose = ip.Results.eigifp_tol_scale_loose;
opts.bleigifp_tol_scale_loose = ip.Results.bleigifp_tol_scale_loose;
opts.disp = ip.Results.disp;
end

