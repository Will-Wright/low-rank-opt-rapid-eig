function [results, data] = solve_noisyimage_saga_vs_exact_SDP_solutions(varargin)
% experiments.exputil.solve_noisyimage_saga_vs_exact_SDP_solutions
%
% Solves noisy phase retrieval problem with saga, then solves primal and 
% dual and gauge dual exactly using `cvx`

ip = inputParser;
ip.addParameter('n', 16);
ip.addParameter('L', 6);
ip.addParameter('y_scale', 1);
ip.addParameter('noise_ratio', 0.20);
ip.addParameter('iterations', 2000);
ip.addParameter('noise_type', 'gaussian');
ip.addParameter('signal', 'gaussian');
ip.addParameter('feaTol', 1e-5);
ip.addParameter('optTol', 1e-5);
ip.addParameter('primalRelErrTol', 1e-6);
ip.addParameter('dualRelErrTol', 1e-4);
ip.addParameter('StopOnRelErr', true);
ip.addParameter('verbosity', false);
ip.addParameter('seed', 0);
ip.addParameter('rank_tol', 1e-6);
ip.parse(varargin{:});

n = ip.Results.n;
L = ip.Results.L;
m = n*L;
y_scale = ip.Results.y_scale;
noise_ratio = ip.Results.noise_ratio;
iterations = ip.Results.iterations;
noise_type = ip.Results.noise_type;
signal = ip.Results.signal;
feaTol = ip.Results.feaTol;
optTol = ip.Results.optTol;
primalRelErrTol = ip.Results.primalRelErrTol;
dualRelErrTol = ip.Results.dualRelErrTol;
StopOnRelErr = ip.Results.StopOnRelErr;
verbosity = ip.Results.verbosity;
seed = ip.Results.seed;
rank_tol = ip.Results.rank_tol;

fprintf('Solving gauge dual phase retrieval problem with saga: ');
pl_opts = struct('n', n, 'L', L, 'noise_ratio', noise_ratio, 'scale', y_scale, ...
   'signal', signal, 'noise_type', noise_type, ...
   'StopOnRelErr', StopOnRelErr, 'iterations', iterations, ...
   'feaTol', feaTol, 'optTol', optTol, ...
   'primalRelErrTol', primalRelErrTol, 'dualRelErrTol', dualRelErrTol, ...
   'verbosity', verbosity, 'seed', seed);
results = experiments.evolvingmatricespl(pl_opts);
fprintf('Solved!\n');
x_saga = vec(results.saga_sd.x);
y_saga = vec(results.saga_sd.solInfo.y);
A = results.gendatapl.A;
b = results.gendatapl.b;
b_vec = vec(b);
b0 = vec(results.gendatapl.genInfo.b0);
x0 = vec(results.gendatapl.x0);
%y0 = vec(results.gendatapl.genInfo.y0);
epsilon = noise_ratio * norm(b_vec);


% Generates matrices A_i for linear matrix expression sum y_i * A*i

A_array = zeros(n, n, m);
for i = 1:m
   ei = zeros(m, 1);
   ei(i, 1) = 1;
   A_array(:, :, i) = experiments.exputil.getfulldualmatrix(results.gendatapl.A, ei);
end
A_mat = reshape(A_array, [n*n, m]);


%{
% This test verifies the transformation of gauge dual phase lift
% problem to cvx is correct
A_test(:,:,1) = [1 2; 3 4]';
A_test(:,:,2) = [5 6; 7 8]';
A_test(:,:,3) = [9 10; 11 12]';
y = [2; 3; 4];
A_test
y
sum(A_test.*reshape(y, [1, 1, 3]), 3)
% USE THIS
Amat_test = reshape(A_test, [2*2, 3]);
reshape(Amat_test*y, [2, 2])
%}


% Solves primal phase retrieval problem with cvx

fprintf('Solving primal phase retrieval problem using cvx: ');
if verbosity
   cvx_begin
else
   cvx_begin quiet
end
   variable X_cvx_primal(n, n) hermitian
   minimize norm_nuc(X_cvx_primal)
   subject to
      norm(A_mat'*vec(X_cvx_primal) - b_vec) <= epsilon;
      X_cvx_primal == hermitian_semidefinite(n);
cvx_end
fprintf('Solved!\n');


% Solves gauge dual with cvx

fprintf('Solving gauge dual phase retrieval problem using cvx: ');
if verbosity
   cvx_begin
else
   cvx_begin quiet
end
   variable y_cvx_dual(m, 1);
   minimize lambda_max(reshape(A_mat*y_cvx_dual, [n, n]));
   subject to
      b_vec'*y_cvx_dual - epsilon*norm(y_cvx_dual, 2) >= 1;
cvx_end
fprintf('Solved!\n');


% Generates primal solution using cvx dual solution

Aty_cvx_dual = experiments.exputil.getfulldualmatrix(A, y_cvx_dual);
Aty_cvx_dual = experiments.exputil.getfulldualmatrix(A, y_cvx_dual);
[V_cvx_dual, D_cvx_dual] = eig(Aty_cvx_dual);
[lam1_cvx_dual, lam1_cvx_idx] = max(real(diag(D_cvx_dual)));
D_cvx_dual_deflated = real(diag(D_cvx_dual));
D_cvx_dual_deflated(lam1_cvx_idx, lam1_cvx_idx) = 0;
lam2_cvx_dual = max(real(D_cvx_dual_deflated));
v_cvx_dual = V_cvx_dual(:, lam1_cvx_idx);
Av = vec(A.forward(v_cvx_dual));
b_eps_cvx = b_vec - epsilon*y_cvx_dual/norm(y_cvx_dual);
s = dot(Av, b_eps_cvx)/ (norm(Av))^2;

%x_cvx_dual = sqrt(s)*v_cvx_dual;
fprintf('Recovering primal solution from dual cvx solution using wflow:');
[~, g, v] = A.gdobjective(reshape(y_cvx_dual, [A.m1, A.m2, A.m3]), 1);
if verbosity
   fprintf('\n');
end
[x_cvx_dual, stats_cvx] = pfd(A, b, epsilon, reshape(y_cvx_dual, [A.m1, A.m2, A.m3]), v, g, 'verbose', verbosity);
%[x_cvx_dual, stats_cvx] = pfd(A, b, 0, reshape(y_cvx_dual, [A.m1, A.m2, A.m3]), v, g, 'verbose', true);
fprintf('Solved!\n');

%{
minFunc_opts.Method = 'sd';
minFunc_opts.MAXFUNEVALS = 10;
minFunc_opts.optTol = 1e-8;
minFunc_opts.progTol = 1e-8;
minFunc_opts.MaxIter = 300;
minFunc_opts.MaxFunEvals = inf;
minFunc_opts.Display = 'iter';

y0 = vec(results.saga_sd.solInfo.y);
objFun = @(y) duObjFun(y, A, b, noise_ratio, y_scale);

[y_true,fval,exitflag,output] = minFunc(objFun, y0, minFunc_opts);
fval
%}


% Measures cvx and saga_sd residuals and gauge dual optimality condition violations

%[U_cvx_primal, S_cvx_primal, V_cvx_primal] = svd(X_cvx_primal);
%S_diag_cvx_primal = sort(diag(S_cvx_primal));

[V_cvx_primal, d_cvx_primal] = eig(0.5*(X_cvx_primal + X_cvx_primal'));
S_diag_cvx_primal = sort(diag(d_cvx_primal));

rank_cvx_primal = sum(S_diag_cvx_primal >= rank_tol);
[sing_val, sing_val_idx] = max(diag(d_cvx_primal));
x_cvx_primal = sqrt(sing_val)*V_cvx_primal(:, sing_val_idx);

%[U_cvx_dual, S_cvx_dual, V_cvx_dual] = svd(Aty_cvx_dual);
%S_diag_cvx_dual = sort(diag(S_cvx_dual));

[V_cvx_dual, d_eig] = eig(0.5*(Aty_cvx_dual + Aty_cvx_dual'));
S_diag_cvx_dual = sort(diag(d_eig));

spectral_opt_cvx = sum( abs( S_diag_cvx_primal.*(lam1_cvx_dual*ones(n,1) - S_diag_cvx_dual ) ) );

spectral_opt_saga = 0; % value is zero since saga enforces rank-1 solution

VtV = V_cvx_primal'*V_cvx_dual;
[~, VtV_idx] = max(abs(VtV));
VtV = VtV(VtV_idx, :);
VtV = diag(sign(real(diag(VtV))))*VtV;

x0_fro_norm = util.hermitianerror(zeros(size(x0)), x0, 'fro');
Ax_cvx_dual = vec(A.forward(x_cvx_dual)); Ax_saga = vec(A.forward(x_saga));
Ax_cvx_primal = vec(A.forward(x_cvx_primal)); AX_cvx_primal = A_mat'*vec(X_cvx_primal);

rank_X_cvx_primal = rank(X_cvx_primal, rank_tol);

pr_obj_cvx_primal = sum(S_diag_cvx_primal);
pr_obj_cvx_dual = norm(x_cvx_dual(:))^2; pr_obj_saga = norm(x_saga(:))^2;
du_obj_cvx_dual = A.gdobjective(y_cvx_dual, 1); du_obj_saga = A.gdobjective(y_saga,1);
fprintf('\n');
if verbosity
   fprintf('Below are the optimality condition violations for the cvx primal/dual pair (X,y) and the saga pair (xx'', y):\n\n');
end
fprintf('Note: cvx obtained optimal %i-by-%i solution matrix X with rank %i.\n\n', n, n, rank_cvx_primal);
fprintf('              complementarity condition                               | cvx-primal/dual pair |   saga_sd\n');
fprintf('-----------------------------------------------------------------------------------------------------------\n');
fprintf(' 1)         primal opt: abs(||A(X) - b|| - epsilon)                   |      %1.4e      | %1.4e\n', ...
   abs(norm(AX_cvx_primal - b_vec) - epsilon), abs(norm(Ax_saga - b_vec) - epsilon));
fprintf(' 2)           dual opt:  abs( b''*y - epsilon*||y|| - 1 )              |      %1.4e      | %1.4e\n', ...
   abs(dot(y_cvx_dual, b_vec) - epsilon*norm(y_cvx_dual) - 1), abs(dot(y_saga, b_vec) - epsilon*norm(y_saga) - 1));
fprintf(' 3) Cauchy-Schwarz opt: abs( (b-A(X))''*y - ||y||*||A(X)-b||           |      %1.4e      | %1.4e\n', ...
   abs( (AX_cvx_primal - b_vec)'*y_cvx_dual + norm(y_cvx_dual)*norm(AX_cvx_primal - b_vec)), ...
   abs( (Ax_saga - b_vec)'*y_saga + norm(y_saga)*norm(Ax_saga - b_vec)));
fprintf(' 4)       spectral opt: sum(abs(lam_i(X)*(lam_1(A''y) - lam_i(A''y))))  |      %1.4e      | %1.4e\n', ...
   spectral_opt_cvx, spectral_opt_saga);
%fprintf(' 5) (X, A''y) simul svd: (1/n^2)*||V_X''*V_y - I||                     |      %1.4e      | %1.4e\n', ...
%   norm(UtU - eye(n,n))/(n^2), 1);
fprintf(' 5)     strong duality: abs(1 - prObj*duObj)                          |      %1.4e      | %1.4e\n', ...
   abs(1-pr_obj_cvx_primal*du_obj_cvx_dual/y_scale), abs(1-pr_obj_saga*du_obj_saga/y_scale));
fprintf('\n');
if verbosity
   fprintf('Next we measure the solution residuals using the true signal x0.\n');
   fprintf('To obtain the approximate signal from the cvx primal and dual solutions, we do the following:\n\n');
   fprintf('  x_cvx_primal: select the leading eigenvector of X_cvx\n');
   fprintf('  x_cvx_dual: apply gradient descent (i.e., wflow) with y_cvx_dual as initial input\n\n');
end
fprintf('Note: rank of X_cvx_primal is %i\n', rank_X_cvx_primal);
fprintf('                  residual value                    | cvx-primal |  cvx-dual  |   saga_sd\n');
fprintf('--------------------------------------------------------------------------------------------\n');
fprintf('  primal rel error: ||xx''-x0x0''||_2 / ||x0x0''||_2   | %1.4e | %1.4e | %1.4e\n', ... 
   util.hermitianerror(x_cvx_primal, x0, 'fro')/x0_fro_norm, util.hermitianerror(x_cvx_dual, x0, 'fro')/x0_fro_norm, ...
   util.hermitianerror(x_saga, x0, 'fro')/x0_fro_norm);
fprintf('primal feasibility: max(0,||A(X) - b|| - epsilon)   | %1.4e | %1.4e | %1.4e\n', ...
   max(0, norm(vec(AX_cvx_primal) - b_vec) - epsilon), max(0, norm(Ax_cvx_dual - b_vec) - epsilon), max(0, norm(Ax_saga - b_vec) - epsilon));
fprintf('measured obs error: ||A(X) - b_obs||                | %1.4e | %1.4e | %1.4e\n', ...
   norm(vec(AX_cvx_primal) - b_vec), norm(Ax_cvx_dual(:) - b_vec), norm(Ax_saga - b_vec));
fprintf('    true obs error: ||A(X) - b_true||               | %1.4e | %1.4e | %1.4e\n', ...
   norm(vec(AX_cvx_primal) - b0), norm(Ax_cvx_dual(:) - b0), norm(Ax_saga - b0));
%fprintf('||y_cvx - y0|| = %1.4e\n', norm(y_cvx - y0));
%fprintf('||y_saga - y0|| = %1.4e\n', norm(y_saga - y0));

data.rank_X_cvx_primal = rank_X_cvx_primal;
data.rank_cvx_primal = rank_cvx_primal;
data.A_mat = A_mat;
data.Aty_cvx_dual = Aty_cvx_dual;
data.X_cvx_primal = X_cvx_primal;
data.x_cvx_primal = x_cvx_primal;
data.x_cvx_dual = x_cvx_dual;
data.y_cvx_dual = y_cvx_dual;
data.x_saga = results.saga_sd.x;
data.y_saga = y_saga;

end


function [f, g] = duObjFun(y_vec, A, b, noise_ratio, y_scale)
   y = reshape(y_vec, [A.m1, A.m2, A.m3]);
   y = project(y, b, noise_ratio, y_scale);
   Aty = experiments.exputil.getfulldualmatrix(A, y);
   [V_eig, D_eig] = eig(Aty);
   [f_eig, f_idx] = max(real(diag(D_eig)));
   g_eig = vec(A.forward(V_eig(:, f_idx)));

%   verifies eig computes correct eigenvalue
%   [f_eigs, g_eigs] = A.gdobjective(y, 1);
%    g_eigs = vec(g_eigs);
%   [f_eig, f_eigs, f_eig-f_eigs]
%   norm(g_eig-g_eigs)
    
   f = f_eig;
   g = g_eig;
end


