function [results] = solve_EMEPs_various_solvers_and_params(varargin)
% function [results] = solve_EMEPs_various_solvers_and_params(varargin)
%
% Solves EMEP with `eigs`, `eigifp`, and `bleigifp` using various solver
% parameters to generate solver performance information for a set of phase
% retrieval problems with various noise ratios and oversampling rates.
%
% Note that this routine can take a very long time if the data has not
% already been cached.

ip = inputParser;
ip.addParameter('signal', 'image'); % can be 'image' or 'gaussian'
ip.addParameter('n', 4096);
ip.addParameter('solver', 'eigs');
ip.addParameter('image_file', 'data/parrot_4k.jpg');
ip.addParameter('disp', 0);
ip.addParameter('folder_root', 'cache/proto/exp_EMEP_various_params/');
ip.parse(varargin{:});

signal = ip.Results.signal;
n = ip.Results.n;
solver = ip.Results.solver;
image_file = ip.Results.image_file;
folder_root = ip.Results.folder_root;
disp = ip.Results.disp;

if strcmp(signal, 'gaussian')
   folder_exp = 'gaussian/';
elseif strcmp(image_file, 'data/parrot_4k.jpg')
   folder_exp = 'parrot_4k/';
elseif strcmp(image_file, 'data/parrot_16k.jpg')
   folder_exp = 'parrot_16k/';
else
   fprintf('Invalid image file provided.');
   return
end

if strcmp(solver, 'eigs')
   file_test = 'eigs.mat';
elseif strcmp(solver, 'eigifp')
   file_test = 'eigifp.mat';
elseif strcmp(solver, 'bleigifp')
   file_test = 'bleigifp.mat';
else
   fprintf('Invalid solver provided.');
   return
end


if strcmp(signal, 'gaussian')
   % Phase model parameters
   L_range = [5, 10];
   noise_ratio_range = [0.05, 0.15, 0.30];

   % EMEP parameters
   EMEP_opts.eigs.block_size_range = 20:20:100;
   EMEP_opts.eigs.num_eigs_range = 2:1:30;

   EMEP_opts.eigifp.inner_iter_range = 1:30;
   EMEP_opts.eigifp.num_eigs_range = 2;

   EMEP_opts.bleigifp.block_size_range = 2:1:5;
   EMEP_opts.bleigifp.inner_iter_range = 4:4:40;
   EMEP_opts.bleigifp.num_eigs_range = 2;
elseif strcmp(image_file, 'data/parrot_4k.jpg')
   % Phase model parameters
   L_range = [5, 10];
   noise_ratio_range = [0.05, 0.15, 0.30];

   % EMEP parameters
   EMEP_opts.eigs.block_size_range = 20:20:100;
   EMEP_opts.eigs.num_eigs_range = 2:1:30;

   EMEP_opts.eigifp.inner_iter_range = 1:30; % fails at 17 for exp with L = 10, noise_ratio = 0.3
   EMEP_opts.eigifp.num_eigs_range = 2;

   EMEP_opts.bleigifp.block_size_range = 2:1:5;
   EMEP_opts.bleigifp.inner_iter_range = 4:4:40;
   EMEP_opts.bleigifp.num_eigs_range = 2;
elseif strcmp(image_file, 'data/parrot_16k.jpg')
   % Phase model parameters
   L_range = [5, 10];
   noise_ratio_range = [0.05, 0.15, 0.30];

   % EMEP parameters
   EMEP_opts.eigs.block_size_range = 20:20:60;
   EMEP_opts.eigs.num_eigs_range = 2:1:15;

   EMEP_opts.eigifp.inner_iter_range = 2:2:20;
   EMEP_opts.eigifp.num_eigs_range = 2;

   EMEP_opts.bleigifp.block_size_range = 2:1:5;
   EMEP_opts.bleigifp.inner_iter_range = 4:4:40;
   EMEP_opts.bleigifp.num_eigs_range = 2;
else
   fprintf('Invalid signal or image file parameters.\n');
   return;
end


% Search for saved results
if exist(folder_root) ~= 7
   mkdir(folder_root)
end
if exist(strcat(folder_root, folder_exp)) ~= 7
   mkdir(strcat(folder_root, folder_exp));
end

if exist(strcat(folder_root, folder_exp, 'results.mat')) ~= 2 % folder not found
   fprintf('Experiment did not identify a set of solved phase problems\n');
   fprintf('Solving 6 phase problems\n');
   results = cell(0,0);
   for L_idx = 1:2
      for noise_ratio_idx = 1:3
         L = L_range(L_idx);
         noise_ratio = noise_ratio_range(noise_ratio_idx);
         results{L_idx, noise_ratio_idx} = experiments.evolvingmatricespl(...
            'signal', signal, 'n', n, 'L', L, 'noise_ratio', noise_ratio, ...
            'image_file', image_file, 'eigIts', 30000, ...
            'cutDim', 20, 'eigs_basis_size', 50);
         
         % Computes largest magnitude eigenvalues for all experiments
         results_temp = results{L_idx, noise_ratio_idx};

         A_data = experiments.exputil.convert_results_to_A_matrix_file(results_temp);
         AClassObj = experiments.exputil.A_class(A_data);

         for iter = 1:results_temp.saga_sd.solInfo.nit

            Aty_fun = @(x) AClassObj.A_fun(iter, x);
            [~, d] = eigs(Aty_fun, A_data.n, 1, 'LM');
            results{L_idx, noise_ratio_idx}.saga_sd.solInfo.iterData.largest_mag_eigval{iter} = d;

            % Verifies largest magnitude eigenvalue comes from same matrix as largest
            % real eigenvalue
            %[~, d] = eigs(Aty_fun, A_data.n, 1, 'LR');
            %fprintf('%1.3e\n', d - results_EMEP_solvers.eigs{L_idx,noise_ratio_idx,1,1}.d{iter}(1,1));
         end
         
      end
   end
   save(strcat(folder_root, folder_exp, 'results.mat'), 'results');
else
   load(strcat(folder_root, folder_exp, 'results.mat'), 'results')
end   

if exist(strcat(folder_root, folder_exp, file_test)) ~= 2 % no test file found
   results_EMEP_solvers = struct;
   results_EMEP_solvers.exp_is_complete = false;
   save(strcat(folder_root, folder_exp, file_test), 'results_EMEP_solvers');
else
   load(strcat(folder_root, folder_exp, file_test), 'results_EMEP_solvers')
end 

results_EMEP_solvers = solve_EMEPs(results, EMEP_opts, L_range, noise_ratio_range, ...
   solver, folder_root, folder_exp, file_test, results_EMEP_solvers, disp);

results_EMEP_solvers.exp_is_complete = true;
save(strcat(folder_root, folder_exp, file_test), 'results_EMEP_solvers')

% Deletes unnecessary cached data
delete(strcat(folder_root, folder_exp, 'results.mat'))

fprintf('\nAll experiments have been solved for signal type ''%s'' with image ''%s'' using solver ''%s''.\n', signal, image_file, solver);

end








function [results_EMEP_solvers] = solve_EMEPs(results, EMEP_opts, L_range, noise_ratio_range, ...
   solver, folder_root, folder_exp, file_test, results_EMEP_solvers, disp)
   iter_max_solver = 5000;  % this is modified dynamically below for 'eigifp'
   mat_vec_max = 5000; % used for 'eigifp' and 'bleigifp'
   iter_start = 1;

   
% Runs eigs tests
if strcmp(solver, 'eigs')
   if ~isfield(results_EMEP_solvers, 'eigs')
      results_EMEP_solvers.eigs = cell(length(L_range), ...
         length(noise_ratio_range), length(EMEP_opts.eigs.block_size_range), ...
         length(EMEP_opts.eigs.num_eigs_range));
   end
   inner_iter = nan; % unused in eigs
   for L_idx = length(L_range):-1:1
      for noise_ratio_idx = length(noise_ratio_range):-1:1
         iter_end = results{L_idx, noise_ratio_idx}.saga_sd.solInfo.nit;
         for block_size_idx = 1:length(EMEP_opts.eigs.block_size_range)
            block_size = EMEP_opts.eigs.block_size_range(block_size_idx);
            for num_eigs_idx = 1:length(EMEP_opts.eigs.num_eigs_range)
               num_eigs = EMEP_opts.eigs.num_eigs_range(num_eigs_idx);
               if num_eigs <= (block_size - 5)
                  
                  if isempty(results_EMEP_solvers.eigs{L_idx, noise_ratio_idx, block_size_idx, num_eigs_idx})

                     exp_params = struct('num_eigs', num_eigs, 'block_size', block_size, ...
                        'inner_iter', inner_iter, 'iter_max_solver', iter_max_solver, ...
                        'iter_start', iter_start, 'iter_end', iter_end, 'disp', disp, ...
                        'solver', solver);
                     results_EMEP_solvers.eigs{L_idx, noise_ratio_idx, block_size_idx, num_eigs_idx} ...
                        = experiments.exputil.solve_EMEP(...
                           'results', results{L_idx, noise_ratio_idx}, exp_params);
                     results_EMEP_solvers.eigs{L_idx, noise_ratio_idx, block_size_idx, num_eigs_idx}.exp_params ...
                        = exp_params;
                     results_EMEP_solvers.eigs{L_idx, noise_ratio_idx, block_size_idx, num_eigs_idx}.exp_params.L = L_range(L_idx);
                     results_EMEP_solvers.eigs{L_idx, noise_ratio_idx, block_size_idx, num_eigs_idx}.exp_params.noise_ratio = noise_ratio_range(noise_ratio_idx);
               
                     save(strcat(folder_root, folder_exp, file_test), 'results_EMEP_solvers');                             
                     
                  end
               end
            end
         end
      end
   end
   L_idx = []; noise_ratio_idx = []; inner_iter_idx = []; block_size_idx = []; num_eigs_idx = [];


% Runs eigifp tests
elseif strcmp(solver, 'eigifp')
   % Solves EMEP with default eigifp parameters
   if ~isfield(results_EMEP_solvers, 'eigifp_default')
      results_EMEP_solvers.eigifp_default = cell(length(L_range), ...
         length(noise_ratio_range));
   end
   for L_idx = length(L_range):-1:1
      for noise_ratio_idx = length(noise_ratio_range):-1:1
         iter_end = results{L_idx, noise_ratio_idx}.saga_sd.solInfo.nit;
         if isempty(results_EMEP_solvers.eigifp_default{L_idx, noise_ratio_idx})
            num_eigs = 2;

            % iter max changed to enforce same number of mat-vecs per exp
            iter_max_solver = round(mat_vec_max / 2, -2);
            inner_iter_max = 32;

            exp_params = struct('num_eigs', num_eigs, ...
               'eigifp_use_default_block_method', true, ...
               'iter_max_solver', iter_max_solver, ...
               'iter_start', iter_start, 'iter_end', iter_end, 'disp', disp, ...
               'solver', solver, 'inner_iter_max', inner_iter_max); 
            fprintf('Solving EMEP with eigifp, L = %i, noise ratio = %1.3f, with default eigifp settings\n', ...
               L_range(L_idx), noise_ratio_range(noise_ratio_idx));                  
            results_EMEP_solvers.eigifp_default{L_idx, noise_ratio_idx} ...
               = experiments.exputil.solve_EMEP(...
                  'results', results{L_idx, noise_ratio_idx}, exp_params);
            results_EMEP_solvers.eigifp_default{L_idx, noise_ratio_idx}.exp_params ...
               = exp_params;
            results_EMEP_solvers.eigifp_default{L_idx, noise_ratio_idx}.exp_params.L = L_range(L_idx);
            results_EMEP_solvers.eigifp_default{L_idx, noise_ratio_idx}.exp_params.noise_ratio = noise_ratio_range(noise_ratio_idx);


            save(strcat(folder_root, folder_exp, file_test), 'results_EMEP_solvers');                             
         end
      end
   end
   
   
   % Solves EMEP with various eigifp parameters
   if ~isfield(results_EMEP_solvers, 'eigifp')
      results_EMEP_solvers.eigifp = cell(length(L_range), ...
         length(noise_ratio_range), length(EMEP_opts.eigifp.inner_iter_range), ...
         length(EMEP_opts.eigifp.num_eigs_range));
   end
   block_size = nan; % unused in eigifp
   for L_idx = 1:length(L_range)
      for noise_ratio_idx = length(noise_ratio_range):-1:1
         iter_end = results{L_idx, noise_ratio_idx}.saga_sd.solInfo.nit;
         for inner_iter_idx = 1:length(EMEP_opts.eigifp.inner_iter_range)
            inner_iter = EMEP_opts.eigifp.inner_iter_range(inner_iter_idx);
            for num_eigs_idx = 1:length(EMEP_opts.eigifp.num_eigs_range)
               if isempty(results_EMEP_solvers.eigifp{L_idx, noise_ratio_idx, inner_iter_idx, num_eigs_idx})
                  num_eigs = EMEP_opts.eigifp.num_eigs_range(num_eigs_idx);
                  
                  % iter max changed to enforce same number of mat-vecs per exp
                  iter_max_solver = round(mat_vec_max / inner_iter, -2);
                  
                  exp_params = struct('num_eigs', num_eigs, 'block_size', block_size, ...
                     'inner_iter', inner_iter, 'iter_max_solver', iter_max_solver, ...
                     'iter_start', iter_start, 'iter_end', iter_end, 'disp', disp, ...
                     'solver', solver); 
                  fprintf('Solving EMEP with eigifp, L = %i, noise ratio = %1.3f, num_eigs = %i, inner_block = %i\n', ...
                     L_range(L_idx), noise_ratio_range(noise_ratio_idx), num_eigs, inner_iter);                  
                  results_EMEP_solvers.eigifp{L_idx, noise_ratio_idx, inner_iter_idx, num_eigs_idx} ...
                     = experiments.exputil.solve_EMEP(...
                        'results', results{L_idx, noise_ratio_idx}, exp_params);
                  results_EMEP_solvers.eigifp{L_idx, noise_ratio_idx, inner_iter_idx, num_eigs_idx}.exp_params ...
                     = exp_params;
                  results_EMEP_solvers.eigifp{L_idx, noise_ratio_idx, inner_iter_idx, num_eigs_idx}.exp_params.L = L_range(L_idx);
                  results_EMEP_solvers.eigifp{L_idx, noise_ratio_idx, inner_iter_idx, num_eigs_idx}.exp_params.noise_ratio = noise_ratio_range(noise_ratio_idx);
                  
                  
                  save(strcat(folder_root, folder_exp, file_test), 'results_EMEP_solvers');                             

               end
            end
         end
      end
   end
   L_idx = []; noise_ratio_idx = []; inner_iter_idx = []; block_size_idx = []; num_eigs_idx = [];


% Runs bleigifp tests
elseif strcmp(solver, 'bleigifp')
   % Solves EMEP with default bleigifp parameters
   if ~isfield(results_EMEP_solvers, 'bleigifp_default')
      results_EMEP_solvers.bleigifp_default = cell(length(L_range), ...
         length(noise_ratio_range));
   end
   for L_idx = length(L_range):-1:1
      for noise_ratio_idx = length(noise_ratio_range):-1:1
         iter_end = results{L_idx, noise_ratio_idx}.saga_sd.solInfo.nit;
         if isempty(results_EMEP_solvers.bleigifp_default{L_idx, noise_ratio_idx})
            num_eigs = 2;

            % iter max changed to enforce same number of mat-vecs per exp
            iter_max_solver = round(mat_vec_max / (1) );


            exp_params = struct('num_eigs', num_eigs, 'iter_max_solver', iter_max_solver, ...
               'iter_start', iter_start, 'iter_end', iter_end, 'disp', disp, ...
               'solver', solver, 'bleigifp_use_default_block_method', true); 
            fprintf('Solving EMEP with bleigifp, L = %i, noise ratio = %1.3f, with default bleigifp settings\n', ...
               L_range(L_idx), noise_ratio_range(noise_ratio_idx));                  
            results_EMEP_solvers.bleigifp_default{L_idx, noise_ratio_idx} ...
               = experiments.exputil.solve_EMEP(...
                  'results', results{L_idx, noise_ratio_idx}, exp_params);
            results_EMEP_solvers.bleigifp_default{L_idx, noise_ratio_idx}.exp_params ...
               = exp_params;
            results_EMEP_solvers.bleigifp_default{L_idx, noise_ratio_idx}.exp_params.L = L_range(L_idx);
            results_EMEP_solvers.bleigifp_default{L_idx, noise_ratio_idx}.exp_params.noise_ratio = noise_ratio_range(noise_ratio_idx);


            save(strcat(folder_root, folder_exp, file_test), 'results_EMEP_solvers');                             
         end
      end
   end
   
   
   
   
   
   % Solves EMEP with various bleigifp parameters
   if ~isfield(results_EMEP_solvers, 'bleigifp')
      results_EMEP_solvers.bleigifp = cell(length(L_range), ...
         length(noise_ratio_range), length(EMEP_opts.bleigifp.block_size_range), ...
         length(EMEP_opts.bleigifp.inner_iter_range), length(EMEP_opts.bleigifp.num_eigs_range));
   end
   for L_idx = length(L_range):-1:1
      for noise_ratio_idx = length(noise_ratio_range):-1:1
         iter_end = results{L_idx, noise_ratio_idx}.saga_sd.solInfo.nit;
         for block_size_idx = 1:length(EMEP_opts.bleigifp.block_size_range)
            block_size = EMEP_opts.bleigifp.block_size_range(block_size_idx);
            for inner_iter_idx = 1:length(EMEP_opts.bleigifp.inner_iter_range)
               inner_iter = EMEP_opts.bleigifp.inner_iter_range(inner_iter_idx);
               for num_eigs_idx = 1:length(EMEP_opts.bleigifp.num_eigs_range)
                  
                  if isempty(results_EMEP_solvers.bleigifp{L_idx, noise_ratio_idx, block_size_idx, inner_iter_idx, num_eigs_idx})
                  
                     % iter max changed to enforce same number of mat-vecs per exp
                     iter_max_solver = round(mat_vec_max / (inner_iter*block_size) )
                     
                     num_eigs = EMEP_opts.bleigifp.num_eigs_range(num_eigs_idx);
                     exp_params = struct('num_eigs', num_eigs, 'block_size', block_size, ...
                        'inner_iter', inner_iter, 'iter_max_solver', iter_max_solver, ...
                        'iter_start', iter_start, 'iter_end', iter_end, 'disp', disp, ...
                        'solver', solver);
                     fprintf('Solving EMEP with bleigifp, L = %i, noise ratio = %1.3f, num eigs = %i, inner = %i, outer = %i\n', ...
                        L_range(L_idx), noise_ratio_range(noise_ratio_idx), num_eigs, inner_iter, block_size); 
                     results_EMEP_solvers.bleigifp{L_idx, noise_ratio_idx, block_size_idx, inner_iter_idx, num_eigs_idx} ...
                        = experiments.exputil.solve_EMEP(...
                           'results', results{L_idx, noise_ratio_idx}, exp_params);
                     results_EMEP_solvers.bleigifp{L_idx, noise_ratio_idx, block_size_idx, inner_iter_idx, num_eigs_idx}.exp_params ...
                        = exp_params;
                     results_EMEP_solvers.bleigifp{L_idx, noise_ratio_idx, block_size_idx, inner_iter_idx, num_eigs_idx}.exp_params.L = L_range(L_idx);
                     results_EMEP_solvers.bleigifp{L_idx, noise_ratio_idx, block_size_idx, inner_iter_idx, num_eigs_idx}.exp_params.noise_ratio = noise_ratio_range(noise_ratio_idx);
                     
                     save(strcat(folder_root, folder_exp, file_test), 'results_EMEP_solvers');                             
                     
                  end
               end
            end
         end
      end
   end
end

end




