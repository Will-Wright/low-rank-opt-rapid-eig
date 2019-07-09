function [data] = noisyrandom_mean_rank_lifted_signal_exact_solution
% experiments.table.noisyrandom_mean_rank_lifted_signal_exact_solution
% 
% Solves noisy phase problems with saga and cvx to compare results
% and determine average rank of optimal matrix X.

folder_name = 'cache/table_noisyrandom_mean_rank_lifted_signal_exact_solution/';

noise_ratio_range = [0.05, 0.15, 0.30];
L_range = [4, 6, 8];
signal_type_range = ["dual", "gaussian"]; % 'dual' is synthetic noise, 'gaussian' is natural noise
n = 16;
iterations = 1000;
noise_type = 'gaussian';
max_num_tests = 100;
% rank1 and rank2 match, redundant measurement from both primal and dual SDPs
mean_rank1 = zeros(length(L_range), length(noise_ratio_range), length(signal_type_range));
mean_rank2 = zeros(length(L_range), length(noise_ratio_range), length(signal_type_range));
rank1 = zeros(length(L_range), length(noise_ratio_range), length(signal_type_range), max_num_tests);
rank2 = zeros(length(L_range), length(noise_ratio_range), length(signal_type_range), max_num_tests);
num_iters_array = zeros(length(L_range), length(noise_ratio_range), length(signal_type_range), max_num_tests);
mean_num_iters = zeros(length(L_range), length(noise_ratio_range), length(signal_type_range));
mGap_array = zeros(length(L_range), length(noise_ratio_range), length(signal_type_range), max_num_tests);
mean_mGap = zeros(length(L_range), length(noise_ratio_range), length(signal_type_range));

L_idx = 1;
noise_ratio_idx = 1;
signal_idx = 1;
seed = 1;



% Loads cached (or solves and caches) results of test problem
if exist(folder_name) ~= 7
   mkdir(folder_name);
end
matFile = fullfile(folder_name, 'data.mat');
if exist(matFile, 'file') ~= 2
   % Solves set of problems to generate 
   for signal = signal_type_range
      for L = L_range
         for noise_ratio = noise_ratio_range
            for i = 1:max_num_tests
               [results, data_temp] = experiments.exputil.solve_noisyimage_saga_vs_exact_SDP_solutions(...
                  'n', n, 'L', L, 'noise_ratio', noise_ratio, ...
                  'noise_type', noise_type, 'signal', signal, 'iterations', iterations, ...
                  'seed', seed);
               rank1(L_idx, noise_ratio_idx, signal_idx, i) = data_temp.rank_X_cvx_primal;
               rank2(L_idx, noise_ratio_idx, signal_idx, i) = data_temp.rank_cvx_primal;      
               num_iters_array(L_idx, noise_ratio_idx, signal_idx, i) = results.saga_sd.solInfo.nit;
               mGap_array(L_idx, noise_ratio_idx, signal_idx, i) = results.saga_sd.solInfo.mGap;
               seed = seed + 1;
            end
            mean_rank1(L_idx, noise_ratio_idx, signal_idx) = mean(rank1(L_idx, noise_ratio_idx, signal_idx, :));
            mean_rank2(L_idx, noise_ratio_idx, signal_idx) = mean(rank2(L_idx, noise_ratio_idx, signal_idx, :));
            mean_num_iters(L_idx, noise_ratio_idx, signal_idx) = mean(num_iters_array(L_idx, noise_ratio_idx, signal_idx, :));
            mean_mGap(L_idx, noise_ratio_idx, signal_idx) = mean(mGap_array(L_idx, noise_ratio_idx, signal_idx, :));
            noise_ratio_idx = noise_ratio_idx + 1;
         end
         noise_ratio_idx = 1;
         L_idx = L_idx + 1;
      end
      L_idx = 1;
      signal_idx = signal_idx + 1;
      % Resets seed so same set of models are used for both experiments
      seed = 1;
   end
   
   data.rank1 = rank1;
   data.rank2 = rank2;
   data.avg_rank1 = mean_rank1;
   data.avg_rank2 = mean_rank2;
   data.num_iters_array = num_iters_array;
   data.mean_num_iters = mean_num_iters;   
   data.mGap_array = mGap_array;
   data.mean_mGap = mean_mGap;   
   save(matFile, 'data');
else
   load(matFile, 'data');
end




% Old test code to verify results above match proper saga results
%{
max_total_tests = 3*3*max_num_tests;
for L = L_range
   for noise_ratio = noise_ratio_range
      for i = 1:max_num_tests
         [results] = experiments.evolvingmatricespl(...
            'n', n, 'L', L, 'noise_ratio', noise_ratio, ...
            'noise_type', noise_type, 'signal', signal, 'iterations', iterations, ...
            'seed', seed, 'StopOnRelErr', false, 'verbosity', false);
         num_iters_array(L_idx, noise_ratio_idx, i) = results.saga_sd.solInfo.nit;
         mGap_array(L_idx, noise_ratio_idx, i) = results.saga_sd.solInfo.mGap;
         seed = seed + 1;
         fprintf('test num %i out of %i\n', seed-1, max_total_tests);
      end
      
      noise_ratio_idx = noise_ratio_idx + 1;
   end
   noise_ratio_idx = 1;
   L_idx = L_idx + 1;
end

data.num_iters_array2 = num_iters_array;
data.mGap_array2 = mGap_array;
%}

% Prints table of results for sequence of solutions
fprintf('\n                               |         Synthetic noise        |         Natural noise\n')
fprintf('              n = %2i           ', n);
for L_idx = 1:length(L_range)
   fprintf('|  L = %2i  ', L_range(L_idx));
end
for L_idx = 1:length(L_range)
   fprintf('|  L = %2i  ', L_range(L_idx));
end
fprintf('\n')
for noise_ratio_idx = 1:length(noise_ratio_range)
   fprintf('---------------------------------------------------------------------------------------------------\n');
   fprintf('                       rank(X) ')
   for L_idx = 1:length(L_range)
      fprintf('|     %2.2f ', data.avg_rank1(L_idx, noise_ratio_idx, 1));
   end
   for L_idx = 1:length(L_range)
      fprintf('|     %2.2f ', data.avg_rank1(L_idx, noise_ratio_idx, 2));
   end   
   fprintf('\n')
   fprintf('noise_ratio = %1.2f    saga its ', noise_ratio_range(noise_ratio_idx));
   for L_idx = 1:length(L_range)
      fprintf('|    %4.2f ', data.mean_num_iters(L_idx, noise_ratio_idx, 1));
   end
   for L_idx = 1:length(L_range)
      fprintf('|    %4.2f ', data.mean_num_iters(L_idx, noise_ratio_idx, 2));
   end
   fprintf('\n')   
   fprintf('                      dual gap ')
	for L_idx = 1:length(L_range)
      fprintf('| %1.2e ', data.mean_mGap(L_idx, noise_ratio_idx, 1));
   end
	for L_idx = 1:length(L_range)
      fprintf('|    %2.2f ', data.mean_mGap(L_idx, noise_ratio_idx, 2));
   end   
   fprintf('\n')   
end








end