function [data_noisy_128] = noisyrandom_Gaussian_var_dugap_stagnate_w_fig()

noise_ratio_range = [0.05, 0.15, 0.30];
L_range = [5, 10, 15];
n = 128;
signal = 'gaussian';
noise_type = 'gaussian';
iterations = 1000;

x_rel_err_array = zeros(length(L_range), length(noise_ratio_range));
num_total_iters_array = zeros(length(L_range), length(noise_ratio_range));
du_gap_array = zeros(length(L_range), length(noise_ratio_range));


% Loads cached (or solves and caches) results of test problem
pathName = 'cache/table_noisyrandom_Gaussian_var_dugap_stagnate_w_fig/';
if exist(pathName) ~= 7
   mkdir(pathName);
end


dataFile = fullfile(pathName, 'data_noisy_128.mat'); 
if exist(dataFile, 'file') ~= 2
   
   % PFD and DFP on
   matFile = fullfile(pathName, 'results_saga_noisy_128.mat');
   if exist(matFile, 'file') ~= 2
      results_saga_noisy_128 = cell(0);
      L_idx = 1;
      noise_ratio_idx = 1;

      for L = L_range
         for noise_ratio = noise_ratio_range
            results_saga_noisy_128{L_idx, noise_ratio_idx} ...
               = experiments.evolvingmatricespl('n', n, 'L', L, ...
                  'noise_ratio', noise_ratio, ...
                  'signal', signal, 'noise_type', noise_type, ...
                  'StopOnRelErr', false, 'iterations', iterations, ...
                  'cutDim', 10, 'eigs_basis_size', 40);
            noise_ratio_idx = noise_ratio_idx + 1;
         end
         noise_ratio_idx = 1;
         L_idx = L_idx + 1;
      end

      save(matFile, 'results_saga_noisy_128');
   else
      load(matFile, 'results_saga_noisy_128');
   end


   for L_idx = 1:length(L_range)
      for noise_ratio_idx = 1:length(noise_ratio_range)
         num_total_iters = results_saga_noisy_128{L_idx, noise_ratio_idx}.saga_sd.solInfo.nit;
         x_rel_err = results_saga_noisy_128{L_idx, noise_ratio_idx}.saga_sd.solInfo.iterData.x_error_rel(end);      
         du_gap = results_saga_noisy_128{L_idx, noise_ratio_idx}.saga_sd.solInfo.mGap;
         num_total_iters_array(L_idx, noise_ratio_idx) = num_total_iters;
         x_rel_err_array(L_idx, noise_ratio_idx) = x_rel_err;
         du_gap_array(L_idx, noise_ratio_idx) = du_gap;
      end
   end
end






dataFile = fullfile(pathName, 'data_noisy_128.mat'); 
if exist(dataFile, 'file') ~= 2
   data_noisy_128.num_total_iters_array = num_total_iters_array;
   data_noisy_128.x_rel_err_array = x_rel_err_array;
   data_noisy_128.du_gap_array = du_gap_array;
   
   L_idx = 2;
   noise_ratio_idx = 2;
   data_noisy_128.norm_b = norm(results_saga_noisy_128{L_idx, noise_ratio_idx}.gendatapl.b(:));
   data_noisy_128.res_primal = results_saga_noisy_128{L_idx, noise_ratio_idx}.saga_sd.solInfo.iterData.res_primal;
   data_noisy_128.nit = results_saga_noisy_128{L_idx, noise_ratio_idx}.saga_sd.solInfo.nit;
   
   data_noisy_128.res_mGap = results_saga_noisy_128{L_idx, noise_ratio_idx}.saga_sd.solInfo.iterData.res_mGap;
   
   save(dataFile, 'data_noisy_128');
   delete(fullfile(pathName, 'results_saga_noisy_128.mat'));
else
   load(dataFile, 'data_noisy_128');
end




% Prints results table and figure
fprintf('Total num iters for various values L, noise.\n')
fprintf('n = %5i     | L = %2i | L = %2i | L = %2i \n', n, L_range);  
fprintf('--------------|--------|--------|--------\n');
fprintf('noise = %1.2f  | %3.2f | %3.2f | %3.2f \n', noise_ratio_range(1), data_noisy_128.du_gap_array(:, 1));
fprintf('noise = %1.2f  | %3.2f | %3.2f | %3.2f \n', noise_ratio_range(2), data_noisy_128.du_gap_array(:, 2));
fprintf('noise = %1.2f  | %3.2f | %3.2f | %3.2f \n\n', noise_ratio_range(3), data_noisy_128.du_gap_array(:, 3));

fprintf('Relative signal error for various values L, noise.\n')
fprintf('n = %5i     | L = %2i    | L = %2i    | L = %2i \n', n, L_range);  
fprintf('--------------|-----------|-----------|------------\n');
fprintf('noise = %1.2f  | %1.3e | %1.3e | %1.3e  \n', noise_ratio_range(1), data_noisy_128.x_rel_err_array(:, 1));
fprintf('noise = %1.2f  | %1.3e | %1.3e | %1.3e  \n', noise_ratio_range(2), data_noisy_128.x_rel_err_array(:, 2));
fprintf('noise = %1.2f  | %1.3e | %1.3e | %1.3e  \n', noise_ratio_range(3), data_noisy_128.x_rel_err_array(:, 3));




L_idx = 2;
noise_ratio_idx = 2;
noise_ratio = noise_ratio_range(noise_ratio_idx);

% Plots relative error of signal obseration A(x_k) against true observation
subplot(1, 2, 1);
%norm_b = norm(results_saga_noisy_128{L_idx, noise_ratio_idx}.gendatapl.b(:));
%semilogx(results_saga_noisy_128{L_idx, noise_ratio_idx}.saga_sd.solInfo.iterData.res_primal / norm_b, 'b');
semilogx(data_noisy_128.res_primal / data_noisy_128.norm_b, 'b');
hold on
%plot(noise_ratio*ones(results_saga_noisy_128{L_idx, noise_ratio_idx}.saga_sd.solInfo.nit + 1, 1), 'r');
plot(noise_ratio*ones(data_noisy_128.nit + 1, 1), 'r');
hold off
title('Primal relative error');
ylabel('Relative error');
xlabel('Iteration');
xMin = 0; xMax = iterations + 1;
yMin = 0.10;%0.9*min(results_saga_noisy_128{L_idx, noise_ratio_idx}.saga_sd.solInfo.iterData.res_primal / norm_b);
yMax = 0.25;%1.1*max(results_saga_noisy_128{L_idx, noise_ratio_idx}.saga_sd.solInfo.iterData.res_primal / norm_b);
axis([xMin xMax yMin yMax]);
set(gca,'XTick',[1, logspace(1,4,4)])


% Plots strong duality gap
subplot(1, 2, 2);
%semilogx(results_saga_noisy_128{L_idx, noise_ratio_idx}.saga_sd.solInfo.iterData.res_mGap, 'b');
semilogx(data_noisy_128.res_mGap, 'b');
title('Duality gap');
ylabel('Error');
xlabel('Iteration');
xMin = 0; xMax = iterations + 1;
yMin = 0;
yMax = 1.05*max(data_noisy_128.res_mGap);
axis([xMin xMax yMin yMax]);
set(gca,'XTick',[1, logspace(1,4,4)])





end