function [data] = noisyrandom_synthetic_dual_variable_converges()

noise_ratio_range = [0.05, 0.15, 0.30];
L_range = [5, 10, 15];
n = 128;
signal = 'gaussian';
noise_type = 'dual';
iterations = 1;

x_rel_err_array = zeros(length(L_range), length(noise_ratio_range));
num_total_iters_array = zeros(length(L_range), length(noise_ratio_range));
x_rel_err_wflow_array = zeros(length(L_range), length(noise_ratio_range));

L_idx = 1;
noise_ratio_idx = 1;

for L = L_range
   for noise_ratio = noise_ratio_range
      results = experiments.evolvingmatricespl('n', n, 'L', L, ...
         'noise_ratio', noise_ratio, ...
         'signal', signal, 'noise_type', noise_type, ...
         'StopOnRelErr', false, 'iterations', iterations, ...
         'skipBBMinStepScale', 1e-10, 'skipBBMaxStepScale', 1e+10, ...
         'wflowOn', true);
      num_total_iters = results.saga_sd.solInfo.nit;
      x_rel_err = results.saga_sd.solInfo.iterData.x_error_rel(end);      
      num_total_iters_array(L_idx, noise_ratio_idx) = num_total_iters;
      x_rel_err_array(L_idx, noise_ratio_idx) = x_rel_err;
      x0 = results.gendatapl.x0;
      x_wflow = results.wflow.x;
      x_rel_err_wflow_array(L_idx, noise_ratio_idx) = util.hermitianerror(x_wflow(:), x0(:), 'fro') / norm(x0(:))^2;
      
      noise_ratio_idx = noise_ratio_idx + 1;
   end
   noise_ratio_idx = 1;
   L_idx = L_idx + 1;
end

data.num_total_iters_array = num_total_iters_array;
data.x_rel_err_array = x_rel_err_array;

fprintf('saga_sd total num iters for various values L, noise.\n')
fprintf('n = %5i     | L = %2i | L = %2i | L = %2i \n', n, L_range);  
fprintf('--------------|--------|--------|--------\n');
fprintf('noise = %1.2f  | %6i | %6i | %6i \n', noise_ratio_range(1), num_total_iters_array(:, 1));
fprintf('noise = %1.2f  | %6i | %6i | %6i \n', noise_ratio_range(2), num_total_iters_array(:, 2));
fprintf('noise = %1.2f  | %6i | %6i | %6i \n\n', noise_ratio_range(3), num_total_iters_array(:, 3));

fprintf('saga_sd relative signal error for various values L, noise.\n')
fprintf('n = %5i     | L = %2i    | L = %2i    | L = %2i \n', n, L_range);  
fprintf('--------------|-----------|-----------|------------\n');
fprintf('noise = %1.2f  | %1.3e | %1.3e | %1.3e  \n', noise_ratio_range(1), x_rel_err_array(:, 1));
fprintf('noise = %1.2f  | %1.3e | %1.3e | %1.3e  \n', noise_ratio_range(2), x_rel_err_array(:, 2));
fprintf('noise = %1.2f  | %1.3e | %1.3e | %1.3e  \n', noise_ratio_range(3), x_rel_err_array(:, 3));

fprintf('wflow relative signal error for various values L, noise.\n')
fprintf('n = %5i     | L = %2i    | L = %2i    | L = %2i \n', n, L_range);  
fprintf('--------------|-----------|-----------|------------\n');
fprintf('noise = %1.2f  | %1.3e | %1.3e | %1.3e  \n', noise_ratio_range(1), x_rel_err_wflow_array(:, 1));
fprintf('noise = %1.2f  | %1.3e | %1.3e | %1.3e  \n', noise_ratio_range(2), x_rel_err_wflow_array(:, 2));
fprintf('noise = %1.2f  | %1.3e | %1.3e | %1.3e  \n', noise_ratio_range(3), x_rel_err_wflow_array(:, 3));

end