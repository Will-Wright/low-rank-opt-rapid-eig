function [data_small_random, data_large_image] = noisyrandom_saga_vs_wflow_rel_errs(varargin)
% function [data] = noisyrandom_saga_vs_wflow_rel_errs(varargin)
%
% Generates phase retrieval results for random gaussian signals or an image
% using saga_sd and wflow.  Results are presented regarding signal relative
% error and probability of approximating true observation within a
% tolerance.

ip = inputParser;

ip.addParameter('data', []); % allows user to input results of previous call to experiment
ip.addParameter('test_type', 'both'); % determines test type, 
                  %valid options: 'small_random', 'large_image', or 'both'
ip.addParameter('print_tex_mode_on', false); % allows user to recreate test font for LaTeX
ip.addParameter('image_file', 'data/parrot_16k.jpg');
ip.addParameter('max_test_num', 100);
ip.addParameter('max_saga_iterations', 100);
ip.addParameter('resize_image', 1/1);
ip.parse(varargin{:});
test_type = ip.Results.test_type;
print_tex_mode_on = ip.Results.print_tex_mode_on;

data_small_random = [];
data_large_image = [];

folder_name = 'cache/table_noisyrandom_saga_vs_wflow_rel_errs';
if exist(folder_name) ~= 7
   mkdir(folder_name);
end  



% Generates test results for large image, wflow vs saga
if strcmp(test_type, 'small_random') || strcmp(test_type, 'both')
   % Searches for saved experiments
   if exist(strcat(folder_name, '/data_small_random.mat'))
      load(strcat(folder_name, '/data_small_random.mat'));
      
      wflow_Ax_rel_err_array = data_small_random.wflow_Ax_rel_err_array;
      saga_Ax_rel_err_array = data_small_random.saga_Ax_rel_err_array;
      wflow_x_rel_err_array = data_small_random.wflow_x_rel_err_array;
      saga_x_rel_err_array = data_small_random.saga_x_rel_err_array;
      saga_cosine_y_eta_angle_array = data_small_random.saga_cosine_y_eta_angle_array;
      n_range = data_small_random.n_range;
      L_range = data_small_random.L_range;
      noise_ratio_range = data_small_random.noise_ratio_range;
      max_test_num = data_small_random.max_test_num;
   
   else
      signal = 'gaussian';
      normalize_signal = false;
      max_test_num = ip.Results.max_test_num;
      image_file = ip.Results.image_file;

      n_range = 128;
      L_range = [4, 6, 8];
      noise_ratio_range = [0.05, 0.15, 0.30];

      max_saga_iterations = ip.Results.max_saga_iterations;

      wflow_x_rel_err_array = zeros(length(n_range), length(L_range), length(noise_ratio_range), max_test_num);
      wflow_Ax_rel_err_array = zeros(length(n_range), length(L_range), length(noise_ratio_range), max_test_num);
      saga_x_rel_err_array = zeros(length(n_range), length(L_range), length(noise_ratio_range), max_test_num);
      saga_Ax_rel_err_array = zeros(length(n_range), length(L_range), length(noise_ratio_range), max_test_num);
      saga_cosine_y_eta_angle_array = zeros(length(n_range), length(L_range), length(noise_ratio_range), max_test_num);


      seed = 1;

      for n = n_range
         for L = L_range
            for noise_ratio = noise_ratio_range
               for test_num = 1:max_test_num
                  results = experiments.evolvingmatricespl('signal', signal, ...
                     'optTol', 0, 'StopOnRelErr', false, ...
                     'cutDim', 5, 'image_file', image_file, ...
                     'wflowOn', true, 'noise_ratio', noise_ratio, 'L', L, 'n', n, ...
                     'iterations', max_saga_iterations, 'seed', seed, ...
                     'mask', 'gaussian', 'normalize_signal', normalize_signal);
                  b = results.gendatapl.b;
                  b0 = results.gendatapl.genInfo.b0;
                  eta = b(:) - b0(:);
                  x0 = results.gendatapl.x0;
                  y = results.saga_sd.solInfo.y(:);

                  n_idx = find(n == n_range);
                  L_idx = find(L == L_range);
                  noise_ratio_idx = find(noise_ratio == noise_ratio_range);

                  Ax_wflow = results.gendatapl.A.forward(results.wflow.x);
                  wflow_Ax_rel_err_array(n_idx, L_idx, noise_ratio_idx, test_num) ...
                     = norm(Ax_wflow(:) - b0(:))/norm(b0(:));
                  Ax_saga = results.gendatapl.A.forward(results.saga_sd.x);
                  saga_Ax_rel_err_array(n_idx, L_idx, noise_ratio_idx, test_num) ...
                     = norm(Ax_saga(:) - b0(:))/norm(b0(:));
                  wflow_x_rel_err_array(n_idx, L_idx, noise_ratio_idx, test_num) ...
                     = util.hermitianerror(results.wflow.x(:), x0(:), 'fro') / norm(x0(:))^2;
                  saga_x_rel_err_array(n_idx, L_idx, noise_ratio_idx, test_num) ...
                     = util.hermitianerror(results.saga_sd.x(:), x0(:), 'fro') / norm(x0(:))^2;
                  saga_cosine_y_eta_angle_array(n_idx, L_idx, noise_ratio_idx, test_num) ...
                     = eta'*y / (norm(eta)*norm(y));

                  seed = seed + 1;
               end
            end
         end
      end

      data_small_random.wflow_Ax_rel_err_array = wflow_Ax_rel_err_array;
      data_small_random.saga_Ax_rel_err_array = saga_Ax_rel_err_array;
      data_small_random.wflow_x_rel_err_array = wflow_x_rel_err_array;
      data_small_random.saga_x_rel_err_array = saga_x_rel_err_array;
      data_small_random.saga_cosine_y_eta_angle_array = saga_cosine_y_eta_angle_array;
      data_small_random.n_range = n_range;
      data_small_random.L_range = L_range;
      data_small_random.noise_ratio_range = noise_ratio_range;
      data_small_random.max_test_num = max_test_num;
      
      save(strcat(folder_name, '/data_small_random.mat'), 'data_small_random');
      
   end


   n_idx = 1;

   fprintf('\n');
   fprintf('           |              saga_sd             |            wflow\n');
   fprintf('  L  noise | cos(eta,y) x rel err  %% w/in tol | x rel err %% w/in tol \n');
   fprintf('     ratio |                         1.0  0.8 |             1.0  0.8\n');
   for L = L_range
      for noise_ratio = noise_ratio_range
         L_idx = find(L == L_range);
         noise_ratio_idx = find(noise_ratio == noise_ratio_range);
         saga_x_rel_err_mean = mean(saga_x_rel_err_array(n_idx, L_idx, noise_ratio_idx,:));
         saga_x_rel_err_std = std(saga_x_rel_err_array(n_idx, L_idx, noise_ratio_idx,:));
         saga_cosine_y_eta_angle_mean = mean(saga_cosine_y_eta_angle_array(n_idx, L_idx, noise_ratio_idx,:));
         saga_cosine_y_eta_angle_std = std(saga_cosine_y_eta_angle_array(n_idx, L_idx, noise_ratio_idx,:));
         tol = 1.0; saga_success_tol_10 = mean(saga_Ax_rel_err_array(n_idx, L_idx, noise_ratio_idx,:) <= tol*noise_ratio);
         tol = 0.8; saga_success_tol_08 = mean(saga_Ax_rel_err_array(n_idx, L_idx, noise_ratio_idx,:) <= tol*noise_ratio);

         wflow_x_rel_err_mean = mean(wflow_x_rel_err_array(n_idx, L_idx, noise_ratio_idx,:));
         wflow_x_rel_err_std = std(wflow_x_rel_err_array(n_idx, L_idx, noise_ratio_idx,:));
         tol = 1.0; wflow_success_tol_10 = mean(wflow_Ax_rel_err_array(n_idx, L_idx, noise_ratio_idx,:) <= tol*noise_ratio);
         tol = 0.8; wflow_success_tol_08 = mean(wflow_Ax_rel_err_array(n_idx, L_idx, noise_ratio_idx,:) <= tol*noise_ratio);
         if print_tex_mode_on
            fprintf('%3i & $%1.3f$ & $%1.2e$ & $%1.2e$ &  %1.2f & %1.2f & $%1.2e$ & %1.2f & %1.2f\n', ...
               L, noise_ratio, saga_cosine_y_eta_angle_mean, ...
               saga_x_rel_err_mean,  ...
               saga_success_tol_10, saga_success_tol_08, ...
               wflow_x_rel_err_mean, ...
               wflow_success_tol_10, wflow_success_tol_08);
         else
            fprintf('%3i  %1.3f |  %1.2e   %1.2e   %1.2f %1.2f | %1.2e   %1.2f %1.2f\n', ...
               L, noise_ratio, saga_cosine_y_eta_angle_mean, ...
               saga_x_rel_err_mean, ...
               saga_success_tol_10, saga_success_tol_08, ...
               wflow_x_rel_err_mean, ...
               wflow_success_tol_10, wflow_success_tol_08);
         end
      end
   end
   
end   
   




% Generates test results for large image, wflow vs saga
if strcmp(test_type, 'large_image') || strcmp(test_type, 'both')
   
   % Searches for saved experiments
   if exist(strcat(folder_name, '/data_large_image.mat'))
      load(strcat(folder_name, '/data_large_image.mat'));

      wflow_Ax_rel_err_array = data_large_image.wflow_Ax_rel_err_array;
      saga_Ax_rel_err_array = data_large_image.saga_Ax_rel_err_array;
      wflow_x_rel_err_array = data_large_image.wflow_x_rel_err_array;
      saga_x_rel_err_array = data_large_image.saga_x_rel_err_array;
      saga_cosine_y_eta_angle_array = data_large_image.saga_cosine_y_eta_angle_array;
      n_range = data_large_image.n_range;
      L_range = data_large_image.L_range;
      noise_ratio_range = data_large_image.noise_ratio_range;
      max_test_num = data_large_image.max_test_num; 
      
   else   
      signal = 'image';
      image_file = ip.Results.image_file;
      normalize_signal = false;
      max_test_num = 1;
      resize_image = ip.Results.resize_image;
      
      n_range = 1;
      max_test_num = 1;


      L_range = [4, 6, 8];
      noise_ratio_range = [0.05, 0.15, 0.30];

      max_saga_iterations = ip.Results.max_saga_iterations;

      wflow_x_rel_err_array = zeros(length(n_range), length(L_range), length(noise_ratio_range), max_test_num);
      wflow_Ax_rel_err_array = zeros(length(n_range), length(L_range), length(noise_ratio_range), max_test_num);
      saga_x_rel_err_array = zeros(length(n_range), length(L_range), length(noise_ratio_range), max_test_num);
      saga_Ax_rel_err_array = zeros(length(n_range), length(L_range), length(noise_ratio_range), max_test_num);
      saga_cosine_y_eta_angle_array = zeros(length(n_range), length(L_range), length(noise_ratio_range), max_test_num);


      seed = 1;

      for n = n_range
         for L = L_range
            for noise_ratio = noise_ratio_range
               for test_num = 1:max_test_num
               
                  n_idx = find(n == n_range);
                  L_idx = find(L == L_range);
                  noise_ratio_idx = find(noise_ratio == noise_ratio_range);
                  
                  data_large_image.results{n_idx, L_idx, noise_ratio_idx} = experiments.evolvingmatricespl('signal', signal, ...
                     'optTol', 0, 'StopOnRelErr', false, ...
                     'cutDim', 5, 'image_file', image_file, ...
                     'wflowOn', true, 'noise_ratio', noise_ratio, 'L', L, 'n', n, ...
                     'iterations', max_saga_iterations, 'seed', seed, ...
                     'mask', 'gaussian', 'normalize_signal', normalize_signal, ...
                     'logIterData', false, 'channel', 'RGB', 'resizeImage', resize_image);

                  b_red = data_large_image.results{n_idx, L_idx, noise_ratio_idx}.red.gendatapl.b;
                  b_green = data_large_image.results{n_idx, L_idx, noise_ratio_idx}.green.gendatapl.b;
                  b_blue = data_large_image.results{n_idx, L_idx, noise_ratio_idx}.blue.gendatapl.b;
                  b0_red = data_large_image.results{n_idx, L_idx, noise_ratio_idx}.red.gendatapl.genInfo.b0;
                  b0_green = data_large_image.results{n_idx, L_idx, noise_ratio_idx}.green.gendatapl.genInfo.b0;
                  b0_blue = data_large_image.results{n_idx, L_idx, noise_ratio_idx}.blue.gendatapl.genInfo.b0;
                  eta_red = b_red(:) - b0_red(:);
                  eta_green = b_green(:) - b0_green(:);
                  eta_blue = b_blue(:) - b0_blue(:);
                  x0_red = data_large_image.results{n_idx, L_idx, noise_ratio_idx}.red.gendatapl.x0;
                  x0_green = data_large_image.results{n_idx, L_idx, noise_ratio_idx}.green.gendatapl.x0;
                  x0_blue = data_large_image.results{n_idx, L_idx, noise_ratio_idx}.blue.gendatapl.x0;
                  y_red = data_large_image.results{n_idx, L_idx, noise_ratio_idx}.red.saga_sd.solInfo.y(:);
                  y_green = data_large_image.results{n_idx, L_idx, noise_ratio_idx}.green.saga_sd.solInfo.y(:);
                  y_blue = data_large_image.results{n_idx, L_idx, noise_ratio_idx}.blue.saga_sd.solInfo.y(:);

                  Ax_wflow_red = data_large_image.results{n_idx, L_idx, noise_ratio_idx}.red.gendatapl.A.forward(...
                     data_large_image.results{n_idx, L_idx, noise_ratio_idx}.red.wflow.x);
                  Ax_wflow_green = data_large_image.results{n_idx, L_idx, noise_ratio_idx}.green.gendatapl.A.forward(...
                     data_large_image.results{n_idx, L_idx, noise_ratio_idx}.green.wflow.x);
                  Ax_wflow_blue = data_large_image.results{n_idx, L_idx, noise_ratio_idx}.blue.gendatapl.A.forward(...
                     data_large_image.results{n_idx, L_idx, noise_ratio_idx}.blue.wflow.x);
                  wflow_Ax_rel_err_array(n_idx, L_idx, noise_ratio_idx, test_num) ...
                     = (1/3)*(norm(Ax_wflow_red(:) - b0_red(:))/norm(b0_red(:)) ...
                     + norm(Ax_wflow_green(:) - b0_green(:))/norm(b0_green(:)) ...
                     + norm(Ax_wflow_blue(:) - b0_blue(:))/norm(b0_blue(:)));
                  
                  Ax_saga_red = data_large_image.results{n_idx, L_idx, noise_ratio_idx}.red.gendatapl.A.forward(...
                     data_large_image.results{n_idx, L_idx, noise_ratio_idx}.red.saga_sd.x);
                  Ax_saga_green = data_large_image.results{n_idx, L_idx, noise_ratio_idx}.green.gendatapl.A.forward(...
                     data_large_image.results{n_idx, L_idx, noise_ratio_idx}.green.saga_sd.x);
                  Ax_saga_blue = data_large_image.results{n_idx, L_idx, noise_ratio_idx}.blue.gendatapl.A.forward(...
                     data_large_image.results{n_idx, L_idx, noise_ratio_idx}.blue.saga_sd.x);
                  saga_Ax_rel_err_array(n_idx, L_idx, noise_ratio_idx, test_num) ...
                     = (1/3)*(norm(Ax_saga_red(:) - b0_red(:))/norm(b0_red(:)) ...
                     + norm(Ax_saga_green(:) - b0_green(:))/norm(b0_green(:)) ...
                     + norm(Ax_saga_blue(:) - b0_blue(:))/norm(b0_blue(:)));
                  
                  wflow_x_rel_err_array(n_idx, L_idx, noise_ratio_idx, test_num) ...
                     = (1/3)*(util.hermitianerror(data_large_image.results{n_idx, L_idx, noise_ratio_idx}.red.wflow.x(:), x0_red(:), 'fro') / norm(x0_red(:))^2 ...
                     + util.hermitianerror(data_large_image.results{n_idx, L_idx, noise_ratio_idx}.green.wflow.x(:), x0_green(:), 'fro') / norm(x0_green(:))^2 ...
                     + util.hermitianerror(data_large_image.results{n_idx, L_idx, noise_ratio_idx}.blue.wflow.x(:), x0_blue(:), 'fro') / norm(x0_blue(:))^2);
                  saga_x_rel_err_array(n_idx, L_idx, noise_ratio_idx, test_num) ...
                     = (1/3)*(util.hermitianerror(data_large_image.results{n_idx, L_idx, noise_ratio_idx}.red.saga_sd.x(:), x0_red(:), 'fro') / norm(x0_red(:))^2 ...
                     + util.hermitianerror(data_large_image.results{n_idx, L_idx, noise_ratio_idx}.green.saga_sd.x(:), x0_green(:), 'fro') / norm(x0_green(:))^2 ...
                     + util.hermitianerror(data_large_image.results{n_idx, L_idx, noise_ratio_idx}.blue.saga_sd.x(:), x0_blue(:), 'fro') / norm(x0_blue(:))^2);
                  
                  saga_cosine_y_eta_angle_array(n_idx, L_idx, noise_ratio_idx, test_num) ...
                     = (1/3)*(eta_red'*y_red / (norm(eta_red)*norm(y_red)) ...
                     + eta_green'*y_green / (norm(eta_green)*norm(y_green)) ...
                     + eta_blue'*y_blue / (norm(eta_blue)*norm(y_blue)));

                  data_large_image.wflow_Ax_rel_err_array = wflow_Ax_rel_err_array;
                  data_large_image.saga_Ax_rel_err_array = saga_Ax_rel_err_array;
                  data_large_image.wflow_x_rel_err_array = wflow_x_rel_err_array;
                  data_large_image.saga_x_rel_err_array = saga_x_rel_err_array;
                  data_large_image.saga_cosine_y_eta_angle_array = saga_cosine_y_eta_angle_array;
                  data_large_image.n_range = n_range;
                  data_large_image.L_range = L_range;
                  data_large_image.noise_ratio_range = noise_ratio_range;
                  data_large_image.max_test_num = max_test_num;

                  
                  data_large_image = rmfield(data_large_image, 'results');
                  
                  save(strcat(folder_name, '/data_large_image.mat'), 'data_large_image');

                  
                  seed = seed + 1;
               end
            end
         end
      end
      
   end

   
   n_idx = 1;

   fprintf('\n');
   fprintf('           |              saga_sd             |            wflow\n');
   fprintf('  L  noise | cos(eta,y) x rel err  %% w/in tol | x rel err %% w/in tol \n');
   fprintf('     ratio |                         1.0  0.8 |             1.0  0.8\n');
   for L = L_range
      for noise_ratio = noise_ratio_range
         L_idx = find(L == L_range);
         noise_ratio_idx = find(noise_ratio == noise_ratio_range);
         saga_x_rel_err_mean = mean(saga_x_rel_err_array(n_idx, L_idx, noise_ratio_idx,:));
         saga_x_rel_err_std = std(saga_x_rel_err_array(n_idx, L_idx, noise_ratio_idx,:));
         saga_cosine_y_eta_angle_mean = mean(saga_cosine_y_eta_angle_array(n_idx, L_idx, noise_ratio_idx,:));
         saga_cosine_y_eta_angle_std = std(saga_cosine_y_eta_angle_array(n_idx, L_idx, noise_ratio_idx,:));
         tol = 1.0; saga_success_tol_10 = mean(saga_Ax_rel_err_array(n_idx, L_idx, noise_ratio_idx,:) <= tol*noise_ratio);
         tol = 0.8; saga_success_tol_08 = mean(saga_Ax_rel_err_array(n_idx, L_idx, noise_ratio_idx,:) <= tol*noise_ratio);

         wflow_x_rel_err_mean = mean(wflow_x_rel_err_array(n_idx, L_idx, noise_ratio_idx,:));
         wflow_x_rel_err_std = std(wflow_x_rel_err_array(n_idx, L_idx, noise_ratio_idx,:));
         tol = 1.0; wflow_success_tol_10 = mean(wflow_Ax_rel_err_array(n_idx, L_idx, noise_ratio_idx,:) <= tol*noise_ratio);
         tol = 0.8; wflow_success_tol_08 = mean(wflow_Ax_rel_err_array(n_idx, L_idx, noise_ratio_idx,:) <= tol*noise_ratio);
         if print_tex_mode_on
            fprintf('%3i & $%1.3f$ & $%1.2e$ & $%1.2e$ &  %1.2f & %1.2f & $%1.2e$ & %1.2f & %1.2f\n', ...
               L, noise_ratio, saga_cosine_y_eta_angle_mean, ...
               saga_x_rel_err_mean,  ...
               saga_success_tol_10, saga_success_tol_08, ...
               wflow_x_rel_err_mean, ...
               wflow_success_tol_10, wflow_success_tol_08);
         else
            fprintf('%3i  %1.3f |  %1.2e   %1.2e   %1.2f %1.2f | %1.2e   %1.2f %1.2f\n', ...
               L, noise_ratio, saga_cosine_y_eta_angle_mean, ...
               saga_x_rel_err_mean, ...
               saga_success_tol_10, saga_success_tol_08, ...
               wflow_x_rel_err_mean, ...
               wflow_success_tol_10, wflow_success_tol_08);
         end
      end
   end
   
else
   fprintf('Invalid test type called.\n');
end



end