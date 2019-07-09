function [data] = noisyimage_comparison_adaptive_vs_orig()
% function [data] = noisyimage_comparison_adaptive_vs_orig()
%
% Plots figures showing performance of `eigs`
% solving EMEP with original eigs settings and adaptive method


folder_root = 'cache/figure_noisyimage_comparison_adaptive_vs_orig/';
if exist(folder_root) ~= 7
   mkdir(folder_root)
end

if ~exist(strcat(folder_root, 'data_gaussian.mat'))
   solve_gaussian_phase_problems(folder_root);
   process_gaussian_data(folder_root);
end

if ~exist(strcat(folder_root, 'data_image.mat'))
   solve_image_phase_problems(folder_root);
   process_image_data(folder_root);
end

data = plot_figures_and_tables(folder_root);


end





function [] = solve_gaussian_phase_problems(folder_root)

% Generates results for gaussian signal for various n
signal = 'gaussian';
folder_exp = 'gaussian_n_range/';
n_range = [128, 256, 512, 1024];
max_num_tests = 10;
noise_ratio = 0.15;

if exist(strcat(folder_root, folder_exp)) ~= 7
   mkdir(strcat(folder_root, folder_exp));
end
if exist(strcat(folder_root, folder_exp, 'results.mat')) == 2
   load(strcat(folder_root, folder_exp, 'results.mat'), 'results')
else
   results = struct;
   results.orig = cell(length(n_range),max_num_tests);
   results.ada40 = cell(length(n_range),max_num_tests);
end

fprintf('Solving phase problems: gaussian signals with various size n.\n');
for i = 1:length(n_range)
   n = n_range(i);
   L = 2*round(log(n));
   for test_num = 1:max_num_tests
      seed = test_num;
      
      if isempty(results.ada40{i, test_num})
         results.ada40{i, test_num} = experiments.evolvingmatricespl(...
            'signal', signal, 'n', n, 'L', L, 'noise_ratio', noise_ratio, ...
            'seed', seed, 'eigIts', 100000, 'maxTime', 48*60*60,...
            'useAdaptiveEigs', true, 'eigs_basis_size', 40, ...
            'logIterData', false);
         save(strcat(folder_root, folder_exp, 'results.mat'), 'results');
      end
      
      if isempty(results.orig{i, test_num})
         results.orig{i, test_num} = experiments.evolvingmatricespl(...
            'signal', signal, 'n', n, 'L', L, 'noise_ratio', noise_ratio, ...
            'seed', seed, 'eigIts', 100000, 'maxTime', 48*60*60,...
            'useAdaptiveEigs', false, 'cutDim', 2, 'eigs_basis_size', 20, ...
            'logIterData', false);
         save(strcat(folder_root, folder_exp, 'results.mat'), 'results');
      end
   end
end 





% Generates results for gaussian signal for various noise_ratio
signal = 'gaussian';
folder_exp = 'gaussian_noise_ratio_range/';
noise_ratio_range = [0.05, 0.10, 0.15, 0.20, 0.25, 0.30];
max_num_tests = 10;
n = 128;
L = 10;

if exist(strcat(folder_root, folder_exp)) ~= 7
   mkdir(strcat(folder_root, folder_exp));
end
if exist(strcat(folder_root, folder_exp, 'results.mat')) == 2
   load(strcat(folder_root, folder_exp, 'results.mat'), 'results')
else
   results = struct;
   results.orig = cell(length(noise_ratio_range),max_num_tests);
   results.ada40 = cell(length(noise_ratio_range),max_num_tests);
end

fprintf('Solving phase problems: gaussian signals with various noise ratio.\n');
for i = 1:length(noise_ratio_range)
   noise_ratio = noise_ratio_range(i);
   for test_num = 1:max_num_tests
      seed = test_num;
      
      if isempty(results.ada40{i, test_num})
         results.ada40{i, test_num} = experiments.evolvingmatricespl(...
            'signal', signal, 'n', n, 'L', L, 'noise_ratio', noise_ratio, ...
            'seed', seed, 'eigIts', 100000, 'maxTime', 48*60*60,...
            'useAdaptiveEigs', true, 'eigs_basis_size', 40, ...
            'logIterData', false);
         save(strcat(folder_root, folder_exp, 'results.mat'), 'results');
      end
      
      if isempty(results.orig{i, test_num})
         results.orig{i, test_num} = experiments.evolvingmatricespl(...
            'signal', signal, 'n', n, 'L', L, 'noise_ratio', noise_ratio, ...
            'seed', seed, 'eigIts', 100000, 'maxTime', 48*60*60,...
            'useAdaptiveEigs', false, 'cutDim', 2, 'eigs_basis_size', 20, ...
            'logIterData', false);
         save(strcat(folder_root, folder_exp, 'results.mat'), 'results');
      end
   end
end 
   





% Generates results for gaussian signal for various L
signal = 'gaussian';
folder_exp = 'gaussian_L_range/';
L_range = [8, 12, 16, 20];
max_num_tests = 10;
n = 128;
noise_ratio = 0.15;

if exist(strcat(folder_root, folder_exp)) ~= 7
   mkdir(strcat(folder_root, folder_exp));
end
if exist(strcat(folder_root, folder_exp, 'results.mat')) == 2
   load(strcat(folder_root, folder_exp, 'results.mat'), 'results')
else
   results = struct;
   results.orig = cell(length(L_range),max_num_tests);
   results.ada40 = cell(length(L_range),max_num_tests);
end

fprintf('Solving phase problems: gaussian signals with various oversampling L.\n');
for i = 1:length(L_range)
   L = L_range(i);
   for test_num = 1:max_num_tests
      seed = test_num;
      
      if isempty(results.ada40{i, test_num})
         results.ada40{i, test_num} = experiments.evolvingmatricespl(...
            'signal', signal, 'n', n, 'L', L, 'noise_ratio', noise_ratio, ...
            'seed', seed, 'eigIts', 100000, 'maxTime', 48*60*60,...
            'useAdaptiveEigs', true, 'eigs_basis_size', 40, ...
            'logIterData', false);
         save(strcat(folder_root, folder_exp, 'results.mat'), 'results');
      end
      
      if isempty(results.orig{i, test_num})
         results.orig{i, test_num} = experiments.evolvingmatricespl(...
            'signal', signal, 'n', n, 'L', L, 'noise_ratio', noise_ratio, ...
            'seed', seed, 'eigIts', 100000, 'maxTime', 48*60*60,...
            'useAdaptiveEigs', false, 'cutDim', 2, 'eigs_basis_size', 20, ...
            'logIterData', false);
         save(strcat(folder_root, folder_exp, 'results.mat'), 'results');
      end
   end
end 

end




function [] = process_gaussian_data(folder_root)

signal = 'gaussian';
folder_exp = 'gaussian_n_range/';
data.n_exp.n_range = [128, 256, 512, 1024];
data.n_exp.max_num_tests = 10;
data.n_exp.noise_ratio = 0.15;

load(strcat(folder_root, folder_exp, 'results.mat'), 'results')
[n, num_tests] = size(results.orig);

for n_idx = 1:n
   for test_idx = 1:num_tests

      mvs_orig = results.orig{n_idx, test_idx}.saga_sd.solInfo.iterData.infoEigsFFTs_y{1,1}.nMatVec;
      nit = length(results.orig{n_idx, test_idx}.saga_sd.solInfo.iterData.infoEigsFFTs_yLS);
      for i = 1:nit
         if ~isempty(results.orig{n_idx, test_idx}.saga_sd.solInfo.iterData.infoEigsFFTs_yLS{i})
            for j = 1:length(results.orig{n_idx, test_idx}.saga_sd.solInfo.iterData.infoEigsFFTs_yLS{i})
               mvs_orig = mvs_orig + results.orig{n_idx, test_idx}.saga_sd.solInfo.iterData.infoEigsFFTs_yLS{i}{1,j}.nMatVec;
            end
         end
      end


      mvs_ada40 = results.ada40{n_idx, test_idx}.saga_sd.solInfo.iterData.infoEigsFFTs_y{1,1}.nMatVec;
      nit = length(results.ada40{n_idx, test_idx}.saga_sd.solInfo.iterData.infoEigsFFTs_yLS);
      for i = 1:nit
         if ~isempty(results.ada40{n_idx, test_idx}.saga_sd.solInfo.iterData.infoEigsFFTs_yLS{i})
            for j = 1:length(results.ada40{n_idx, test_idx}.saga_sd.solInfo.iterData.infoEigsFFTs_yLS{i})
               mvs_ada40 = mvs_ada40 + results.ada40{n_idx, test_idx}.saga_sd.solInfo.iterData.infoEigsFFTs_yLS{i}{1,j}.nMatVec;
            end
         end
      end

      data.n_exp.mvs_orig(n_idx, test_idx) = mvs_orig;
      data.n_exp.mvs_ada40(n_idx, test_idx) = mvs_ada40;
   end
end





signal = 'gaussian';
folder_exp = 'gaussian_noise_ratio_range/';
data.nr_exp.noise_ratio_range = [0.05, 0.10, 0.15, 0.20, 0.25, 0.30];
data.nr_exp.max_num_tests = 10;
data.nr_exp.n = 128;
data.nr_exp.L = 10;

load(strcat(folder_root, folder_exp, 'results.mat'), 'results')
[nr, num_tests] = size(results.orig);

for nr_idx = 1:nr
   for test_idx = 1:num_tests

      mvs_orig = results.orig{nr_idx, test_idx}.saga_sd.solInfo.iterData.infoEigsFFTs_y{1,1}.nMatVec;
      nit = length(results.orig{nr_idx, test_idx}.saga_sd.solInfo.iterData.infoEigsFFTs_yLS);
      for i = 1:nit
         if ~isempty(results.orig{nr_idx, test_idx}.saga_sd.solInfo.iterData.infoEigsFFTs_yLS{i})
            for j = 1:length(results.orig{nr_idx, test_idx}.saga_sd.solInfo.iterData.infoEigsFFTs_yLS{i})
               mvs_orig = mvs_orig + results.orig{nr_idx, test_idx}.saga_sd.solInfo.iterData.infoEigsFFTs_yLS{i}{1,j}.nMatVec;
            end
         end
      end


      mvs_ada40 = results.ada40{nr_idx, test_idx}.saga_sd.solInfo.iterData.infoEigsFFTs_y{1,1}.nMatVec;
      nit = length(results.ada40{nr_idx, test_idx}.saga_sd.solInfo.iterData.infoEigsFFTs_yLS);
      for i = 1:nit
         if ~isempty(results.ada40{nr_idx, test_idx}.saga_sd.solInfo.iterData.infoEigsFFTs_yLS{i})
            for j = 1:length(results.ada40{nr_idx, test_idx}.saga_sd.solInfo.iterData.infoEigsFFTs_yLS{i})
               mvs_ada40 = mvs_ada40 + results.ada40{nr_idx, test_idx}.saga_sd.solInfo.iterData.infoEigsFFTs_yLS{i}{1,j}.nMatVec;
            end
         end
      end

      data.nr_exp.mvs_orig(nr_idx, test_idx) = mvs_orig;
      data.nr_exp.mvs_ada40(nr_idx, test_idx) = mvs_ada40;
   end
end






signal = 'gaussian';
folder_exp = 'gaussian_L_range/';
data.L_exp.L_range = [8, 12, 16, 20];
data.L_exp.max_num_tests = 10;
data.L_exp.n = 128;
data.L_exp.noise_ratio = 0.15;

load(strcat(folder_root, folder_exp, 'results.mat'), 'results')
[L, num_tests] = size(results.orig);

for L_idx = 1:L
   for test_idx = 1:num_tests

      mvs_orig = results.orig{L_idx, test_idx}.saga_sd.solInfo.iterData.infoEigsFFTs_y{1,1}.nMatVec;
      nit = length(results.orig{L_idx, test_idx}.saga_sd.solInfo.iterData.infoEigsFFTs_yLS);
      for i = 1:nit
         if ~isempty(results.orig{L_idx, test_idx}.saga_sd.solInfo.iterData.infoEigsFFTs_yLS{i})
            for j = 1:length(results.orig{L_idx, test_idx}.saga_sd.solInfo.iterData.infoEigsFFTs_yLS{i})
               mvs_orig = mvs_orig + results.orig{L_idx, test_idx}.saga_sd.solInfo.iterData.infoEigsFFTs_yLS{i}{1,j}.nMatVec;
            end
         end
      end

      mvs_ada40 = results.ada40{L_idx, test_idx}.saga_sd.solInfo.iterData.infoEigsFFTs_y{1,1}.nMatVec;
      nit = length(results.ada40{L_idx, test_idx}.saga_sd.solInfo.iterData.infoEigsFFTs_yLS);
      for i = 1:nit
         if ~isempty(results.ada40{L_idx, test_idx}.saga_sd.solInfo.iterData.infoEigsFFTs_yLS{i})
            for j = 1:length(results.ada40{L_idx, test_idx}.saga_sd.solInfo.iterData.infoEigsFFTs_yLS{i})
               mvs_ada40 = mvs_ada40 + results.ada40{L_idx, test_idx}.saga_sd.solInfo.iterData.infoEigsFFTs_yLS{i}{1,j}.nMatVec;
            end
         end
      end
      
      data.L_exp.mvs_orig(L_idx, test_idx) = mvs_orig;
      data.L_exp.mvs_ada40(L_idx, test_idx) = mvs_ada40;
   end
end


% Deletes large `results` files and saves `data` file
save(strcat(folder_root, 'data_gaussian.mat'), 'data');
delete(strcat(folder_root, 'gaussian_n_range/', 'results.mat'))
delete(strcat(folder_root, 'gaussian_noise_ratio_range/', 'results.mat'))
delete(strcat(folder_root, 'gaussian_L_range/', 'results.mat'))
rmdir(strcat(folder_root, 'gaussian_n_range/'))
rmdir(strcat(folder_root, 'gaussian_noise_ratio_range/'))
rmdir(strcat(folder_root, 'gaussian_L_range/'))

end





function [] = solve_image_phase_problems(folder_root)


% Generates results for experiment using 1st image
signal = 'image';
image_file = 'data/jul_and_me_small.jpg';
folder_exp = 'jul_and_me_small/';
L_range = [10, 15];
noise_ratio_range = [0.15];

if exist(strcat(folder_root, folder_exp)) ~= 7
   mkdir(strcat(folder_root, folder_exp));
end
if exist(strcat(folder_root, folder_exp, 'results.mat')) == 2
   load(strcat(folder_root, folder_exp, 'results.mat'), 'results')
else
   results = struct;
   results.orig = cell(length(L_range),length(noise_ratio_range));
   results.ada40 = cell(length(L_range),length(noise_ratio_range));
   results.ada80 = cell(length(L_range),length(noise_ratio_range));
end

fprintf('Solving phase problems: image signals 1 of 2.\n');
for noise_ratio_idx = 1:length(noise_ratio_range)
   for L_idx = 1:length(L_range)
   L = L_range(L_idx);
      noise_ratio = noise_ratio_range(noise_ratio_idx);
      
      if isempty(results.ada40{L_idx, noise_ratio_idx})
         results.ada40{L_idx, noise_ratio_idx} = experiments.evolvingmatricespl(...
            'signal', signal, 'L', L, 'noise_ratio', noise_ratio, ...
            'image_file', image_file, 'eigIts', 100000, 'maxTime', 48*60*60,...
            'useAdaptiveEigs', true, 'eigs_basis_size', 40, 'resizeImage', 1/1, ...
            'channel', 'RGB', 'logIterData', false);
         save(strcat(folder_root, folder_exp, 'results.mat'), 'results');
      end
      
      if isempty(results.ada80{L_idx, noise_ratio_idx})
         results.ada80{L_idx, noise_ratio_idx} = experiments.evolvingmatricespl(...
            'signal', signal, 'L', L, 'noise_ratio', noise_ratio, ...
            'image_file', image_file, 'eigIts', 100000, 'maxTime', 48*60*60,...
            'useAdaptiveEigs', true, 'eigs_basis_size', 80, 'resizeImage', 1/1, ...
            'channel', 'RGB', 'logIterData', false);
         save(strcat(folder_root, folder_exp, 'results.mat'), 'results');
      end
      
      if isempty(results.orig{L_idx, noise_ratio_idx})
         results.orig{L_idx, noise_ratio_idx} = experiments.evolvingmatricespl(...
            'signal', signal, 'L', L, 'noise_ratio', noise_ratio, ...
            'image_file', image_file, 'eigIts', 100000, 'maxTime', 48*60*60,...
            'useAdaptiveEigs', false, 'cutDim', 2, 'eigs_basis_size', 20, 'resizeImage', 1/1, ...
            'channel', 'RGB', 'logIterData', false);
         save(strcat(folder_root, folder_exp, 'results.mat'), 'results');
      end
      
   end
end 




% Generates results for experiment using 2nd image
signal = 'image';
image_file = 'data/UC_Davis_small.jpeg';
folder_exp = 'UC_Davis_small/';
L_range = [10, 15];
noise_ratio_range = [0.15];

if exist(strcat(folder_root, folder_exp)) ~= 7
   mkdir(strcat(folder_root, folder_exp));
end
if exist(strcat(folder_root, folder_exp, 'results.mat')) == 2
   load(strcat(folder_root, folder_exp, 'results.mat'), 'results')
else
   results = struct;
   results.orig = cell(length(L_range),length(noise_ratio_range));
   results.ada40 = cell(length(L_range),length(noise_ratio_range));
   results.ada80 = cell(length(L_range),length(noise_ratio_range));
end

fprintf('Solving phase problems: image signals 2 of 2.\n');
for noise_ratio_idx = 1:length(noise_ratio_range)
   for L_idx = 1:length(L_range)
   L = L_range(L_idx);
      noise_ratio = noise_ratio_range(noise_ratio_idx);
      
      if isempty(results.ada40{L_idx, noise_ratio_idx})
         results.ada40{L_idx, noise_ratio_idx} = experiments.evolvingmatricespl(...
            'signal', signal, 'L', L, 'noise_ratio', noise_ratio, ...
            'image_file', image_file, 'eigIts', 100000, 'maxTime', 48*60*60,...
            'useAdaptiveEigs', true, 'eigs_basis_size', 40, 'resizeImage', 1/1, ...
            'channel', 'RGB', 'logIterData', false);
         save(strcat(folder_root, folder_exp, 'results.mat'), 'results');
      end
      
      if isempty(results.ada80{L_idx, noise_ratio_idx})
         results.ada80{L_idx, noise_ratio_idx} = experiments.evolvingmatricespl(...
            'signal', signal, 'L', L, 'noise_ratio', noise_ratio, ...
            'image_file', image_file, 'eigIts', 100000, 'maxTime', 48*60*60,...
            'useAdaptiveEigs', true, 'eigs_basis_size', 80, 'resizeImage', 1/1, ...
            'channel', 'RGB', 'logIterData', false);
         save(strcat(folder_root, folder_exp, 'results.mat'), 'results');
      end
      
      if isempty(results.orig{L_idx, noise_ratio_idx})
         results.orig{L_idx, noise_ratio_idx} = experiments.evolvingmatricespl(...
            'signal', signal, 'L', L, 'noise_ratio', noise_ratio, ...
            'image_file', image_file, 'eigIts', 100000, 'maxTime', 48*60*60,...
            'useAdaptiveEigs', false, 'cutDim', 2, 'eigs_basis_size', 20, 'resizeImage', 1/1, ...
            'channel', 'RGB', 'logIterData', false);
         save(strcat(folder_root, folder_exp, 'results.mat'), 'results');
      end
      
   end
end 

end




function [] = process_image_data(folder_root)

folder_exp = 'jul_and_me_small/';
load(strcat(folder_root, folder_exp, 'results.mat'), 'results')

data.jul_and_me_exp.L_range = [10, 15];
[L, nr] = size(results.orig);

for L_idx = 1:L
   for nr_idx = 1:nr

      mvs_orig_red = results.orig{L_idx, nr_idx}.red.saga_sd.solInfo.iterData.infoEigsFFTs_y{1,1}.nMatVec;
      nit = length(results.orig{L_idx, nr_idx}.red.saga_sd.solInfo.iterData.infoEigsFFTs_yLS);
      for i = 1:nit
         if ~isempty(results.orig{L_idx, nr_idx}.red.saga_sd.solInfo.iterData.infoEigsFFTs_yLS{i})
            for j = 1:length(results.orig{L_idx, nr_idx}.red.saga_sd.solInfo.iterData.infoEigsFFTs_yLS{i})
               mvs_orig_red = mvs_orig_red + results.orig{L_idx, nr_idx}.red.saga_sd.solInfo.iterData.infoEigsFFTs_yLS{i}{1,j}.nMatVec;
            end
         end
      end

      mvs_orig_green = results.orig{L_idx, nr_idx}.green.saga_sd.solInfo.iterData.infoEigsFFTs_y{1,1}.nMatVec;
      nit = length(results.orig{L_idx, nr_idx}.green.saga_sd.solInfo.iterData.infoEigsFFTs_yLS);
      for i = 1:nit
         if ~isempty(results.orig{L_idx, nr_idx}.green.saga_sd.solInfo.iterData.infoEigsFFTs_yLS{i})
            for j = 1:length(results.orig{L_idx, nr_idx}.green.saga_sd.solInfo.iterData.infoEigsFFTs_yLS{i})
               mvs_orig_green = mvs_orig_green + results.orig{L_idx, nr_idx}.green.saga_sd.solInfo.iterData.infoEigsFFTs_yLS{i}{1,j}.nMatVec;
            end
         end
      end

      mvs_orig_blue = results.orig{L_idx, nr_idx}.blue.saga_sd.solInfo.iterData.infoEigsFFTs_y{1,1}.nMatVec;
      nit = length(results.orig{L_idx, nr_idx}.blue.saga_sd.solInfo.iterData.infoEigsFFTs_yLS);
      for i = 1:nit
         if ~isempty(results.orig{L_idx, nr_idx}.blue.saga_sd.solInfo.iterData.infoEigsFFTs_yLS{i})
            for j = 1:length(results.orig{L_idx, nr_idx}.blue.saga_sd.solInfo.iterData.infoEigsFFTs_yLS{i})
               mvs_orig_blue = mvs_orig_blue + results.orig{L_idx, nr_idx}.blue.saga_sd.solInfo.iterData.infoEigsFFTs_yLS{i}{1,j}.nMatVec;
            end
         end
      end


      mvs_ada40_red = results.ada40{L_idx, nr_idx}.red.saga_sd.solInfo.iterData.infoEigsFFTs_y{1,1}.nMatVec;
      nit = length(results.ada40{L_idx, nr_idx}.red.saga_sd.solInfo.iterData.infoEigsFFTs_yLS);
      for i = 1:nit
         if ~isempty(results.ada40{L_idx, nr_idx}.red.saga_sd.solInfo.iterData.infoEigsFFTs_yLS{i})
            for j = 1:length(results.ada40{L_idx, nr_idx}.red.saga_sd.solInfo.iterData.infoEigsFFTs_yLS{i})
               mvs_ada40_red = mvs_ada40_red + results.ada40{L_idx, nr_idx}.red.saga_sd.solInfo.iterData.infoEigsFFTs_yLS{i}{1,j}.nMatVec;
            end
         end
      end

      mvs_ada40_green = results.ada40{L_idx, nr_idx}.green.saga_sd.solInfo.iterData.infoEigsFFTs_y{1,1}.nMatVec;
      nit = length(results.ada40{L_idx, nr_idx}.green.saga_sd.solInfo.iterData.infoEigsFFTs_yLS);
      for i = 1:nit
         if ~isempty(results.ada40{L_idx, nr_idx}.green.saga_sd.solInfo.iterData.infoEigsFFTs_yLS{i})
            for j = 1:length(results.ada40{L_idx, nr_idx}.green.saga_sd.solInfo.iterData.infoEigsFFTs_yLS{i})
               mvs_ada40_green = mvs_ada40_green + results.ada40{L_idx, nr_idx}.green.saga_sd.solInfo.iterData.infoEigsFFTs_yLS{i}{1,j}.nMatVec;
            end
         end
      end

      mvs_ada40_blue = results.ada40{L_idx, nr_idx}.blue.saga_sd.solInfo.iterData.infoEigsFFTs_y{1,1}.nMatVec;
      nit = length(results.ada40{L_idx, nr_idx}.blue.saga_sd.solInfo.iterData.infoEigsFFTs_yLS);
      for i = 1:nit
         if ~isempty(results.ada40{L_idx, nr_idx}.blue.saga_sd.solInfo.iterData.infoEigsFFTs_yLS{i})
            for j = 1:length(results.ada40{L_idx, nr_idx}.blue.saga_sd.solInfo.iterData.infoEigsFFTs_yLS{i})
               mvs_ada40_blue = mvs_ada40_blue + results.ada40{L_idx, nr_idx}.blue.saga_sd.solInfo.iterData.infoEigsFFTs_yLS{i}{1,j}.nMatVec;
            end
         end
      end


      mvs_ada80_red = results.ada80{L_idx, nr_idx}.red.saga_sd.solInfo.iterData.infoEigsFFTs_y{1,1}.nMatVec;
      nit = length(results.ada80{L_idx, nr_idx}.red.saga_sd.solInfo.iterData.infoEigsFFTs_yLS);
      for i = 1:nit
         if ~isempty(results.ada80{L_idx, nr_idx}.red.saga_sd.solInfo.iterData.infoEigsFFTs_yLS{i})
            for j = 1:length(results.ada80{L_idx, nr_idx}.red.saga_sd.solInfo.iterData.infoEigsFFTs_yLS{i})
               mvs_ada80_red = mvs_ada80_red + results.ada80{L_idx, nr_idx}.red.saga_sd.solInfo.iterData.infoEigsFFTs_yLS{i}{1,j}.nMatVec;
            end
         end
      end

      mvs_ada80_green = results.ada80{L_idx, nr_idx}.green.saga_sd.solInfo.iterData.infoEigsFFTs_y{1,1}.nMatVec;
      nit = length(results.ada80{L_idx, nr_idx}.green.saga_sd.solInfo.iterData.infoEigsFFTs_yLS);
      for i = 1:nit
         if ~isempty(results.ada80{L_idx, nr_idx}.green.saga_sd.solInfo.iterData.infoEigsFFTs_yLS{i})
            for j = 1:length(results.ada80{L_idx, nr_idx}.green.saga_sd.solInfo.iterData.infoEigsFFTs_yLS{i})
               mvs_ada80_green = mvs_ada80_green + results.ada80{L_idx, nr_idx}.green.saga_sd.solInfo.iterData.infoEigsFFTs_yLS{i}{1,j}.nMatVec;
            end
         end
      end

      mvs_ada80_blue = results.ada80{L_idx, nr_idx}.blue.saga_sd.solInfo.iterData.infoEigsFFTs_y{1,1}.nMatVec;
      nit = length(results.ada80{L_idx, nr_idx}.blue.saga_sd.solInfo.iterData.infoEigsFFTs_yLS);
      for i = 1:nit
         if ~isempty(results.ada80{L_idx, nr_idx}.blue.saga_sd.solInfo.iterData.infoEigsFFTs_yLS{i})
            for j = 1:length(results.ada80{L_idx, nr_idx}.blue.saga_sd.solInfo.iterData.infoEigsFFTs_yLS{i})
               mvs_ada80_blue = mvs_ada80_blue + results.ada80{L_idx, nr_idx}.blue.saga_sd.solInfo.iterData.infoEigsFFTs_yLS{i}{1,j}.nMatVec;
            end
         end
      end

      data.jul_and_me_exp.mvs_orig(L_idx, nr_idx) = mvs_orig_red + mvs_orig_green + mvs_orig_blue;
      data.jul_and_me_exp.mvs_ada40(L_idx, nr_idx) = mvs_ada40_red + mvs_ada40_green + mvs_ada40_blue;
      data.jul_and_me_exp.mvs_ada80(L_idx, nr_idx) = mvs_ada80_red + mvs_ada80_green + mvs_ada80_blue;
   end
end


% Deletes large `results` files and saves `data` file
results.orig{1,1}.red.saga_sd.solInfo = rmfield(results.orig{1,1}.red.saga_sd.solInfo, 'y');
results.orig{1,1}.red.saga_sd = rmfield(results.orig{1,1}.red.saga_sd, 'r');
results.orig{1,1}.red.gendatapl = rmfield(results.orig{1,1}.red.gendatapl, 'A');
results.orig{1,1}.red.gendatapl = rmfield(results.orig{1,1}.red.gendatapl, 'b');
results.orig{1,1}.red.gendatapl.genInfo = rmfield(results.orig{1,1}.red.gendatapl.genInfo, 'b0');
results.orig{2,1}.red.saga_sd.solInfo = rmfield(results.orig{2,1}.red.saga_sd.solInfo, 'y');
results.orig{2,1}.red.saga_sd = rmfield(results.orig{2,1}.red.saga_sd, 'r');
results.orig{2,1}.red.gendatapl = rmfield(results.orig{2,1}.red.gendatapl, 'A');
results.orig{2,1}.red.gendatapl = rmfield(results.orig{2,1}.red.gendatapl, 'b');
results.orig{2,1}.red.gendatapl.genInfo = rmfield(results.orig{2,1}.red.gendatapl.genInfo, 'b0');

results.orig{1,1}.green.saga_sd.solInfo = rmfield(results.orig{1,1}.green.saga_sd.solInfo, 'y');
results.orig{1,1}.green.saga_sd = rmfield(results.orig{1,1}.green.saga_sd, 'r');
results.orig{1,1}.green.gendatapl = rmfield(results.orig{1,1}.green.gendatapl, 'A');
results.orig{1,1}.green.gendatapl = rmfield(results.orig{1,1}.green.gendatapl, 'b');
results.orig{1,1}.green.gendatapl.genInfo = rmfield(results.orig{1,1}.green.gendatapl.genInfo, 'b0');
results.orig{2,1}.green.saga_sd.solInfo = rmfield(results.orig{2,1}.green.saga_sd.solInfo, 'y');
results.orig{2,1}.green.saga_sd = rmfield(results.orig{2,1}.green.saga_sd, 'r');
results.orig{2,1}.green.gendatapl = rmfield(results.orig{2,1}.green.gendatapl, 'A');
results.orig{2,1}.green.gendatapl = rmfield(results.orig{2,1}.green.gendatapl, 'b');
results.orig{2,1}.green.gendatapl.genInfo = rmfield(results.orig{2,1}.green.gendatapl.genInfo, 'b0');

results.orig{1,1}.blue.saga_sd.solInfo = rmfield(results.orig{1,1}.blue.saga_sd.solInfo, 'y');
results.orig{1,1}.blue.saga_sd = rmfield(results.orig{1,1}.blue.saga_sd, 'r');
results.orig{1,1}.blue.gendatapl = rmfield(results.orig{1,1}.blue.gendatapl, 'A');
results.orig{1,1}.blue.gendatapl = rmfield(results.orig{1,1}.blue.gendatapl, 'b');
results.orig{1,1}.blue.gendatapl.genInfo = rmfield(results.orig{1,1}.blue.gendatapl.genInfo, 'b0');
results.orig{2,1}.blue.saga_sd.solInfo = rmfield(results.orig{2,1}.blue.saga_sd.solInfo, 'y');
results.orig{2,1}.blue.saga_sd = rmfield(results.orig{2,1}.blue.saga_sd, 'r');
results.orig{2,1}.blue.gendatapl = rmfield(results.orig{2,1}.blue.gendatapl, 'A');
results.orig{2,1}.blue.gendatapl = rmfield(results.orig{2,1}.blue.gendatapl, 'b');
results.orig{2,1}.blue.gendatapl.genInfo = rmfield(results.orig{2,1}.blue.gendatapl.genInfo, 'b0');


results.ada40{1,1}.red.saga_sd.solInfo = rmfield(results.ada40{1,1}.red.saga_sd.solInfo, 'y');
results.ada40{1,1}.red.saga_sd = rmfield(results.ada40{1,1}.red.saga_sd, 'r');
results.ada40{1,1}.red.gendatapl = rmfield(results.ada40{1,1}.red.gendatapl, 'A');
results.ada40{1,1}.red.gendatapl = rmfield(results.ada40{1,1}.red.gendatapl, 'b');
results.ada40{1,1}.red.gendatapl.genInfo = rmfield(results.ada40{1,1}.red.gendatapl.genInfo, 'b0');
results.ada40{2,1}.red.saga_sd.solInfo = rmfield(results.ada40{2,1}.red.saga_sd.solInfo, 'y');
results.ada40{2,1}.red.saga_sd = rmfield(results.ada40{2,1}.red.saga_sd, 'r');
results.ada40{2,1}.red.gendatapl = rmfield(results.ada40{2,1}.red.gendatapl, 'A');
results.ada40{2,1}.red.gendatapl = rmfield(results.ada40{2,1}.red.gendatapl, 'b');
results.ada40{2,1}.red.gendatapl.genInfo = rmfield(results.ada40{2,1}.red.gendatapl.genInfo, 'b0');

results.ada40{1,1}.green.saga_sd.solInfo = rmfield(results.ada40{1,1}.green.saga_sd.solInfo, 'y');
results.ada40{1,1}.green.saga_sd = rmfield(results.ada40{1,1}.green.saga_sd, 'r');
results.ada40{1,1}.green.gendatapl = rmfield(results.ada40{1,1}.green.gendatapl, 'A');
results.ada40{1,1}.green.gendatapl = rmfield(results.ada40{1,1}.green.gendatapl, 'b');
results.ada40{1,1}.green.gendatapl.genInfo = rmfield(results.ada40{1,1}.green.gendatapl.genInfo, 'b0');
results.ada40{2,1}.green.saga_sd.solInfo = rmfield(results.ada40{2,1}.green.saga_sd.solInfo, 'y');
results.ada40{2,1}.green.saga_sd = rmfield(results.ada40{2,1}.green.saga_sd, 'r');
results.ada40{2,1}.green.gendatapl = rmfield(results.ada40{2,1}.green.gendatapl, 'A');
results.ada40{2,1}.green.gendatapl = rmfield(results.ada40{2,1}.green.gendatapl, 'b');
results.ada40{2,1}.green.gendatapl.genInfo = rmfield(results.ada40{2,1}.green.gendatapl.genInfo, 'b0');

results.ada40{1,1}.blue.saga_sd.solInfo = rmfield(results.ada40{1,1}.blue.saga_sd.solInfo, 'y');
results.ada40{1,1}.blue.saga_sd = rmfield(results.ada40{1,1}.blue.saga_sd, 'r');
results.ada40{1,1}.blue.gendatapl = rmfield(results.ada40{1,1}.blue.gendatapl, 'A');
results.ada40{1,1}.blue.gendatapl = rmfield(results.ada40{1,1}.blue.gendatapl, 'b');
results.ada40{1,1}.blue.gendatapl.genInfo = rmfield(results.ada40{1,1}.blue.gendatapl.genInfo, 'b0');
results.ada40{2,1}.blue.saga_sd.solInfo = rmfield(results.ada40{2,1}.blue.saga_sd.solInfo, 'y');
results.ada40{2,1}.blue.saga_sd = rmfield(results.ada40{2,1}.blue.saga_sd, 'r');
results.ada40{2,1}.blue.gendatapl = rmfield(results.ada40{2,1}.blue.gendatapl, 'A');
results.ada40{2,1}.blue.gendatapl = rmfield(results.ada40{2,1}.blue.gendatapl, 'b');
results.ada40{2,1}.blue.gendatapl.genInfo = rmfield(results.ada40{2,1}.blue.gendatapl.genInfo, 'b0');


results.ada80{1,1}.red.saga_sd.solInfo = rmfield(results.ada80{1,1}.red.saga_sd.solInfo, 'y');
results.ada80{1,1}.red.saga_sd = rmfield(results.ada80{1,1}.red.saga_sd, 'r');
results.ada80{1,1}.red.gendatapl = rmfield(results.ada80{1,1}.red.gendatapl, 'A');
results.ada80{1,1}.red.gendatapl = rmfield(results.ada80{1,1}.red.gendatapl, 'b');
results.ada80{1,1}.red.gendatapl.genInfo = rmfield(results.ada80{1,1}.red.gendatapl.genInfo, 'b0');
results.ada80{2,1}.red.saga_sd.solInfo = rmfield(results.ada80{2,1}.red.saga_sd.solInfo, 'y');
results.ada80{2,1}.red.saga_sd = rmfield(results.ada80{2,1}.red.saga_sd, 'r');
results.ada80{2,1}.red.gendatapl = rmfield(results.ada80{2,1}.red.gendatapl, 'A');
results.ada80{2,1}.red.gendatapl = rmfield(results.ada80{2,1}.red.gendatapl, 'b');
results.ada80{2,1}.red.gendatapl.genInfo = rmfield(results.ada80{2,1}.red.gendatapl.genInfo, 'b0');

results.ada80{1,1}.green.saga_sd.solInfo = rmfield(results.ada80{1,1}.green.saga_sd.solInfo, 'y');
results.ada80{1,1}.green.saga_sd = rmfield(results.ada80{1,1}.green.saga_sd, 'r');
results.ada80{1,1}.green.gendatapl = rmfield(results.ada80{1,1}.green.gendatapl, 'A');
results.ada80{1,1}.green.gendatapl = rmfield(results.ada80{1,1}.green.gendatapl, 'b');
results.ada80{1,1}.green.gendatapl.genInfo = rmfield(results.ada80{1,1}.green.gendatapl.genInfo, 'b0');
results.ada80{2,1}.green.saga_sd.solInfo = rmfield(results.ada80{2,1}.green.saga_sd.solInfo, 'y');
results.ada80{2,1}.green.saga_sd = rmfield(results.ada80{2,1}.green.saga_sd, 'r');
results.ada80{2,1}.green.gendatapl = rmfield(results.ada80{2,1}.green.gendatapl, 'A');
results.ada80{2,1}.green.gendatapl = rmfield(results.ada80{2,1}.green.gendatapl, 'b');
results.ada80{2,1}.green.gendatapl.genInfo = rmfield(results.ada80{2,1}.green.gendatapl.genInfo, 'b0');

results.ada80{1,1}.blue.saga_sd.solInfo = rmfield(results.ada80{1,1}.blue.saga_sd.solInfo, 'y');
results.ada80{1,1}.blue.saga_sd = rmfield(results.ada80{1,1}.blue.saga_sd, 'r');
results.ada80{1,1}.blue.gendatapl = rmfield(results.ada80{1,1}.blue.gendatapl, 'A');
results.ada80{1,1}.blue.gendatapl = rmfield(results.ada80{1,1}.blue.gendatapl, 'b');
results.ada80{1,1}.blue.gendatapl.genInfo = rmfield(results.ada80{1,1}.blue.gendatapl.genInfo, 'b0');
results.ada80{2,1}.blue.saga_sd.solInfo = rmfield(results.ada80{2,1}.blue.saga_sd.solInfo, 'y');
results.ada80{2,1}.blue.saga_sd = rmfield(results.ada80{2,1}.blue.saga_sd, 'r');
results.ada80{2,1}.blue.gendatapl = rmfield(results.ada80{2,1}.blue.gendatapl, 'A');
results.ada80{2,1}.blue.gendatapl = rmfield(results.ada80{2,1}.blue.gendatapl, 'b');
results.ada80{2,1}.blue.gendatapl.genInfo = rmfield(results.ada80{2,1}.blue.gendatapl.genInfo, 'b0');

data.jul_and_me_results = results;

[~, ~, imOrig, imSaga] = experiments.exputil.viewrecoveredimageandresidual(data.jul_and_me_results.ada40{2,1});
imwrite(imOrig, strcat(folder_root, 'jul_and_me_orig.png'));
imwrite(imSaga, strcat(folder_root, 'jul_and_me_L_15_ada40.png'));
save(strcat(folder_root, 'data_image.mat'), 'data');
delete(strcat(folder_root, 'jul_and_me_small/', 'results.mat'))
rmdir(strcat(folder_root, 'jul_and_me_small/'))




folder_exp = 'UC_Davis_small/';
load(strcat(folder_root, folder_exp, 'results.mat'), 'results')

data.UCD_exp.L_range = [10, 15];
[L, nr] = size(results.orig);


for L_idx = 1:L
   for nr_idx = 1:nr

      mvs_orig_red = results.orig{L_idx, nr_idx}.red.saga_sd.solInfo.iterData.infoEigsFFTs_y{1,1}.nMatVec;
      nit = length(results.orig{L_idx, nr_idx}.red.saga_sd.solInfo.iterData.infoEigsFFTs_yLS);
      for i = 1:nit
         if ~isempty(results.orig{L_idx, nr_idx}.red.saga_sd.solInfo.iterData.infoEigsFFTs_yLS{i})
            for j = 1:length(results.orig{L_idx, nr_idx}.red.saga_sd.solInfo.iterData.infoEigsFFTs_yLS{i})
               mvs_orig_red = mvs_orig_red + results.orig{L_idx, nr_idx}.red.saga_sd.solInfo.iterData.infoEigsFFTs_yLS{i}{1,j}.nMatVec;
            end
         end
      end

      mvs_orig_green = results.orig{L_idx, nr_idx}.green.saga_sd.solInfo.iterData.infoEigsFFTs_y{1,1}.nMatVec;
      nit = length(results.orig{L_idx, nr_idx}.green.saga_sd.solInfo.iterData.infoEigsFFTs_yLS);
      for i = 1:nit
         if ~isempty(results.orig{L_idx, nr_idx}.green.saga_sd.solInfo.iterData.infoEigsFFTs_yLS{i})
            for j = 1:length(results.orig{L_idx, nr_idx}.green.saga_sd.solInfo.iterData.infoEigsFFTs_yLS{i})
               mvs_orig_green = mvs_orig_green + results.orig{L_idx, nr_idx}.green.saga_sd.solInfo.iterData.infoEigsFFTs_yLS{i}{1,j}.nMatVec;
            end
         end
      end

      mvs_orig_blue = results.orig{L_idx, nr_idx}.blue.saga_sd.solInfo.iterData.infoEigsFFTs_y{1,1}.nMatVec;
      nit = length(results.orig{L_idx, nr_idx}.blue.saga_sd.solInfo.iterData.infoEigsFFTs_yLS);
      for i = 1:nit
         if ~isempty(results.orig{L_idx, nr_idx}.blue.saga_sd.solInfo.iterData.infoEigsFFTs_yLS{i})
            for j = 1:length(results.orig{L_idx, nr_idx}.blue.saga_sd.solInfo.iterData.infoEigsFFTs_yLS{i})
               mvs_orig_blue = mvs_orig_blue + results.orig{L_idx, nr_idx}.blue.saga_sd.solInfo.iterData.infoEigsFFTs_yLS{i}{1,j}.nMatVec;
            end
         end
      end


      mvs_ada40_red = results.ada40{L_idx, nr_idx}.red.saga_sd.solInfo.iterData.infoEigsFFTs_y{1,1}.nMatVec;
      nit = length(results.ada40{L_idx, nr_idx}.red.saga_sd.solInfo.iterData.infoEigsFFTs_yLS);
      for i = 1:nit
         if ~isempty(results.ada40{L_idx, nr_idx}.red.saga_sd.solInfo.iterData.infoEigsFFTs_yLS{i})
            for j = 1:length(results.ada40{L_idx, nr_idx}.red.saga_sd.solInfo.iterData.infoEigsFFTs_yLS{i})
               mvs_ada40_red = mvs_ada40_red + results.ada40{L_idx, nr_idx}.red.saga_sd.solInfo.iterData.infoEigsFFTs_yLS{i}{1,j}.nMatVec;
            end
         end
      end

      mvs_ada40_green = results.ada40{L_idx, nr_idx}.green.saga_sd.solInfo.iterData.infoEigsFFTs_y{1,1}.nMatVec;
      nit = length(results.ada40{L_idx, nr_idx}.green.saga_sd.solInfo.iterData.infoEigsFFTs_yLS);
      for i = 1:nit
         if ~isempty(results.ada40{L_idx, nr_idx}.green.saga_sd.solInfo.iterData.infoEigsFFTs_yLS{i})
            for j = 1:length(results.ada40{L_idx, nr_idx}.green.saga_sd.solInfo.iterData.infoEigsFFTs_yLS{i})
               mvs_ada40_green = mvs_ada40_green + results.ada40{L_idx, nr_idx}.green.saga_sd.solInfo.iterData.infoEigsFFTs_yLS{i}{1,j}.nMatVec;
            end
         end
      end

      mvs_ada40_blue = results.ada40{L_idx, nr_idx}.blue.saga_sd.solInfo.iterData.infoEigsFFTs_y{1,1}.nMatVec;
      nit = length(results.ada40{L_idx, nr_idx}.blue.saga_sd.solInfo.iterData.infoEigsFFTs_yLS);
      for i = 1:nit
         if ~isempty(results.ada40{L_idx, nr_idx}.blue.saga_sd.solInfo.iterData.infoEigsFFTs_yLS{i})
            for j = 1:length(results.ada40{L_idx, nr_idx}.blue.saga_sd.solInfo.iterData.infoEigsFFTs_yLS{i})
               mvs_ada40_blue = mvs_ada40_blue + results.ada40{L_idx, nr_idx}.blue.saga_sd.solInfo.iterData.infoEigsFFTs_yLS{i}{1,j}.nMatVec;
            end
         end
      end


      mvs_ada80_red = results.ada80{L_idx, nr_idx}.red.saga_sd.solInfo.iterData.infoEigsFFTs_y{1,1}.nMatVec;
      nit = length(results.ada80{L_idx, nr_idx}.red.saga_sd.solInfo.iterData.infoEigsFFTs_yLS);
      for i = 1:nit
         if ~isempty(results.ada80{L_idx, nr_idx}.red.saga_sd.solInfo.iterData.infoEigsFFTs_yLS{i})
            for j = 1:length(results.ada80{L_idx, nr_idx}.red.saga_sd.solInfo.iterData.infoEigsFFTs_yLS{i})
               mvs_ada80_red = mvs_ada80_red + results.ada80{L_idx, nr_idx}.red.saga_sd.solInfo.iterData.infoEigsFFTs_yLS{i}{1,j}.nMatVec;
            end
         end
      end

      mvs_ada80_green = results.ada80{L_idx, nr_idx}.green.saga_sd.solInfo.iterData.infoEigsFFTs_y{1,1}.nMatVec;
      nit = length(results.ada80{L_idx, nr_idx}.green.saga_sd.solInfo.iterData.infoEigsFFTs_yLS);
      for i = 1:nit
         if ~isempty(results.ada80{L_idx, nr_idx}.green.saga_sd.solInfo.iterData.infoEigsFFTs_yLS{i})
            for j = 1:length(results.ada80{L_idx, nr_idx}.green.saga_sd.solInfo.iterData.infoEigsFFTs_yLS{i})
               mvs_ada80_green = mvs_ada80_green + results.ada80{L_idx, nr_idx}.green.saga_sd.solInfo.iterData.infoEigsFFTs_yLS{i}{1,j}.nMatVec;
            end
         end
      end

      mvs_ada80_blue = results.ada80{L_idx, nr_idx}.blue.saga_sd.solInfo.iterData.infoEigsFFTs_y{1,1}.nMatVec;
      nit = length(results.ada80{L_idx, nr_idx}.blue.saga_sd.solInfo.iterData.infoEigsFFTs_yLS);
      for i = 1:nit
         if ~isempty(results.ada80{L_idx, nr_idx}.blue.saga_sd.solInfo.iterData.infoEigsFFTs_yLS{i})
            for j = 1:length(results.ada80{L_idx, nr_idx}.blue.saga_sd.solInfo.iterData.infoEigsFFTs_yLS{i})
               mvs_ada80_blue = mvs_ada80_blue + results.ada80{L_idx, nr_idx}.blue.saga_sd.solInfo.iterData.infoEigsFFTs_yLS{i}{1,j}.nMatVec;
            end
         end
      end

      data.UCD_exp.mvs_orig(L_idx, nr_idx) = mvs_orig_red + mvs_orig_green + mvs_orig_blue;
      data.UCD_exp.mvs_ada40(L_idx, nr_idx) = mvs_ada40_red + mvs_ada40_green + mvs_ada40_blue;
      data.UCD_exp.mvs_ada80(L_idx, nr_idx) = mvs_ada80_red + mvs_ada80_green + mvs_ada80_blue;
   end
end


% Deletes large `results` files and saves `data` file
results.orig{1,1}.red.saga_sd.solInfo = rmfield(results.orig{1,1}.red.saga_sd.solInfo, 'y');
results.orig{1,1}.red.saga_sd = rmfield(results.orig{1,1}.red.saga_sd, 'r');
results.orig{1,1}.red.gendatapl = rmfield(results.orig{1,1}.red.gendatapl, 'A');
results.orig{1,1}.red.gendatapl = rmfield(results.orig{1,1}.red.gendatapl, 'b');
results.orig{1,1}.red.gendatapl.genInfo = rmfield(results.orig{1,1}.red.gendatapl.genInfo, 'b0');
results.orig{2,1}.red.saga_sd.solInfo = rmfield(results.orig{2,1}.red.saga_sd.solInfo, 'y');
results.orig{2,1}.red.saga_sd = rmfield(results.orig{2,1}.red.saga_sd, 'r');
results.orig{2,1}.red.gendatapl = rmfield(results.orig{2,1}.red.gendatapl, 'A');
results.orig{2,1}.red.gendatapl = rmfield(results.orig{2,1}.red.gendatapl, 'b');
results.orig{2,1}.red.gendatapl.genInfo = rmfield(results.orig{2,1}.red.gendatapl.genInfo, 'b0');

results.orig{1,1}.green.saga_sd.solInfo = rmfield(results.orig{1,1}.green.saga_sd.solInfo, 'y');
results.orig{1,1}.green.saga_sd = rmfield(results.orig{1,1}.green.saga_sd, 'r');
results.orig{1,1}.green.gendatapl = rmfield(results.orig{1,1}.green.gendatapl, 'A');
results.orig{1,1}.green.gendatapl = rmfield(results.orig{1,1}.green.gendatapl, 'b');
results.orig{1,1}.green.gendatapl.genInfo = rmfield(results.orig{1,1}.green.gendatapl.genInfo, 'b0');
results.orig{2,1}.green.saga_sd.solInfo = rmfield(results.orig{2,1}.green.saga_sd.solInfo, 'y');
results.orig{2,1}.green.saga_sd = rmfield(results.orig{2,1}.green.saga_sd, 'r');
results.orig{2,1}.green.gendatapl = rmfield(results.orig{2,1}.green.gendatapl, 'A');
results.orig{2,1}.green.gendatapl = rmfield(results.orig{2,1}.green.gendatapl, 'b');
results.orig{2,1}.green.gendatapl.genInfo = rmfield(results.orig{2,1}.green.gendatapl.genInfo, 'b0');

results.orig{1,1}.blue.saga_sd.solInfo = rmfield(results.orig{1,1}.blue.saga_sd.solInfo, 'y');
results.orig{1,1}.blue.saga_sd = rmfield(results.orig{1,1}.blue.saga_sd, 'r');
results.orig{1,1}.blue.gendatapl = rmfield(results.orig{1,1}.blue.gendatapl, 'A');
results.orig{1,1}.blue.gendatapl = rmfield(results.orig{1,1}.blue.gendatapl, 'b');
results.orig{1,1}.blue.gendatapl.genInfo = rmfield(results.orig{1,1}.blue.gendatapl.genInfo, 'b0');
results.orig{2,1}.blue.saga_sd.solInfo = rmfield(results.orig{2,1}.blue.saga_sd.solInfo, 'y');
results.orig{2,1}.blue.saga_sd = rmfield(results.orig{2,1}.blue.saga_sd, 'r');
results.orig{2,1}.blue.gendatapl = rmfield(results.orig{2,1}.blue.gendatapl, 'A');
results.orig{2,1}.blue.gendatapl = rmfield(results.orig{2,1}.blue.gendatapl, 'b');
results.orig{2,1}.blue.gendatapl.genInfo = rmfield(results.orig{2,1}.blue.gendatapl.genInfo, 'b0');


results.ada40{1,1}.red.saga_sd.solInfo = rmfield(results.ada40{1,1}.red.saga_sd.solInfo, 'y');
results.ada40{1,1}.red.saga_sd = rmfield(results.ada40{1,1}.red.saga_sd, 'r');
results.ada40{1,1}.red.gendatapl = rmfield(results.ada40{1,1}.red.gendatapl, 'A');
results.ada40{1,1}.red.gendatapl = rmfield(results.ada40{1,1}.red.gendatapl, 'b');
results.ada40{1,1}.red.gendatapl.genInfo = rmfield(results.ada40{1,1}.red.gendatapl.genInfo, 'b0');
results.ada40{2,1}.red.saga_sd.solInfo = rmfield(results.ada40{2,1}.red.saga_sd.solInfo, 'y');
results.ada40{2,1}.red.saga_sd = rmfield(results.ada40{2,1}.red.saga_sd, 'r');
results.ada40{2,1}.red.gendatapl = rmfield(results.ada40{2,1}.red.gendatapl, 'A');
results.ada40{2,1}.red.gendatapl = rmfield(results.ada40{2,1}.red.gendatapl, 'b');
results.ada40{2,1}.red.gendatapl.genInfo = rmfield(results.ada40{2,1}.red.gendatapl.genInfo, 'b0');

results.ada40{1,1}.green.saga_sd.solInfo = rmfield(results.ada40{1,1}.green.saga_sd.solInfo, 'y');
results.ada40{1,1}.green.saga_sd = rmfield(results.ada40{1,1}.green.saga_sd, 'r');
results.ada40{1,1}.green.gendatapl = rmfield(results.ada40{1,1}.green.gendatapl, 'A');
results.ada40{1,1}.green.gendatapl = rmfield(results.ada40{1,1}.green.gendatapl, 'b');
results.ada40{1,1}.green.gendatapl.genInfo = rmfield(results.ada40{1,1}.green.gendatapl.genInfo, 'b0');
results.ada40{2,1}.green.saga_sd.solInfo = rmfield(results.ada40{2,1}.green.saga_sd.solInfo, 'y');
results.ada40{2,1}.green.saga_sd = rmfield(results.ada40{2,1}.green.saga_sd, 'r');
results.ada40{2,1}.green.gendatapl = rmfield(results.ada40{2,1}.green.gendatapl, 'A');
results.ada40{2,1}.green.gendatapl = rmfield(results.ada40{2,1}.green.gendatapl, 'b');
results.ada40{2,1}.green.gendatapl.genInfo = rmfield(results.ada40{2,1}.green.gendatapl.genInfo, 'b0');

results.ada40{1,1}.blue.saga_sd.solInfo = rmfield(results.ada40{1,1}.blue.saga_sd.solInfo, 'y');
results.ada40{1,1}.blue.saga_sd = rmfield(results.ada40{1,1}.blue.saga_sd, 'r');
results.ada40{1,1}.blue.gendatapl = rmfield(results.ada40{1,1}.blue.gendatapl, 'A');
results.ada40{1,1}.blue.gendatapl = rmfield(results.ada40{1,1}.blue.gendatapl, 'b');
results.ada40{1,1}.blue.gendatapl.genInfo = rmfield(results.ada40{1,1}.blue.gendatapl.genInfo, 'b0');
results.ada40{2,1}.blue.saga_sd.solInfo = rmfield(results.ada40{2,1}.blue.saga_sd.solInfo, 'y');
results.ada40{2,1}.blue.saga_sd = rmfield(results.ada40{2,1}.blue.saga_sd, 'r');
results.ada40{2,1}.blue.gendatapl = rmfield(results.ada40{2,1}.blue.gendatapl, 'A');
results.ada40{2,1}.blue.gendatapl = rmfield(results.ada40{2,1}.blue.gendatapl, 'b');
results.ada40{2,1}.blue.gendatapl.genInfo = rmfield(results.ada40{2,1}.blue.gendatapl.genInfo, 'b0');


results.ada80{1,1}.red.saga_sd.solInfo = rmfield(results.ada80{1,1}.red.saga_sd.solInfo, 'y');
results.ada80{1,1}.red.saga_sd = rmfield(results.ada80{1,1}.red.saga_sd, 'r');
results.ada80{1,1}.red.gendatapl = rmfield(results.ada80{1,1}.red.gendatapl, 'A');
results.ada80{1,1}.red.gendatapl = rmfield(results.ada80{1,1}.red.gendatapl, 'b');
results.ada80{1,1}.red.gendatapl.genInfo = rmfield(results.ada80{1,1}.red.gendatapl.genInfo, 'b0');
results.ada80{2,1}.red.saga_sd.solInfo = rmfield(results.ada80{2,1}.red.saga_sd.solInfo, 'y');
results.ada80{2,1}.red.saga_sd = rmfield(results.ada80{2,1}.red.saga_sd, 'r');
results.ada80{2,1}.red.gendatapl = rmfield(results.ada80{2,1}.red.gendatapl, 'A');
results.ada80{2,1}.red.gendatapl = rmfield(results.ada80{2,1}.red.gendatapl, 'b');
results.ada80{2,1}.red.gendatapl.genInfo = rmfield(results.ada80{2,1}.red.gendatapl.genInfo, 'b0');

results.ada80{1,1}.green.saga_sd.solInfo = rmfield(results.ada80{1,1}.green.saga_sd.solInfo, 'y');
results.ada80{1,1}.green.saga_sd = rmfield(results.ada80{1,1}.green.saga_sd, 'r');
results.ada80{1,1}.green.gendatapl = rmfield(results.ada80{1,1}.green.gendatapl, 'A');
results.ada80{1,1}.green.gendatapl = rmfield(results.ada80{1,1}.green.gendatapl, 'b');
results.ada80{1,1}.green.gendatapl.genInfo = rmfield(results.ada80{1,1}.green.gendatapl.genInfo, 'b0');
results.ada80{2,1}.green.saga_sd.solInfo = rmfield(results.ada80{2,1}.green.saga_sd.solInfo, 'y');
results.ada80{2,1}.green.saga_sd = rmfield(results.ada80{2,1}.green.saga_sd, 'r');
results.ada80{2,1}.green.gendatapl = rmfield(results.ada80{2,1}.green.gendatapl, 'A');
results.ada80{2,1}.green.gendatapl = rmfield(results.ada80{2,1}.green.gendatapl, 'b');
results.ada80{2,1}.green.gendatapl.genInfo = rmfield(results.ada80{2,1}.green.gendatapl.genInfo, 'b0');

results.ada80{1,1}.blue.saga_sd.solInfo = rmfield(results.ada80{1,1}.blue.saga_sd.solInfo, 'y');
results.ada80{1,1}.blue.saga_sd = rmfield(results.ada80{1,1}.blue.saga_sd, 'r');
results.ada80{1,1}.blue.gendatapl = rmfield(results.ada80{1,1}.blue.gendatapl, 'A');
results.ada80{1,1}.blue.gendatapl = rmfield(results.ada80{1,1}.blue.gendatapl, 'b');
results.ada80{1,1}.blue.gendatapl.genInfo = rmfield(results.ada80{1,1}.blue.gendatapl.genInfo, 'b0');
results.ada80{2,1}.blue.saga_sd.solInfo = rmfield(results.ada80{2,1}.blue.saga_sd.solInfo, 'y');
results.ada80{2,1}.blue.saga_sd = rmfield(results.ada80{2,1}.blue.saga_sd, 'r');
results.ada80{2,1}.blue.gendatapl = rmfield(results.ada80{2,1}.blue.gendatapl, 'A');
results.ada80{2,1}.blue.gendatapl = rmfield(results.ada80{2,1}.blue.gendatapl, 'b');
results.ada80{2,1}.blue.gendatapl.genInfo = rmfield(results.ada80{2,1}.blue.gendatapl.genInfo, 'b0');

data.UC_Davis_results = results;

[~, ~, imOrig, imSaga] = experiments.exputil.viewrecoveredimageandresidual(data.UC_Davis_results.ada40{2,1});
imwrite(imOrig, strcat(folder_root, 'UCD_orig.png'));
imwrite(imSaga, strcat(folder_root, 'UCD_L_15_ada40.png'));
save(strcat(folder_root, 'data_image.mat'), 'data');
delete(strcat(folder_root, 'UC_Davis_small/', 'results.mat'))
rmdir(strcat(folder_root, 'UC_Davis_small/'))


end





function data = plot_figures_and_tables(folder_root)

%{

load(strcat(folder_root, 'data_gaussian.mat'), 'data')

mvs_n_orig = mean(data.n_exp.mvs_orig, 2);
mvs_n_ada40 = mean(data.n_exp.mvs_ada40, 2);
mvs_n_percent_decr = 1 - mvs_n_ada40./mvs_n_orig;

mvs_L_orig = mean(data.L_exp.mvs_orig, 2);
mvs_L_ada40 = mean(data.L_exp.mvs_ada40, 2);
mvs_L_percent_decr = 1 - mvs_L_ada40./mvs_L_orig;

mvs_nr_orig = mean(data.nr_exp.mvs_orig, 2);
mvs_nr_ada40 = mean(data.nr_exp.mvs_ada40, 2);
mvs_nr_percent_decr = 1 - mvs_nr_ada40./mvs_nr_orig;

figure
subplot(1,3,1)
plot(mvs_n_orig)
hold on
plot(mvs_n_ada40, '--')
xticks([1,2,3,4])
xticklabels({'128','256','512','1028'})
xlabel('Signal size n')
ylim([0, round(max(mvs_n_orig))+10000])
%yticklabels({'0', '5000', '10000', '15000', '20000'})
%curtick = get(gca, 'YTick');
%set(gca, 'YTickLabel', cellstr(num2str(curtick(:))));
ylabel('Number mat-vecs')
hold off

subplot(1,3,2)
plot(mvs_L_orig)
hold on
plot(mvs_L_ada40, '--')
xticks([1,2,3,4])
xticklabels({'8', '12', '16', '20'})
xlabel('Oversampling rate L')
ylim([0, round(max(mvs_L_orig))+10000])
%yticklabels({'0', '5000', '10000', '15000'})
%curtick = get(gca, 'YTick');
%set(gca, 'YTickLabel', cellstr(num2str(curtick(:))));
ylabel('Number mat-vecs')
hold off


subplot(1,3,3)
plot(mvs_nr_orig)
hold on
plot(mvs_nr_ada40, '--')
xticks([1,2,3,4,5,6])
xticklabels({'0.05','0.10','0.15','0.20','0.25','0.30'})
xlabel('Noise ratio')
ylim([0, round(max(mvs_nr_orig))+10000])
%yticklabels({'0', '5000', '10000', '15000'})
%curtick = get(gca, 'YTick');
%set(gca, 'YTickLabel', cellstr(num2str(curtick(:))));
ylabel('Number mat-vecs')
hold off

%}



load(strcat(folder_root, 'data_image.mat'), 'data');

mvs_jul_and_me_percent_decr_ada40 = 1 - data.jul_and_me_exp.mvs_ada40./data.jul_and_me_exp.mvs_orig;
mvs_jul_and_me_percent_decr_ada80 = 1 - data.jul_and_me_exp.mvs_ada80./data.jul_and_me_exp.mvs_orig;
[n1, n2] = size(data.jul_and_me_results.orig{2,1}.red.saga_sd.x);
n_jul_and_me = n1*n2;
mvs_UCD_percent_decr_ada40 = 1 - data.UCD_exp.mvs_ada40./data.UCD_exp.mvs_orig;
mvs_UCD_percent_decr_ada80 = 1 - data.UCD_exp.mvs_ada80./data.UCD_exp.mvs_orig;
[n1, n2] = size(data.UC_Davis_results.orig{2,1}.red.saga_sd.x);
n_UC_Davis = n1*n2;


fprintf('\n')
fprintf('Number of matrix-vector products for image experiments\n');
fprintf('                   |   Original    |  Algorithm 8  |  Algorithm 8\n')
fprintf(' image    n     L  | r = 2, m = 20 |    m = 40     |    m = 80\n')
fprintf(' Family %6i  %2i |   %8i    | %8i %1.2f | %8i %1.2f \n', n_jul_and_me, 10, data.jul_and_me_exp.mvs_orig(1,1), ...
   data.jul_and_me_exp.mvs_ada40(1,1), mvs_jul_and_me_percent_decr_ada40(1,1), ...
   data.jul_and_me_exp.mvs_ada80(1,1), mvs_jul_and_me_percent_decr_ada80(1,1));
fprintf(' Family %6i  %2i |   %8i    | %8i %1.2f | %8i %1.2f \n', n_jul_and_me, 15, data.jul_and_me_exp.mvs_orig(2,1), ...
   data.jul_and_me_exp.mvs_ada40(2,1), mvs_jul_and_me_percent_decr_ada40(2,1), ...
   data.jul_and_me_exp.mvs_ada80(2,1), mvs_jul_and_me_percent_decr_ada80(2,1));
fprintf(' UCD    %6i  %2i |   %8i    | %8i %1.2f | %8i %1.2f \n', n_UC_Davis, 10, data.UCD_exp.mvs_orig(1,1), ...
   data.UCD_exp.mvs_ada40(1,1), mvs_UCD_percent_decr_ada40(1,1), ...
   data.UCD_exp.mvs_ada80(1,1), mvs_UCD_percent_decr_ada80(1,1));
fprintf(' UCD    %6i  %2i |   %8i    | %8i %1.2f | %8i %1.2f \n', n_UC_Davis, 15, data.UCD_exp.mvs_orig(2,1), ...
   data.UCD_exp.mvs_ada40(2,1), mvs_UCD_percent_decr_ada40(2,1), ...
   data.UCD_exp.mvs_ada80(2,1), mvs_UCD_percent_decr_ada80(2,1));


figure
subplot(1,2,1);
img = imread(strcat(folder_root, 'UCD_orig.png'));
image(img);
set(gca,'visible','off')

subplot(1,2,2);
img = imread(strcat(folder_root, 'UCD_L_15_ada40.png'));
image(img);
set(gca,'visible','off')

figure
subplot(1,2,1);
img = imread(strcat(folder_root, 'jul_and_me_orig.png'));
image(img);
set(gca,'visible','off')

subplot(1,2,2);
img = imread(strcat(folder_root, 'jul_and_me_L_15_ada40.png'));
image(img);
set(gca,'visible','off')

end

