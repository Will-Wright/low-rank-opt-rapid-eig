function [results] = noisyimage_adaptive_eig_demo
% experiments.figure.noisyimage_adaptive_eig_demo



noise_ratio = 0.3;
L = 8;
optsNoisyModel = struct('wflowOn', false, ... 
                        'L', L, 'noise_ratio', noise_ratio, ...
                        'noise_type', 'gaussian', ...
                        'resizeImage', 1/1, ...
                        'optTol', 0, ...
                        'feaTol', 0, ...
                        'logIterData', false, ...
                        'StopOnRelErr', false, ...
                        'cutDim', 2);
folder_name = 'cache/figure_noisyimage_adaptive_eig_demo';
if exist(folder_name) ~= 7
   mkdir(folder_name);
end
exists_default_experiment = exist(strcat(folder_name, '/results.mat'));
if exists_default_experiment
   load(strcat(folder_name, '/results.mat'));
else
   results = cell(0);
   fprintf('Solving model with %f noise ratio.\n', noise_ratio);
   results.adaptive_4k = experiments.evolvingmatricespl(optsNoisyModel, 'useAdaptiveEigs', true, 'eigs_basis_size', 40, 'image_file', 'data/parrot_4k.jpg');
   results.orig_4k = experiments.evolvingmatricespl(optsNoisyModel, 'image_file', 'data/parrot_4k.jpg');
   results.adaptive_16k = experiments.evolvingmatricespl(optsNoisyModel, 'useAdaptiveEigs', true, 'eigs_basis_size', 40, 'image_file', 'data/parrot_16k.jpg');
   results.orig_16k = experiments.evolvingmatricespl(optsNoisyModel, 'image_file', 'data/parrot_16k.jpg');
   save(strcat(folder_name, '/results.mat'), 'results');
end
   

nMatVec_ada_4k = [];
nMatVec_ada_4k(1,1) = results.adaptive_4k.saga_sd.solInfo.iterData.infoEigsFFTs_y{1,1}.nMatVec;
for i = 1:length(results.adaptive_4k.saga_sd.solInfo.iterData.infoEigsFFTs_yLS)
   if ~isempty(results.adaptive_4k.saga_sd.solInfo.iterData.infoEigsFFTs_yLS{1,i})
      for j = 1:length(results.adaptive_4k.saga_sd.solInfo.iterData.infoEigsFFTs_yLS{1,i})
         nMatVec_ada_4k = [nMatVec_ada_4k; results.adaptive_4k.saga_sd.solInfo.iterData.infoEigsFFTs_yLS{1,i}{j}.nMatVec];
      end
   end
end

nMatVec_orig_4k = [];
nMatVec_orig_4k(1,1) = results.orig_4k.saga_sd.solInfo.iterData.infoEigsFFTs_y{1,1}.nMatVec;
for i = 1:length(results.orig_4k.saga_sd.solInfo.iterData.infoEigsFFTs_yLS)
   if ~isempty(results.orig_4k.saga_sd.solInfo.iterData.infoEigsFFTs_yLS{1,i})
      for j = 1:length(results.orig_4k.saga_sd.solInfo.iterData.infoEigsFFTs_yLS{1,i})
         nMatVec_orig_4k = [nMatVec_orig_4k; results.orig_4k.saga_sd.solInfo.iterData.infoEigsFFTs_yLS{1,i}{j}.nMatVec];
      end
   end
end

nMatVec_ada_16k = [];
nMatVec_ada_16k(1,1) = results.adaptive_16k.saga_sd.solInfo.iterData.infoEigsFFTs_y{1,1}.nMatVec;
for i = 1:length(results.adaptive_16k.saga_sd.solInfo.iterData.infoEigsFFTs_yLS)
   if ~isempty(results.adaptive_16k.saga_sd.solInfo.iterData.infoEigsFFTs_yLS{1,i})
      for j = 1:length(results.adaptive_16k.saga_sd.solInfo.iterData.infoEigsFFTs_yLS{1,i})
         nMatVec_ada_16k = [nMatVec_ada_16k; results.adaptive_16k.saga_sd.solInfo.iterData.infoEigsFFTs_yLS{1,i}{j}.nMatVec];
      end
   end
end

nMatVec_orig_16k = [];
nMatVec_orig_16k(1,1) = results.orig_16k.saga_sd.solInfo.iterData.infoEigsFFTs_y{1,1}.nMatVec;
for i = 1:length(results.orig_16k.saga_sd.solInfo.iterData.infoEigsFFTs_yLS)
   if ~isempty(results.orig_16k.saga_sd.solInfo.iterData.infoEigsFFTs_yLS{1,i})
      for j = 1:length(results.orig_16k.saga_sd.solInfo.iterData.infoEigsFFTs_yLS{1,i})
         nMatVec_orig_16k = [nMatVec_orig_16k; results.orig_16k.saga_sd.solInfo.iterData.infoEigsFFTs_yLS{1,i}{j}.nMatVec];
      end
   end
end

nMatVec_ada_4k = nMatVec_ada_4k(1:100);
nMatVec_orig_4k = nMatVec_orig_4k(1:100);
nMatVec_ada_16k = nMatVec_ada_16k(1:100);
nMatVec_orig_16k = nMatVec_orig_16k(1:100);


figure;
plot(nMatVec_ada_4k, 'o');
hold on
plot(nMatVec_orig_4k, '*');

figure;
plot(nMatVec_ada_16k, 'o');
hold on
plot(nMatVec_orig_16k, '*');

fprintf('\n   Total Mat-Vec Mults\n')
fprintf('          Default   Adaptive eigs  Decrease   \n')
fprintf('  4k    %9i       %9i     %1.3f\n', sum(nMatVec_orig_4k), sum(nMatVec_ada_4k), 1-sum(nMatVec_ada_4k)/sum(nMatVec_orig_4k));
fprintf(' 16k    %9i       %9i     %1.3f\n', sum(nMatVec_orig_16k), sum(nMatVec_ada_16k), 1-sum(nMatVec_ada_16k)/sum(nMatVec_orig_16k));
end









