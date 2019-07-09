function [data] = noisyimage_costs()
% experiments.figure.noisyimage_costs
%
% This script determines if synthetic and natural noise 
% create different behaviors for eig problem in saga_sd

L_range = [5, 10];
noise_ratio_range = [0.05, 0.15, 0.30];
im_file_range = ["data/parrot_4k.jpg", "data/parrot_16k.jpg"];
max_num_tests = 1;
TeX_print_mode_on = true;

folder_name = 'cache/figure_noisyimage_costs';
if exist(folder_name) ~= 7
   mkdir(folder_name);
end
exists_default_experiment = exist(strcat(folder_name, '/data.mat'));
if exists_default_experiment
   load(strcat(folder_name, '/data.mat'));
else
   data = cell(length(L_range), length(noise_ratio_range), length(im_file_range), max_num_tests);
end



for i = 1:length(L_range)
   for j = 1:length(noise_ratio_range)
      for k = 1:length(im_file_range)
         for test_num = 1:max_num_tests
            if isempty(data{i,j,k,test_num})
               L = L_range(i);
               noise_ratio = noise_ratio_range(j);
               image_file = im_file_range(k);
%[i,j,k,test_num]
               evol_mats_opts = struct('noise_ratio', noise_ratio, 'L', L, ...
                  'image_file', image_file{1}, 'logIterData', false, 'noise_type', 'gaussian', ...
                  'seed', test_num, 'maxTime', 24*60*60);
               results = experiments.evolvingmatricespl(evol_mats_opts);
               data{i,j,k,test_num} = get_data(results, L, noise_ratio, image_file);

               save(strcat(folder_name, '/data.mat'), 'data');
            end
         end
      end
   end
end


runtimeTotalMean = zeros(length(L_range), length(noise_ratio_range), length(im_file_range));
runtimeEigsMean = zeros(length(L_range), length(noise_ratio_range), length(im_file_range));
runtimePrimalRecMean = zeros(length(L_range), length(noise_ratio_range), length(im_file_range));
runtimeOtherOpsMean = zeros(length(L_range), length(noise_ratio_range), length(im_file_range));

numDFTsTotalMean = zeros(length(L_range), length(noise_ratio_range), length(im_file_range));
numDFTsEigsMean = zeros(length(L_range), length(noise_ratio_range), length(im_file_range));
numDFTsPrimalRecMean = zeros(length(L_range), length(noise_ratio_range), length(im_file_range));
numDFTsOtherOpsMean = zeros(length(L_range), length(noise_ratio_range), length(im_file_range));

for i = 1:length(L_range)
   for j = 1:length(noise_ratio_range)
      for k = 1:length(im_file_range)
         for test_num = 1:max_num_tests
            runtimeTotalMean(i, j, k) = runtimeTotalMean(i, j, k) + data{i, j, k, test_num}.runtimeTotal/max_num_tests;
            runtimeEigsMean(i, j, k) = runtimeEigsMean(i, j, k) + data{i, j, k, test_num}.runtimeEigs/max_num_tests;
            runtimePrimalRecMean(i, j, k) = runtimePrimalRecMean(i, j, k) + data{i, j, k, test_num}.runtimePrimalRec/max_num_tests;
            runtimeOtherOpsMean(i, j, k) = runtimeOtherOpsMean(i, j, k) + data{i, j, k, test_num}.runtimeOtherOps/max_num_tests;

            numDFTsTotalMean(i, j, k) = numDFTsTotalMean(i, j, k) + data{i, j, k, test_num}.numDFTsTotal/max_num_tests;
            numDFTsEigsMean(i, j, k) = numDFTsEigsMean(i, j, k) + data{i, j, k, test_num}.numDFTsEigs/max_num_tests;
            numDFTsPrimalRecMean(i, j, k) = numDFTsPrimalRecMean(i, j, k) + data{i, j, k, test_num}.numDFTsPrimalRec/max_num_tests;
            numDFTsOtherOpsMean(i, j, k) = numDFTsOtherOpsMean(i, j, k) + data{i, j, k, test_num}.numDFTsOtherOps/max_num_tests;
         end
      end
   end
end


if TeX_print_mode_on
   fprintf('\n')
   fprintf('Runtimes and number of DFTs for various models\n');
   fprintf('           noise| Num |             EMEP           |      Primal refinement    |      All other ops\n')
   fprintf('    n   L  ratio| eigs|     min          DFTs      |      min          DFTs    |     min         DFTs  \n')
   for i = 1:length(L_range)
      k = 1;
      L = L_range(i);
      for j = 1:length(noise_ratio_range)
         noise_ratio = noise_ratio_range(j);
         fprintf('%5i & %2i & %1.2f & %3i &', 4096, L, noise_ratio, length(data{i, j, k, 1}.nMatVec))
         fprintf(' %4.2f  (%1.2f) &', runtimeEigsMean(i, j, k)/60, runtimeEigsMean(i, j, k)/runtimeTotalMean(i, j, k));
         fprintf(' %6.0f  (%1.2f) &', numDFTsEigsMean(i, j, k)/60, numDFTsEigsMean(i, j, k)/numDFTsTotalMean(i, j, k));
         fprintf(' %4.2f  (%1.2f) &', runtimePrimalRecMean(i, j, k)/60, runtimePrimalRecMean(i, j, k)/runtimeTotalMean(i, j, k));
         fprintf(' %6.0f  (%1.2f) &', numDFTsPrimalRecMean(i, j, k)/60, numDFTsPrimalRecMean(i, j, k)/numDFTsTotalMean(i, j, k));
         fprintf(' %4.2f &', runtimeOtherOpsMean(i, j, k)/60);
         fprintf(' %6.0f ', numDFTsOtherOpsMean(i, j, k)/60);   
         fprintf('\n')
      end
   end
   for i = 1:length(L_range)
      k = 2;
      L = L_range(i);
      for j = 1:length(noise_ratio_range)
         noise_ratio = noise_ratio_range(j);
         fprintf('%5i & %2i & %1.2f & %3i &', 16384, L, noise_ratio, length(data{i, j, k, 1}.nMatVec))
         fprintf(' %4.2f  (%1.2f) &', runtimeEigsMean(i, j, k)/60, runtimeEigsMean(i, j, k)/runtimeTotalMean(i, j, k));
         fprintf(' %6.0f  (%1.2f) &', numDFTsEigsMean(i, j, k)/60, numDFTsEigsMean(i, j, k)/numDFTsTotalMean(i, j, k));
         fprintf(' %4.2f  (%1.2f) &', runtimePrimalRecMean(i, j, k)/60, runtimePrimalRecMean(i, j, k)/runtimeTotalMean(i, j, k));
         fprintf(' %6.0f  (%1.2f) &', numDFTsPrimalRecMean(i, j, k)/60, numDFTsPrimalRecMean(i, j, k)/numDFTsTotalMean(i, j, k));
         fprintf(' %4.2f &', runtimeOtherOpsMean(i, j, k)/60);
         fprintf(' %6.0f', numDFTsOtherOpsMean(i, j, k)/60);   
         fprintf('\n')
      end
   end
else
   fprintf('\n')
   fprintf('Runtimes and number of DFTs for various models\n');
   fprintf('           noise| Num |             EMEP           |      Primal refinement    |      All other ops\n')
   fprintf('    n   L  ratio| eigs|      min          DFTs     |      min          DFTs    |     min           DFTs  \n')
   for i = 1:length(L_range)
      k = 1;
      L = L_range(i);
      for j = 1:length(noise_ratio_range)
         noise_ratio = noise_ratio_range(j);
         fprintf('%5i  %2i  %1.2f | %3i |', 4096, L, noise_ratio, length(data{i, j, k, 1}.nMatVec))
         fprintf(' %5.2f (%1.2f)', runtimeEigsMean(i, j, k)/60, runtimeEigsMean(i, j, k)/runtimeTotalMean(i, j, k));
         fprintf(' %6.0f (%1.2f) |', numDFTsEigsMean(i, j, k)/60, numDFTsEigsMean(i, j, k)/numDFTsTotalMean(i, j, k));
         fprintf(' %4.2f (%1.2f)', runtimePrimalRecMean(i, j, k)/60, runtimePrimalRecMean(i, j, k)/runtimeTotalMean(i, j, k));
         fprintf(' %6.0f (%1.2f) |', numDFTsPrimalRecMean(i, j, k)/60, numDFTsPrimalRecMean(i, j, k)/numDFTsTotalMean(i, j, k));
         fprintf(' %4.2f (%1.2f)', runtimeOtherOpsMean(i, j, k)/60, runtimeOtherOpsMean(i, j, k)/runtimeTotalMean(i, j, k));
         fprintf(' %6.0f (%1.2f) |', numDFTsOtherOpsMean(i, j, k)/60, numDFTsOtherOpsMean(i, j, k)/numDFTsTotalMean(i, j, k));   
         fprintf('\n')
      end
   end
   for i = 1:length(L_range)
      k = 2;
      L = L_range(i);
      for j = 1:length(noise_ratio_range)
         noise_ratio = noise_ratio_range(j);
         fprintf('%5i  %2i  %1.2f | %3i |', 16384, L, noise_ratio, length(data{i, j, k, 1}.nMatVec))
         fprintf(' %5.2f (%1.2f)', runtimeEigsMean(i, j, k)/60, runtimeEigsMean(i, j, k)/runtimeTotalMean(i, j, k));
         fprintf(' %6.0f (%1.2f) |', numDFTsEigsMean(i, j, k)/60, numDFTsEigsMean(i, j, k)/numDFTsTotalMean(i, j, k));
         fprintf(' %4.2f (%1.2f)', runtimePrimalRecMean(i, j, k)/60, runtimePrimalRecMean(i, j, k)/runtimeTotalMean(i, j, k));
         fprintf(' %6.0f (%1.2f) |', numDFTsPrimalRecMean(i, j, k)/60, numDFTsPrimalRecMean(i, j, k)/numDFTsTotalMean(i, j, k));
         fprintf(' %4.2f (%1.2f)', runtimeOtherOpsMean(i, j, k)/60, runtimeOtherOpsMean(i, j, k)/runtimeTotalMean(i, j, k));
         fprintf(' %6.0f (%1.2f) |', numDFTsOtherOpsMean(i, j, k)/60, numDFTsOtherOpsMean(i, j, k)/numDFTsTotalMean(i, j, k));   
         fprintf('\n')
      end
   end
end
%}

figure;
subplot_idx = 1;
k = 1;
for L_idx = 1:length(L_range)
   for noise_ratio_idx = 1:length(noise_ratio_range)
      L = L_range(L_idx);
      noise_ratio = noise_ratio_range(noise_ratio_idx);
      subplot(2, 3, subplot_idx)
      plot(data{L_idx,noise_ratio_idx,k,1}.nMatVec, 'r*')
      xlim([1, length(data{L_idx,noise_ratio_idx,k,1}.nMatVec)])
      title_str = strcat("L = ", num2str(L), ", noise ratio = ", num2str(noise_ratio));
      title(title_str)
      xlabel('EMEP iteration')
      ylabel('Number mat-vecs')
      subplot_idx = subplot_idx + 1;
   end
end






end




function [data] = get_data(results, L, noise_ratio, image_file)

data = struct;
nMatVec = [];
nMatVec(1,1) = results.saga_sd.solInfo.iterData.infoEigsFFTs_y{1,1}.nMatVec;
for i = 1:length(results.saga_sd.solInfo.iterData.infoEigsFFTs_yLS)
   for j = 1:length(results.saga_sd.solInfo.iterData.infoEigsFFTs_yLS{1,i})
      nMatVec = [nMatVec; results.saga_sd.solInfo.iterData.infoEigsFFTs_yLS{1,i}{j}.nMatVec];
   end
end
data.nMatVec = nMatVec;
data.L = L;
data.noise_ratio = noise_ratio;
data.image_file = image_file;

data.runtimeTotal = results.saga_sd.solInfo.main.tMain ...
   + results.saga_sd.solInfo.main.tObjective + results.saga_sd.solInfo.main.tPrimalRecovery + ...
   + results.saga_sd.solInfo.main.tDualRecovery + results.saga_sd.solInfo.main.tLinesearch;
data.runtimeEigs = results.saga_sd.solInfo.hop.tEigs;
data.runtimePrimalRec = results.saga_sd.solInfo.main.tPrimalRecovery;
data.runtimeOtherOps = data.runtimeTotal - (data.runtimeEigs + data.runtimePrimalRec);

data.numDFTsTotal = results.saga_sd.solInfo.nfft;
data.numDFTsEigs = results.saga_sd.solInfo.hop.nFFTsInEigs;
data.numDFTsPrimalRec = results.saga_sd.solInfo.hop.nPFDObjective*2*L;
data.numDFTsOtherOps = results.saga_sd.solInfo.nfft - (data.numDFTsEigs + data.numDFTsPrimalRec);

end
