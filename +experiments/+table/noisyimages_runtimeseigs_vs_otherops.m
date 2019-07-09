function [data] = noisyimages_runtimeseigs_vs_otherops()
% experiments.table.noisyimages_runtimeseigs_vs_otherops
%
% Generates table of runtime results for recovery of noisy natural image
% with resizing to 4096 pixels and 16384 pixels.

resizeImage = 1/1;
iterations = 100;
noise_ratio = 0.15;
L_range = [6, 9, 12];
image_file_range = {'data/parrot_4k.jpg', 'data/parrot_16k.jpg'};
noise_type = 'gaussian';
max_num_tests = 10;
TeX_print_mode_on = false;


% Loads cached (or solves and caches) results of test problem
folder_name = 'cache/table_noisyimages_runtimeseigs_vs_otherops';
if exist(folder_name) ~= 7
   mkdir(folder_name);
end
exists_default_experiment = exist(strcat(folder_name, '/data.mat'));
if exists_default_experiment
	load(strcat(folder_name, '/data.mat'));
else
	data = [];
end


if isempty(data)
   for image_idx = 1:2
      image_file = image_file_range{image_idx};
      for L_idx = 1:3
         L = L_range(L_idx);
         for test_num = 1:max_num_tests
            results = experiments.evolvingmatricespl('L', L, 'seed', test_num, ...
                  'image_file', image_file, ...
                  'noise_type', noise_type, ...
                  'resizeImage', resizeImage, ...
                  'noise_ratio', noise_ratio, 'iterations', iterations, ...
                  'verbosity', 1, 'logIterData', false, ...
                  'printRelErr', true, ...
                  'skipDFP', true, 'StopOnFeasible', false, ...
                  'returnEigsIterData', true, 'optTol', 0, ...
                  'eigs_basis_size', [], 'cutDim', 2);
            solInfo = results.saga_sd.solInfo;
            runtimeTotal(image_idx, L_idx, test_num) ...
               = solInfo.main.tMain + solInfo.main.tObjective + solInfo.main.tPrimalRecovery + ...
               + solInfo.main.tDualRecovery + solInfo.main.tLinesearch;
            runtimeEigs(image_idx, L_idx, test_num) = solInfo.hop.tEigs;
            runtimePrimalRec(image_idx, L_idx, test_num) = solInfo.main.tPrimalRecovery;
            runtimeOtherOps(image_idx, L_idx, test_num) = runtimeTotal(image_idx, L_idx, test_num) ...
               - (runtimeEigs(image_idx, L_idx, test_num) + runtimePrimalRec(image_idx, L_idx, test_num));

            numDFTsTotal(image_idx, L_idx, test_num) = solInfo.nfft;
            numDFTsEigs(image_idx, L_idx, test_num) = solInfo.hop.nFFTsInEigs;
            numDFTsPrimalRec(image_idx, L_idx, test_num) = solInfo.hop.nPFDObjective*2*L;
            numDFTsOtherOps(image_idx, L_idx, test_num) = solInfo.nfft ...
               - (numDFTsEigs(image_idx, L_idx, test_num) + numDFTsPrimalRec(image_idx, L_idx, test_num));
         end   
         data.runtimeTotalMean(image_idx, L_idx) = mean(runtimeTotal(image_idx, L_idx, :));
         data.runtimeEigsMean(image_idx, L_idx) = mean(runtimeEigs(image_idx, L_idx, :));
         data.runtimePrimalRecMean(image_idx, L_idx) = mean(runtimePrimalRec(image_idx, L_idx, :));
         data.runtimeOtherOpsMean(image_idx, L_idx) = mean(runtimeOtherOps(image_idx, L_idx, :));

         data.numDFTsTotalMean(image_idx, L_idx) = mean(numDFTsTotal(image_idx, L_idx, :));
         data.numDFTsEigsMean(image_idx, L_idx) = mean(numDFTsEigs(image_idx, L_idx, :));
         data.numDFTsPrimalRecMean(image_idx, L_idx) = mean(numDFTsPrimalRec(image_idx, L_idx, :));
         data.numDFTsOtherOpsMean(image_idx, L_idx) = mean(numDFTsOtherOps(image_idx, L_idx, :));
      end
   end
   save(strcat(folder_name, '/data.mat'), 'data');
end


if TeX_print_mode_on
   fprintf('\n')
   fprintf('Runtimes and number of DFTs for various models\n');
   fprintf('           |             EMEP           |      Primal refinement    |      All other ops\n')
   fprintf('    n   L  |     min          DFTs      |      min          DFTs    |     min         DFTs  \n')
   for L_idx = 1:3
      image_idx = 1;
      L = L_range(L_idx);
      fprintf('%5i & %2i &', 4096, L)
      fprintf(' %5.2f  (%1.2f) &', data.runtimeEigsMean(image_idx, L_idx)/60, data.runtimeEigsMean(image_idx, L_idx)/data.runtimeTotalMean(image_idx, L_idx));
      fprintf(' %6.0f  (%1.2f) &', data.numDFTsEigsMean(image_idx, L_idx)/60, data.numDFTsEigsMean(image_idx, L_idx)/data.numDFTsTotalMean(image_idx, L_idx));
      fprintf(' %4.2f  (%1.2f) &', data.runtimePrimalRecMean(image_idx, L_idx)/60, data.runtimePrimalRecMean(image_idx, L_idx)/data.runtimeTotalMean(image_idx, L_idx));
      fprintf(' %6.0f  (%1.2f) &', data.numDFTsPrimalRecMean(image_idx, L_idx)/60, data.numDFTsPrimalRecMean(image_idx, L_idx)/data.numDFTsTotalMean(image_idx, L_idx));
      fprintf(' %4.2f &', data.runtimeOtherOpsMean(image_idx, L_idx)/60);
      fprintf(' %6.0f ', data.numDFTsOtherOpsMean(image_idx, L_idx)/60);   
      fprintf('\n')
   end
   for L_idx = 1:3
      image_idx = 2;
      L = L_range(L_idx);
      fprintf('%5i & %2i &', 16384, L)
      fprintf(' %4.2f  (%1.2f) &', data.runtimeEigsMean(image_idx, L_idx)/60, data.runtimeEigsMean(image_idx, L_idx)/data.runtimeTotalMean(image_idx, L_idx));
      fprintf(' %6.0f  (%1.2f) &', data.numDFTsEigsMean(image_idx, L_idx)/60, data.numDFTsEigsMean(image_idx, L_idx)/data.numDFTsTotalMean(image_idx, L_idx));
      fprintf(' %4.2f  (%1.2f) &', data.runtimePrimalRecMean(image_idx, L_idx)/60, data.runtimePrimalRecMean(image_idx, L_idx)/data.runtimeTotalMean(image_idx, L_idx));
      fprintf(' %6.0f  (%1.2f) &', data.numDFTsPrimalRecMean(image_idx, L_idx)/60, data.numDFTsPrimalRecMean(image_idx, L_idx)/data.numDFTsTotalMean(image_idx, L_idx));
      fprintf(' %4.2f &', data.runtimeOtherOpsMean(image_idx, L_idx)/60);
      fprintf(' %6.0f', data.numDFTsOtherOpsMean(image_idx, L_idx)/60);   
      fprintf('\n')
   end
else
   fprintf('\n')
   fprintf('Runtimes and number of DFTs for various models\n');
   fprintf('           |             EMEP           |      Primal refinement    |      All other ops\n')
   fprintf('    n   L  |     min          DFTs      |      min          DFTs    |     min         DFTs  \n')
   for L_idx = 1:3
      image_idx = 1;
      L = L_range(L_idx);
      fprintf('%5i  %2i  |', 4096, L)
      fprintf(' %5.2f (%1.2f)', data.runtimeEigsMean(image_idx, L_idx)/60, data.runtimeEigsMean(image_idx, L_idx)/data.runtimeTotalMean(image_idx, L_idx));
      fprintf(' %6.0f (%1.2f) |', data.numDFTsEigsMean(image_idx, L_idx)/60, data.numDFTsEigsMean(image_idx, L_idx)/data.numDFTsTotalMean(image_idx, L_idx));
      fprintf(' %4.2f (%1.2f)', data.runtimePrimalRecMean(image_idx, L_idx)/60, data.runtimePrimalRecMean(image_idx, L_idx)/data.runtimeTotalMean(image_idx, L_idx));
      fprintf(' %6.0f (%1.2f) |', data.numDFTsPrimalRecMean(image_idx, L_idx)/60, data.numDFTsPrimalRecMean(image_idx, L_idx)/data.numDFTsTotalMean(image_idx, L_idx));
      fprintf(' %4.2f (%1.2f)', data.runtimeOtherOpsMean(image_idx, L_idx)/60, data.runtimeOtherOpsMean(image_idx, L_idx)/data.runtimeTotalMean(image_idx, L_idx));
      fprintf(' %6.0f (%1.2f) |', data.numDFTsOtherOpsMean(image_idx, L_idx)/60, data.numDFTsOtherOpsMean(image_idx, L_idx)/data.numDFTsTotalMean(image_idx, L_idx));   
      fprintf('\n')
   end
   for L_idx = 1:3
      image_idx = 2;
      L = L_range(L_idx);
      fprintf('%5i  %2i  |', 16384, L)
      fprintf(' %4.2f (%1.2f)', data.runtimeEigsMean(image_idx, L_idx)/60, data.runtimeEigsMean(image_idx, L_idx)/data.runtimeTotalMean(image_idx, L_idx));
      fprintf(' %6.0f (%1.2f) |', data.numDFTsEigsMean(image_idx, L_idx)/60, data.numDFTsEigsMean(image_idx, L_idx)/data.numDFTsTotalMean(image_idx, L_idx));
      fprintf(' %4.2f (%1.2f)', data.runtimePrimalRecMean(image_idx, L_idx)/60, data.runtimePrimalRecMean(image_idx, L_idx)/data.runtimeTotalMean(image_idx, L_idx));
      fprintf(' %6.0f (%1.2f) |', data.numDFTsPrimalRecMean(image_idx, L_idx)/60, data.numDFTsPrimalRecMean(image_idx, L_idx)/data.numDFTsTotalMean(image_idx, L_idx));
      fprintf(' %4.2f (%1.2f)', data.runtimeOtherOpsMean(image_idx, L_idx)/60, data.runtimeOtherOpsMean(image_idx, L_idx)/data.runtimeTotalMean(image_idx, L_idx));
      fprintf(' %6.0f (%1.2f) |', data.numDFTsOtherOpsMean(image_idx, L_idx)/60, data.numDFTsOtherOpsMean(image_idx, L_idx)/data.numDFTsTotalMean(image_idx, L_idx));   
      fprintf('\n')
   end
end


end