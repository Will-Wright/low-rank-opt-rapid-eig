function [] = noiselessrandom_refinement_comparison()
% experiments.table.noiselessrandom_refinement_comparison
%
% Generates table with number of DFTs and runtimes for recovery of 
% noiseless random gaussian signal with and without refinement subroutines.
% DFP is dual from primal (dual refinement) and PFD is primal from dual


% Loads cached (or solves and caches) results of test problem
pathName = 'cache/table_noiselessrandom_refinement_comparison/';
if exist(pathName) ~= 7
   mkdir(pathName);
end


dataFile = fullfile(pathName, 'data_noiseless_10_128.mat');
if exist(dataFile, 'file') ~= 2

   % PFD and DFP on
   matFile = fullfile(pathName, 'results_saga_noiseless_10_128.mat');
   if exist(matFile, 'file') ~= 2
      [results_saga_noiseless_10_128] ...
         = experiments.evolvingmatricespl('L', 10, 'n', 128, ...
                                          'signal', 'gaussian', 'scale', [], ...
                                          'noise_ratio', 0.00, 'iterations', 1000, ...
                                          'verbosity', 1, 'logIterData', true, ...
                                          'printRelErr', true, ...
                                          'skipPFD', false, 'skipDFP', false, 'StopOnFeasible', false, ...
                                          'returnEigsIterData', true, ...
                                          'feaTol', 1e-5, 'optTol', 1e-5, ...
                                          'StopOnRelErr', false);
      save(matFile, 'results_saga_noiseless_10_128');
   else
      load(matFile, 'results_saga_noiseless_10_128');
   end

   % PFD on, DFP off
   matFile = fullfile(pathName, 'results_saga_noiseless_10_128_noDFP.mat');
   if exist(matFile, 'file') ~= 2
      [results_saga_noiseless_10_128_noDFP] ...
         = experiments.evolvingmatricespl('L', 10, 'n', 128, ...
                                             'signal', 'gaussian', 'scale', [], ...
                                             'noise_ratio', 0.00, 'iterations', 1000, ...
                                             'verbosity', 1, 'logIterData', true, ...
                                             'printRelErr', true, ...
                                             'skipPFD', false, 'skipDFP', true, 'StopOnFeasible', false, ...
                                             'returnEigsIterData', true, ...
                                             'feaTol', 1e-5, 'optTol', 1e-5, ...
                                             'StopOnRelErr', false);
      save(matFile, 'results_saga_noiseless_10_128_noDFP');
   else
      load(matFile, 'results_saga_noiseless_10_128_noDFP');
   end

   % PFD off, DFP on
   matFile = fullfile(pathName, 'results_saga_noiseless_10_128_noPFD.mat');
   if exist(matFile, 'file') ~= 2
      [results_saga_noiseless_10_128_noPFD] ...
         = experiments.evolvingmatricespl('L', 10, 'n', 128, ...
                                             'signal', 'gaussian', 'scale', [], ...
                                             'noise_ratio', 0.00, 'iterations', 1000, ...
                                             'verbosity', 1, 'logIterData', true, ...
                                             'printRelErr', true, ...
                                             'skipPFD', true, 'skipDFP', false, 'StopOnFeasible', false, ...
                                             'returnEigsIterData', true, ...
                                             'feaTol', 1e-5, 'optTol', 1e-5, ...
                                             'StopOnRelErr', false);
      save(matFile, 'results_saga_noiseless_10_128_noPFD');
   else
      load(matFile, 'results_saga_noiseless_10_128_noPFD');
   end

   % both PFD and DFP off
   matFile = fullfile(pathName, 'results_saga_noiseless_10_128_noPFD_noDFP.mat');
   if exist(matFile, 'file') ~= 2
      [results_saga_noiseless_10_128_noPFD_noDFP] ...
         = experiments.evolvingmatricespl('L', 10, 'n', 128, ...
                                          'signal', 'gaussian', 'scale', [], ...
                                          'noise_ratio', 0.00, 'iterations', 1000, ...
                                          'verbosity', 1, 'logIterData', true, ...
                                          'printRelErr', true, ...
                                          'skipPFD', true, 'skipDFP', true, 'StopOnFeasible', false, ...
                                          'returnEigsIterData', true, ...
                                          'feaTol', 1e-5, 'optTol', 1e-5, ...
                                          'StopOnRelErr', false);
      save(matFile, 'results_saga_noiseless_10_128_noPFD_noDFP');
   else
      load(matFile, 'results_saga_noiseless_10_128_noPFD_noDFP');
   end

end



% Loads cached data (or caches data and deletes unnecessary data)
matFile = fullfile(pathName, 'data_noiseless_10_128.mat');
if exist(matFile, 'file') ~= 2
   data_noiseless_10_128 = struct;
   
   data_noiseless_10_128.nit = results_saga_noiseless_10_128.saga_sd.solInfo.nit;
   data_noiseless_10_128.nit_noDFP = results_saga_noiseless_10_128_noDFP.saga_sd.solInfo.nit;
   data_noiseless_10_128.nit_noPFD = results_saga_noiseless_10_128_noPFD.saga_sd.solInfo.nit;
   data_noiseless_10_128.nit_noPFD_noDFP = results_saga_noiseless_10_128_noPFD_noDFP.saga_sd.solInfo.nit;

   data_noiseless_10_128.nfft = results_saga_noiseless_10_128.saga_sd.solInfo.nfft;
   data_noiseless_10_128.nfft_noDFP = results_saga_noiseless_10_128_noDFP.saga_sd.solInfo.nfft;
   data_noiseless_10_128.nfft_noPFD = results_saga_noiseless_10_128_noPFD.saga_sd.solInfo.nfft;
   data_noiseless_10_128.nfft_noPFD_noDFP = results_saga_noiseless_10_128_noPFD_noDFP.saga_sd.solInfo.nfft;   
   
   data_noiseless_10_128.time = results_saga_noiseless_10_128.saga_sd.solInfo.time;
   data_noiseless_10_128.time_noDFP = results_saga_noiseless_10_128_noDFP.saga_sd.solInfo.time;
   data_noiseless_10_128.time_noPFD = results_saga_noiseless_10_128_noPFD.saga_sd.solInfo.time;
   data_noiseless_10_128.time_noPFD_noDFP = results_saga_noiseless_10_128_noPFD_noDFP.saga_sd.solInfo.time;   
   
   x0 = results_saga_noiseless_10_128.gendatapl.x0;
   % x = results_saga_noiseless_10_128.saga_sd.x;
   x = results_saga_noiseless_10_128.saga_sd.solInfo.iterData.x{1,end};
   data_noiseless_10_128.rel_err = util.hermitianerror(x(:), x0(:), 'fro') / norm(x0(:))^2;

   x0 = results_saga_noiseless_10_128_noDFP.gendatapl.x0;
   % x = results_saga_noiseless_10_128_noDFP.saga_sd.x;
   x = results_saga_noiseless_10_128_noDFP.saga_sd.solInfo.iterData.x{1,end};
   data_noiseless_10_128.rel_err_noDFP = util.hermitianerror(x(:), x0(:), 'fro') / norm(x0(:))^2;

   x0 = results_saga_noiseless_10_128_noPFD.gendatapl.x0;
   % x = results_saga_noiseless_10_128_noPFD.saga_sd.x;
   x = results_saga_noiseless_10_128_noPFD.saga_sd.solInfo.iterData.x{1,end-1};
   data_noiseless_10_128.rel_err_noPFD = util.hermitianerror(x(:), x0(:), 'fro') / norm(x0(:))^2;

   x0 = results_saga_noiseless_10_128_noPFD_noDFP.gendatapl.x0;
   % x = results_saga_noiseless_10_128_noPFD_noDFP.saga_sd.x;
   x = results_saga_noiseless_10_128_noPFD_noDFP.saga_sd.solInfo.iterData.x{1,end-1};
   data_noiseless_10_128.rel_err_noPFD_noDFP = util.hermitianerror(x(:), x0(:), 'fro') / norm(x0(:))^2;
   save(matFile, 'data_noiseless_10_128');
   delete(fullfile(pathName, 'results_saga_noiseless_10_128.mat'));
   delete(fullfile(pathName, 'results_saga_noiseless_10_128_noDFP.mat'));
   delete(fullfile(pathName, 'results_saga_noiseless_10_128_noPFD.mat'));
   delete(fullfile(pathName, 'results_saga_noiseless_10_128_noPFD_noDFP.mat'));
else
   load(matFile, 'data_noiseless_10_128');
end








% Prints runtime and computation data for solve routine
fprintf('\nRefinement comparison results\n');
fprintf('   n = 128  | primal and dual ref |   only primal ref   |    only dual ref    |     no refinement\n');
fprintf('------------------------------------------------------------------------------------------------------\n');
fprintf(' Iterations | %15i     | %15i     | %15i     | %15i \n', data_noiseless_10_128.nit, data_noiseless_10_128.nit_noDFP, data_noiseless_10_128.nit_noPFD, data_noiseless_10_128.nit_noPFD_noDFP);
fprintf('       FFTs | %15i     | %15i     | %15i     | %15i \n', data_noiseless_10_128.nfft, data_noiseless_10_128.nfft_noDFP, data_noiseless_10_128.nfft_noPFD, data_noiseless_10_128.nfft_noPFD_noDFP);
fprintf('Sig rel err | %15.5e     | %15.5e     | %15.5e     | %15.5e\n', data_noiseless_10_128.rel_err, data_noiseless_10_128.rel_err_noDFP, data_noiseless_10_128.rel_err_noPFD, data_noiseless_10_128.rel_err_noPFD_noDFP);
fprintf('   Runtimes | %15.4f sec | %15.4f sec | %15.4f sec | %15.4f sec \n', data_noiseless_10_128.time, data_noiseless_10_128.time_noDFP, data_noiseless_10_128.time_noPFD, data_noiseless_10_128.time_noPFD_noDFP);



end