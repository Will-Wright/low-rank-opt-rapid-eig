function [] = noisyrandom_noisypl_vs_evolvingmatricespl()
% experiments.test.noisyrandom_noisypl_vs_evolvingmatricespl
%
% Runs tests on noisy random phase retrieval problems using
% two scripts: original script 'noisypl' (from Friedlander, 2016) and
% new script 'evolvingmatricespl'.  Verifies results match for both runs.
%
% Note: results for number of ffts in `noisypl` in this directory differ 
% slightly from the version in original `low-rank-opt` directory 
% due to random number generator behavior causing `eigs` to seed 
% differently for each call to `noisypl`.  This is an inconsequential
% difference and unavoidable.


% Sets generated data parameters
n = 128;
resizeImage = 1;
itersEta = [0.005, 0.05, 0.5];
itersL = [6, 9, 12];
signal = 'dual'; % 'gaussian' DRAMATICALLY INCREASES NUMBER OF ITERATIONS
mask = 'octanary'; % 'gaussian' minimal
normalize_signal = true; % false has minimal change
normalize_mask = true; % false has minimal change
seed = 1;

% Sets solver parameters
StopOnRelErr = false;
eigIts = 300;
differentiableObjTol = 0;
logIterData = false;
verbosity = 0;
wflowVerbosity = 0;

% Runs tests using saga_sd and wflow based on (Friedlander, 2016) script experiments.pl
[T, Tsumm] = experiments.noisypl('test', true);
fprintf('\n');

% Runs tests using saga_sd
T_row = 1;
PASSED_ALL_TESTS = true;
for i = 1:length(itersEta)
   for j = 1:length(itersL)
      noise_ratio = itersEta(i);
      L = itersL(j);
      [results] = experiments.evolvingmatricespl('noise_ratio', noise_ratio, 'n', n, 'L', L, ...
                                                 'signal', signal, 'mask', mask, ...
                                                 'resizeImage', resizeImage, ...
                                                 'iterations', n^2, ...
                                                 'normalize_signal', normalize_signal, ...
                                                 'normalize_mask', normalize_mask, ...
                                                 'StopOnRelErr', StopOnRelErr, ...
                                                 'differentiableObjTol', differentiableObjTol, ...
                                                 'eigTol', 1e-2, ...
                                                 'eigTolFact', 0.1, ...
                                                 'eigIts', eigIts, ...
                                                 'seed', seed, ...
                                                 'logIterData', logIterData, 'verbosity', verbosity);
      nfft = results.saga_sd.solInfo.nfft;
      if T{T_row, 9} ~= nfft
         fprintf('FAIL TEST: Incorrect nFFT for L = %3i, noise_ratio = %1.4d for saga_sd\n', itersL(j), itersEta(i));
         [T{T_row, 9}, nfft]
         PASSED_ALL_TESTS = false;
      end
      T_row = T_row + 1;
   end
end   
if PASSED_ALL_TESTS
   fprintf('PASS TEST: nFFT identical for saga_sd for all tests in experiments.noisypl and this script.\n');
end

% Runts tests using wflow
T_row = 10;
PASSED_ALL_TESTS = true;
for i = 1:length(itersEta)
   for j = 1:length(itersL)
      noise_ratio = itersEta(i);
      L = itersL(j);
      [results] = experiments.evolvingmatricespl('wflowOn', true, 'saga_sdOn', false, ...
                                                 'wflowVerbosity', wflowVerbosity, ...
                                                 'noise_ratio', noise_ratio, 'n', n, 'L', L, ...
                                                 'signal', signal, 'mask', mask, ...
                                                 'resizeImage', resizeImage, ...
                                                 'normalize_signal', normalize_signal, ...
                                                 'normalize_mask', normalize_mask, ...
                                                 'StopOnRelErr', StopOnRelErr, ...
                                                 'differentiableObjTol', differentiableObjTol, ...
                                                 'seed', seed, ...
                                                 'logIterData', logIterData, 'verbosity', verbosity);
      nfft = results.wflow.solInfo.nfft;
      if T{T_row, 9} ~= nfft
         fprintf('FAIL TEST: Incorrect nFFT for L = %3i, noise_ratio = %1.4d for wflow\n', itersL(j), itersEta(i));
         PASSED_ALL_TESTS = false;
      end
      T_row = T_row + 1;
   end
end
if PASSED_ALL_TESTS
   fprintf('PASS TEST: nFFT identical for wflow for all tests in experiments.noisypl and this script.\n');
end

end