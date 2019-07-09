function [results] = evolvingmatricespl(varargin)
% experiments.evolvingmatricespl
%
% Generates noiseless/noisy observations of natural image, then solves the PhaseLift
% model and returns the sequence of dual iterates yk.  Default setting is noisy model.
% Default setting is to solve with saga_sd and not with wflow.
%
% Note: Evolving matrices Ak = A'(yk) are represented as objects of the class
% +hop/pl.m

addpath(genpath('external'));

% Resets RNG.
rng('default');

ip = inputParser;

% Sets generated data parameters
ip.addParameter('resizeImage', 1/1); 
ip.addParameter('L', 6);
ip.addParameter('n', 128);    % This is the default setting for 2016 random results.  Does not correspond to image size.
ip.addParameter('channel', 'gray'); % Can be 'gray', 'red', 'green', 'blue', or 'RGB' to run all 3 experiments
ip.addParameter('noise_ratio', 0.30);
ip.addParameter('signal', 'image');    % Can be 'image' for specific image, or 'gaussian' for generic random problem
ip.addParameter('image_file', 'data/parrot_16k.jpg');
ip.addParameter('noise_type', 'gaussian');   % Can be 'gaussian', 'dual', or 'exponential'
ip.addParameter('mu', 1); % parameter for exponential distribution of noise term
ip.addParameter('mask', 'gaussian');
ip.addParameter('normalize_signal', false);
ip.addParameter('normalize_mask', false);
ip.addParameter('seed', 0);

% Sets saga_sd solver parameters
ip.addParameter('saga_sdOn', true);
ip.addParameter('verbosity', 1);
ip.addParameter('logIterData', true);
ip.addParameter('returnEigsIterData', true);
ip.addParameter('printRelErr', true);
ip.addParameter('trueSignalx0', []);    % x0 saved below after test problem is generated
ip.addParameter('iterations', 300);
ip.addParameter('maxTime', 3*60*60);
ip.addParameter('numPrevItsForBestSignal', 20);
ip.addParameter('eigTol', 1e-2); % initial e-val tol for eigs
ip.addParameter('eigTolFact', 0.1); % factor to reduct eigTol each itn
ip.addParameter('eigIts', 10000); % Default setting was originally 300
ip.addParameter('LSmaxit', 16); % Default setting: 16 iterations
ip.addParameter('scale', 1);  % Default settings: for noisy models scale = 1, for noiseless scale = norm(b) - noise_ratio, with noise_ratio = norm(b - b*)/norm(b)
ip.addParameter('skipPFD', false); % Modification by WW, original default is false
ip.addParameter('skipDFP', true); % Default settings: true for noisy models, false for noiseless
ip.addParameter('skipLS', false); % Modification by WW, original default is false
ip.addParameter('skipBB', false); % Modification by WW, original default is false
ip.addParameter('skipBBMaxStepScale', 200);
ip.addParameter('skipBBMinStepScale', 1);

ip.addParameter('useAdaptiveEigs', false);
ip.addParameter('num_eigs_min', 2);
ip.addParameter('num_eigs_max', 30);
ip.addParameter('num_lin_interp_pts', 4);
ip.addParameter('eigs_basis_size', 40);
ip.addParameter('cutDim', 2);

ip.addParameter('differentiableObjTol', 1e-5); % Modification by WW, turns off Barzilai-Borwein linesearch if objective function deemed nondifferentiable
ip.addParameter('StopOnFeasible', false); % Default settings; false for noisy models, true or false for noiseless
ip.addParameter('StopOnRelErr', true); % Modification by WW, original default is false
ip.addParameter('optTol', 2e-4); % Default settings: for noisy models optTol = 2e-4, for noiseless optTol = 1e-5
ip.addParameter('feaTol', 2e-4); % Default settings: for noisy models optTol = 2e-4, for noiseless optTol = 1e-5
ip.addParameter('primalRelErrTol', 1e-5); % Added by WW
ip.addParameter('dualRelErrTol', 1e-4); % Added by WW
ip.addParameter('pfdOpts', struct('gtol', 2e-7, 'ftol', 0, ...
											 'maxit', 1e3, 'solver', 'minFunc' )); % options for the PFD solver

% Sets wflow parameters
ip.addParameter('wflowOn', false);
ip.addParameter('wflowIterations', 1000);
ip.addParameter('wflowVerbosity', 1);
ip.addParameter('wflowUseOptimalStep', false);
ip.parse(varargin{:});


% Solves all three color channel problems if 'RGB' selected
if strcmp(ip.Results.channel, 'RGB')
   channel_set = ["red", "green", "blue"];
   results.red.params = ip.Results;
   results.green.params = ip.Results;
   results.blue.params = ip.Results;
elseif strcmp(ip.Results.channel, 'gray')
   channel_set = "gray";
   results.params = ip.Results;
elseif strcmp(ip.Results.channel, 'red')
   channel_set = "red";
   results.red.params = ip.Results;
elseif strcmp(ip.Results.channel, 'green')
   channel_set = "green";
   results.green.params = ip.Results;
elseif strcmp(ip.Results.channel, 'blue')
   channel_set = "blue";
   results.blue.params = ip.Results;
end

for channel = channel_set
   if strcmp(ip.Results.signal, 'image')
      fprintf('Solving phase retrieval image recovery with color channel %s.\n', channel);
   end
   
   % Generates experiment data: sensing operator A, noisy observation b, true signal x0
   genOpts = struct('signal', ip.Results.signal, ...
                     'image_file', ip.Results.image_file, 'noise_type', ip.Results.noise_type, ...
                     'resizeImage', ip.Results.resizeImage, ...
                     'L', ip.Results.L, 'n', ip.Results.n, ...
                     'channel', channel, ...
                     'noise_ratio', ip.Results.noise_ratio, 'mu', ip.Results.mu, ...
                     'mask', ip.Results.mask, 'normalize_signal', ip.Results.normalize_signal, ...
                     'normalize_mask', ip.Results.normalize_mask, 'seed', ip.Results.seed);
   [A, b, x0, genInfo] = experiments.gendatapl(genOpts);

   % Sets solver parameters
   epsilon = genInfo.epsilon;
   if strcmp(channel, 'gray')
      results.gendatapl = struct('A', A, 'b', b, 'x0', x0, 'genInfo', genInfo);
      results.params.epsilon = epsilon;
   elseif strcmp(channel, 'red')
      results.red.gendatapl = struct('A', A, 'b', b, 'x0', x0, 'genInfo', genInfo);
      results.red.params.epsilon = epsilon;
   elseif strcmp(channel, 'green')
      results.green.gendatapl = struct('A', A, 'b', b, 'x0', x0, 'genInfo', genInfo);
      results.green.params.epsilon = epsilon;
   elseif strcmp(channel, 'blue')
      results.blue.gendatapl = struct('A', A, 'b', b, 'x0', x0, 'genInfo', genInfo);
      results.blue.params.epsilon = epsilon;
   end
   
   solverOpts = struct('epsilon', epsilon, ...
                        'noise_type', ip.Results.noise_type, ...
                        'verbosity', ip.Results.verbosity, ...
                        'logIterData', ip.Results.logIterData, ...
                        'printRelErr', ip.Results.printRelErr, ...
                        'returnEigsIterData', ip.Results.returnEigsIterData, ...
                        'trueSignalx0', x0, ...
                        'iterations', ip.Results.iterations, ...
                        'maxTime', ip.Results.maxTime, ...
                        'eigTol', ip.Results.eigTol, ...
                        'eigTolFact', ip.Results.eigTolFact, ...
                        'eigIts', ip.Results.eigIts, ...
                        'scale', ip.Results.scale, ...
                        'skipPFD', ip.Results.skipPFD, ...
                        'skipDFP', ip.Results.skipDFP, ...
                        'skipLS', ip.Results.skipLS, ...
                        'skipBB', ip.Results.skipBB, ...
                        'cutDim', ip.Results.cutDim, ...
                        'useAdaptiveEigs', ip.Results.useAdaptiveEigs, ...
                        'eigs_basis_size', ip.Results.eigs_basis_size, ...
                        'differentiableObjTol', ip.Results.differentiableObjTol, ...
                        'optTol', ip.Results.optTol, ...
                        'feaTol', ip.Results.feaTol, ...
                        'primalRelErrTol', ip.Results.primalRelErrTol, ...
                        'dualRelErrTol', ip.Results.dualRelErrTol, ...
                        'StopOnFeasible', ip.Results.StopOnFeasible, ...
                        'StopOnRelErr', ip.Results.StopOnRelErr, ...
                        'pfdOpts', ip.Results.pfdOpts);
%{
   % Original saga_sd options
   solverOpts = struct('epsilon', epsilon, ...
                        'verbosity', ip.Results.verbosity, ...
                        'iterations', ip.Results.iterations, ...
                        'maxTime', ip.Results.maxTime, ...
                        'eigTol', ip.Results.eigTol, ...
                        'eigTolFact', ip.Results.eigTolFact, ...
                        'eigIts', ip.Results.eigIts, ...
                        'scale', ip.Results.scale, ...
                        'skipDFP', ip.Results.skipDFP, ...
                        'cutDim', ip.Results.cutDim, ...
                        'optTol', ip.Results.optTol, ...
                        'feaTol', ip.Results.feaTol, ...
                        'StopOnFeasible', ip.Results.StopOnFeasible, ...
                        'pfdOpts', ip.Results.pfdOpts);   
%}   
   
   wflowSolverOpts.epsilon = solverOpts.epsilon;
   wflowSolverOpts.iterations = ip.Results.wflowIterations;
   wflowSolverOpts.verbosity = ip.Results.wflowVerbosity;
   wflowUseOptimalStep = ip.Results.wflowUseOptimalStep;
   if wflowSolverOpts.epsilon > 0
      wflowSolverOpts.Y = b*A.n;
   end

   
   % Solves PhaseLift model with saga_sd
   if ip.Results.saga_sdOn
      % Original saga_sd implementation
      %[x, r, solInfo] = saga_sd_orig(A, b, solverOpts);
      [x, r, solInfo] = saga_sd(A, b, solverOpts);      
      if strcmp(channel, 'red')
         results.red.saga_sd = struct('x', x, 'r', r, 'solInfo', solInfo);
      elseif strcmp(channel, 'green')
         results.green.saga_sd = struct('x', x, 'r', r, 'solInfo', solInfo);
      elseif strcmp(channel, 'blue')
         results.blue.saga_sd = struct('x', x, 'r', r, 'solInfo', solInfo);
      else
         results.saga_sd = struct('x', x, 'r', r, 'solInfo', solInfo);
      end
   end

   % Solves PhaseLift model with wflow
   if ip.Results.wflowOn
      % The `wflow` routine does its own measuring. Otherwise, we have
      % to redefine its internal operators, which isn't advisable!
      X0 = reshape(x0, A.n1, A.n2);
      if wflowUseOptimalStep
         [X, solInfo] = wflow_opt(conj(A.masks), X0, wflowSolverOpts, 'gTol', 1e-8);
      else
         [X, solInfo] = wflow(conj(A.masks), X0, wflowSolverOpts);
      end

      % Corrects phase of recovered signal.  The wflow implementation is not
      % coded to project iterates onto reals.
      if strcmp('image', ip.Results.signal)
         X = abs(X);
      end
      r = b - A.forward(X);
      x = X(:);
      solInfo.mGap = nan; % not relevant for wflow
      if strcmp(channel, 'red')
         results.red.wflow = struct('x', x, 'r', r, 'solInfo', solInfo);
      elseif strcmp(channel, 'green')
         results.green.wflow = struct('x', x, 'r', r, 'solInfo', solInfo);
      elseif strcmp(channel, 'blue')
         results.blue.wflow = struct('x', x, 'r', r, 'solInfo', solInfo);
      else
         results.wflow = struct('x', x, 'r', r, 'solInfo', solInfo);
      end
   end
end
   
end % function evolvingmatricespl
