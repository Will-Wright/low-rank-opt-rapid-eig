function [x, r, info] = saga_sd(A, b, varargin)
% low-rank-opt main algorithm solver by Michael Friedlander, 2016. 
%
% This solver has been modified in the following ways:
%		 Parameter 				Use
%		 printRelErr		   prints true relative error when user passes true signal 'x0'
%		 logIterData			logs iterate data: 'y', 'x', and eigen-solver tols
%      tMeas, tMain, etc.	logs subroutine runtimes
% Modifications are commented with ~WW (Will Wright, 2017)

import util.*;

% Save and reset global FFT counters
global nfft2 nifft2 ndwt nidwt;
if isempty(nfft2 ), nfft2  = 0; end
if isempty(nifft2), nifft2 = 0; end
if isempty(ndwt  ), ndwt   = 0; end
if isempty(nidwt ), nidwt  = 0; end
[dftPrevious,idftPrevious] = deal(nfft2,nifft2);
[dwtPrevious,idwtPrevious] = deal(ndwt ,nidwt );
[nfft2,nifft2,ndwt,nidwt]  = deal(0);
% Track subroutine runtimes. ~WW
global tMeas tMain tObjective tPrimalRecovery tDualRecovery tLinesearch;
[tMain,tObjective,tPrimalRecovery,tDualRecovery,tLinesearch] = deal(0);

n = size(A, 2);

% Parse input parameters.
p = inputParser;
p.addParameter('epsilon', 0);     % noise parameter
p.addParameter('noise_type', 'guassian'); % noise types: 'gaussian', 'dual', 'exponential'
                                          % Function: changes norms
                                          % ('exponential' calls inf/1 norms)
p.addParameter('iterations', 300);
p.addParameter('maxStep', 1e+10);
p.addParameter('minStep', 1e-10);
p.addParameter('feaTol', 1e-5); % feasibility tol
p.addParameter('optTol', 1e-5); % optimality tol
p.addParameter('primalRelErrTol', 1e-5); % primal relative error, for noisy model convergence ~WW
p.addParameter('dualRelErrTol', 1e-4); % dual relative error, for noisy model convergence ~WW
p.addParameter('eigTol', 1e-2); % initial e-val tol for eigs
p.addParameter('eigTolFact', 0.1); % factor to reduct eigTol each itn
p.addParameter('eigIts', 10000); % max number of iterations per eigs
p.addParameter('eigTolMin', eps); % min eigenvalue tol
p.addParameter('eigPlus', []);  % num of Lanczos vecs above cutDim 
p.addParameter('verbosity', 1);
p.addParameter('useAdaptiveEigs', true); % adaptive `eigs` parameter strategy for rapid computation ~WW
p.addParameter('cutDim', 2); % default number of requested eigenvalues
p.addParameter('num_eigs_min', 2);
p.addParameter('num_eigs_max', 30);
p.addParameter('eigs_basis_size', 40);
p.addParameter('num_lin_interp_pts', 4);
p.addParameter('fid', 1);
p.addParameter('scale', []); % RHS of dual constraint: <b,y> = scale
p.addParameter('LSetamin', 0);  % Linesearch parameters
p.addParameter('LSetamax', 1);
p.addParameter('LSdelta', 1e-4);
p.addParameter('LSsigma', .9);
p.addParameter('LSrho', (sqrt(5)+1)/2);
% p.addParameter('LSeta', .85);
p.addParameter('LSeta', 1);
p.addParameter('LSmaxit', 16);   % 0 == no backtracks
p.addParameter('numPrevItsForBestSignal', 20);
p.addParameter('maxTime', 1800);% max number of seconds allowed for run
p.addParameter('StopOnFeasible', false);
p.addParameter('StopOnRelErr', true);
p.addParameter('skipPFD', false); % skips primal recovery until converged based on eigenvalues ~WW
p.addParameter('skipDFP', false);
p.addParameter('skipLS', false); % skips linesearch. changed to 'true' if skipping BB step ~WW
p.addParameter('pfdOpts', struct()); % options for the PfD solver
p.addParameter('dfpOpts', struct()); % options for the DFP solver
% The following parameters are for testing purposes. ~WW
p.addParameter('logIterData', false);
p.addParameter('returnEigsIterData', false);
p.addParameter('printRelErr', false);
p.addParameter('relErrNormType', 'Frobenius'); % options are '1' and 'Frobenius'
p.addParameter('trueSignalx0', []);
p.addParameter('y0', []);
p.addParameter('differentiableObjTol', 1e-5); % Turns off BB step if |lam1 - lam2|/lam1 < tol
p.addParameter('skipBB', false);
p.addParameter('skipBBMaxStepScale', 200);
p.addParameter('skipBBMinStepScale', 1);
p.parse(varargin{:});

epsilon = p.Results.epsilon;
% Sets norms for primal and gauge dual constraints
if strcmp(p.Results.noise_type, 'exponential')
   pr_norm = 1; % du_norm = Inf;
else
   pr_norm = 2; % du_norm = 2;
end
fid = p.Results.fid;
verbosity = p.Results.verbosity;
maxItns = p.Results.iterations;
maxStep = p.Results.maxStep;
minStep = p.Results.minStep;
feaTol = p.Results.feaTol;
optTol = p.Results.optTol;
primalRelErrTol = p.Results.primalRelErrTol;
dualRelErrTol = p.Results.dualRelErrTol;
eigTol = p.Results.eigTol;
eigTolMin = p.Results.eigTolMin;
eigTolFact = p.Results.eigTolFact;
eigIts = p.Results.eigIts;
numPrevItsForBestSignal = p.Results.numPrevItsForBestSignal;
maxTime = p.Results.maxTime;
scale = p.Results.scale;
bNorm = normv(b, pr_norm);
if isempty(scale)
   scale = bNorm - epsilon;
end
StopOnFeasible = p.Results.StopOnFeasible;
StopOnRelErr = p.Results.StopOnRelErr;
skipPFD = p.Results.skipPFD;
skipDFP = p.Results.skipDFP;
skipLS = p.Results.skipLS;
if epsilon > 0 && (skipDFP == false)
   skipDFP = true; % Turns off dual refinement if noise in model. ~WW
   warning('Dual refinement step has been turned off due to noise in model.');
end
pfdOpts = p.Results.pfdOpts;
dfpOpts = p.Results.dfpOpts;
logIterData = p.Results.logIterData;
iterData = struct('d', cell(1), 'dLinesearch', cell(1), 'dRefined', cell(1), ...
                  'eigsTol', cell(1)', 'eigsTries', cell(1), 'eigsTolLinesearch', cell(1), ...
                  'eigsTriesLinesearch', cell(1), 'eigsTolRefined', cell(1), ...
                  'eigsTriesRefined', cell(1), 'infoEigsFFTs', cell(1), ...
                  'stepInit', cell(1), 'stepLS', cell(1), ...
                  'skipBB', cell(1), 'lam_diff', cell(1));
if logIterData
   iterData.x = cell(1);
   iterData.y = cell(1);
   iterData.yLinesearch = cell(1);
   iterData.yRefined = cell(1);
%							 'iterDataDescription', 'Cell 1: y0 = proj(b), no linesearch. order: LS, y, refined. last LS = y'
% Saves previous signal iterates to select best at termination
else
   xPrevIts = zeros(A.n1, A.n2, numPrevItsForBestSignal);
end
returnEigsIterData = p.Results.returnEigsIterData;
printRelErr = p.Results.printRelErr;
relErrNormType = p.Results.relErrNormType;
if printRelErr && ~isempty(p.Results.trueSignalx0)
	x0 = p.Results.trueSignalx0;
else
	printRelErr = false;
   x0 = [];
end
y0 = p.Results.y0;
differentiableObjTol = p.Results.differentiableObjTol;
skipBB = p.Results.skipBB;
if skipBB
   skipLS = true; % if skipping BB step, linesearch is not a valid option.
end
skipBBMaxStepScale = p.Results.skipBBMaxStepScale;
skipBBMinStepScale = p.Results.skipBBMinStepScale;

% Quiet subsolvers if saga_sd is quiet.
if verbosity >= 2
   pfdOpts.verbose = true;
   dfpOpts.verbose = true;
end

% Linesearch parameters.
LS.etamin = p.Results.LSetamin;
LS.etamax = p.Results.LSetamax;
LS.delta = p.Results.LSdelta;
LS.sigma = p.Results.LSsigma;
LS.rho = p.Results.LSrho;
LS.mu = p.Results.maxStep;
LS.eta = p.Results.LSeta;
LS.maxit = p.Results.LSmaxit;
LS.nls = 0;

% Exit conditions (constants).
stat = 0;
EXIT_OPTIMAL = 1;
EXIT_ITERATIONS = 2;
EXIT_INFEASIBLE_ITERATE = 3;
EXIT_TIME_OUT = 4;
EXIT_PRIMAL_INFEASIBLE = 5;
EXIT_OBJECTIVE_ERROR = 6;
EXIT_PRIMAL_FEASIBLE = 7;
EXIT_RELATIVE_ERRORS_CONVERGED = 8;
EXIT_MSG = {
   'Unknown termination condition!'
   'Optimal solution found'
   'Too many iterations'
   'Infeasible dual iterate'
   'Time out'
   'Primal infeasible'
   'Objective error'
   'Feasible solution found (and requested)'
   'Primal and dual relative errors converged'
};

% Initialize counters.
nit = 0; % no. of iterations
nfe = 0; % no. of function evals

% Initialize eigenvalue options.
if isempty(p.Results.eigPlus)
   eigPlus = double(~A.isrealin);
else
   eigPlus = p.Results.eigPlus;
end
eigsOpts = struct('issym', true, ...
                  'isreal', A.isrealin, ...
                  'maxit', eigIts, ...
                  'k', p.Results.cutDim, ...
                  'kPrev', p.Results.cutDim, ...
                  'returnEigsIterData', returnEigsIterData, ...
                  'useAdaptiveEigs', p.Results.useAdaptiveEigs);
if ~isempty(p.Results.eigs_basis_size)
   eigsOpts.p = p.Results.eigs_basis_size;
else
   eigsOpts.p = min(n, max(2*p.Results.cutDim+eigPlus, 20));
end
% Sets parameters for adaptive eigs method
if eigsOpts.useAdaptiveEigs
   eigsOpts.adaptiveEigsOpts.num_eigs_min = p.Results.num_eigs_min;
   eigsOpts.adaptiveEigsOpts.num_eigs_max = p.Results.num_eigs_max;
   eigsOpts.adaptiveEigsOpts.num_lin_interp_pts = p.Results.num_lin_interp_pts;
   eigsOpts.adaptiveEigsData.num_eigs = p.Results.num_eigs_min;
   eigsOpts.adaptiveEigsData.num_matvecs = [];
   eigsOpts.adaptiveEigsData.EMEP_iter = 1;
   eigsOpts.k = p.Results.num_eigs_min;
   eigsOpts.kPrev = p.Results.num_eigs_min;
end
time.objective = 0;
time.solver = 0;

% Get on your bikes and ride!
tStart = tic;

%----------------------------------------------------------------------
% Initialize primal dual iterates, and associated gradient and step.
% 1. Construct a feasible point y, ie, <b, y> = scale.
%----------------------------------------------------------------------
tMeas = tic;
if isempty(y0)
   y = projection(b);
else
   y = y0;
end
tMain = tMain + toc(tMeas); tMeas = tic;
[duObj, g, v, lamDiff, d, eTries, eigsOpts, infoEigsFFTs] = objective(y, A, [], eigsOpts, eigTol); nfe=nfe+1;
duObjPrev = nan;
gNorm = normv(g, inf);
tObjective = tObjective + toc(tMeas); tMeas = tic;

iterData.eigsTol{nit+1} = eigTol;
iterData.eigsTries{nit+1} = eTries;
iterData.eigsTolLinesearch{nit+1} = [];
iterData.eigsTriesLinesearch{nit+1} = [];
iterData.d{nit+1} = d;
iterData.infoEigsFFTs_y{nit+1} = infoEigsFFTs;
iterData.skipBB{nit+1} = skipBB;
if logIterData
	iterData.y{nit+1} = y;
   iterData.yLinesearch{nit+1} = [];
end

% Initialize linesearch quantities.
if gNorm < (1 / maxStep)
   step = maxStep;
else
   step = min( maxStep, max(minStep, 1/gNorm) );
end
LS.C = duObj;
LS.Q = 1;

%----------------------------------------------------------------------
% Log header.
%----------------------------------------------------------------------
if printRelErr
   logB = '%4i  %12.6e  %12.6e  %8.1e  %8.1e  %8.1e  %7.1e  %10.4e  %10.4e  %8.1e  %7.1e  %8.1e  %8.1e  %3s%8.1e    %4i  %3i    %6i  %6i\n';
   logH = '%4s  %12s  %12s  %8s %8s %8s  %7s  %10s  %10s   %7s  %7s  %8s  %8s %3s %8s %6s  %3s    %6s  %6s\n';
else
   logB = '%4i  %12.6e  %12.6e  %8.1e  %8.1e  %8.1e  %7.1e  %8.1e  %7.1e  %8.1e  %8.1e %3s %8.1e    %4i  %3i  %6i  %6i\n';
   logH = '%4s  %12s  %12s  %8s %8s %8s  %7s  %7s  %8s  %8s  %8s %3s %8s %6s  %3s  %6s  %6s\n';
end
if verbosity
   fprintf(fid,'\n');
   fprintf(fid,' %-20s: %6ix%6i %5s'   ,'matrix size'      ,n,n  ,'');
   fprintf(fid,' %-20s: %13.2e %5s\n'  ,'feasibility tol'  ,feaTol,'');
   fprintf(fid,' %-20s: %13i %5s'      ,'no. observations' ,numel(b),'');
   fprintf(fid,' %-20s: %13.2e %5s\n'  ,'optimality tol'   ,optTol   ,'');
   fprintf(fid,' %-20s: %13.2e %5s'    ,'max step'         ,maxStep,'');
   fprintf(fid,' %-20s: %13.2e %5s\n'  ,'scaling'          ,scale,'');
   fprintf(fid,' %-20s: %13.2e %5s'    ,'min step'         ,minStep,'');
   fprintf(fid,' %-20s: %13.2e %5s\n'  ,'norm b'           ,bNorm,'');
   if eigsOpts.useAdaptiveEigs
      fprintf(fid,' %-20s: %13s %5s'      ,'no. of eigs'      ,'adaptive','');
   else
      fprintf(fid,' %-20s: %13i %5s'      ,'no. of eigs'      ,eigsOpts.kPrev,'');
   end
   fprintf(fid,' %-20s: %13.2e %5s\n'  ,'epsilon'          ,epsilon,'');
   fprintf(fid,' %-20s: %13i %5s'      ,'max its per eigs' ,eigIts,'');
   fprintf(fid,' %-20s: %13.2e %5s\n'  ,'eigTol decrease fact',eigTolFact,'');
   fprintf(fid,'\n');
end

%----------------------------------------------------------------------
% MAIN LOOP.
%----------------------------------------------------------------------

tMain = tMain + toc(tMeas); tMeas = tic;
% Primal recovery.
if skipPFD
   brhs = b-epsilon*y/normv(y, pr_norm);
   nx = sqrt(scale)*sqrt(max(0, rdot(g,brhs))) / normv(g); % scaling based on WFLOW method
   % nx = sqrt(scale/duObj); % scaling based on primal-gauge dual optimality
   x = nx*v;
   r = b - A.forward(x);
   pfdstats.pgNorms = nan;
   pfdstats.iterations = 0;
else
   [x, r, pfdstats] = primal_recovery(A, b, epsilon, y, v, g);
   tPrimalRecovery = tPrimalRecovery + toc(tMeas); tMeas = tic;
end
if logIterData
   iterData.x{nit+1} = x;
else
   if isequal(numel(x),A.n)
      xPrevIts(:, :, nit+1) = reshape(x, [A.n1, A.n2]);
   else
      xPrevIts(:, :, nit+1) = nan(A.n1, A.n2);
   end
end


% Improve dual estimate.
spaced = '';
if skipDFP
   dfpstats.fOut = 0;
   dfpstats.iterations = 0;
else
   tMain = tMain + toc(tMeas); tMeas = tic;
   [yT, duObjT, gT, vT, lamDiffT, dfpstats, dT, eTol, eTries, infoEigsFFTs] ...
         = dual_recovery(A, b, x, y, v);
	tDualRecovery = tDualRecovery + toc(tMeas); tMeas = tic;
   if duObjT < LS.C
      yPrev     = y;
      LS.C      = (duObjT-duObj)/LS.Q+LS.C;
      y         = yT;
      d         = dT;
      duObjPrev = duObj;
      duObj     = duObjT;
      g         = gT;
      v         = vT;
      lamDiff   = lamDiffT;
      spaced    = 'Y';
      iterData.dRefined{nit+1} = d;
		iterData.eigsTolRefined{nit+1} = eTol;
      iterData.eigsTriesRefined{nit+1} = eTries;
      iterData.infoEigsFFTs_yRefined{nit+1} = infoEigsFFTs;
      if logIterData
	      iterData.yRefined{nit+1} = y;
      end
   end
end

while true

   % Computes approximate dual objective if eigenvalue solver didn't converge.
   if isnan(duObj)
      duObj = real(vec(A.adjoint(y, v(:)))'*v(:)) / real(v(:)'*v(:));
      % Avoids termination if dual objective value is invalid.
      if duObj < 0
         duObj = NaN;
      end
   end

   % Quantities needed for log, test exit conditions, and determine if BB step.
   primalRes    = normv(r);
   duFeas       = log(min(1, rdot(y,b) - epsilon*normv(y)));
   prFeas       = max(0,primalRes-epsilon);
   prObj        = normv(x)^2;
   mGap         = abs(1-prObj*duObj/scale);
   eigVal = duObj/scale; eigValPrev = duObjPrev/scale;
   duObjRelErr  = abs(eigVal - eigValPrev)/abs(eigVal);
   testPrFeas   = prFeas <= feaTol*(1+normv(b));
   testOptiml   = abs(mGap) <= optTol;
   if nit == 0
      prRelErr = NaN;
      duRelErr = NaN;
      testPrimalRelErr = false;
      testDualRelErr = false;
   else
      prRelErr = abs(primalRes - primalResPrev) / abs(primalRes);
      duRelErr = norm(y(:) - yPrev(:)) / norm(y(:));
      testPrimalRelErr = prRelErr <= primalRelErrTol;
      testDualRelErr = duRelErr <= dualRelErrTol;
   end

   iterData.res_duRelErr(nit+1) = duRelErr;
   iterData.res_prRelErr(nit+1) = prRelErr;
   iterData.res_primal(nit+1) = primalRes;
   iterData.res_duFeas(nit+1) = duFeas;
   iterData.res_prFeas(nit+1) = prFeas;
   iterData.res_mGap(nit+1) = mGap;
   iterData.res_duObjRelErr(nit+1) = duObjRelErr;


   % Determines if BB step is not applicable (i.e., obj is not differentiable)
   % or linesearch is giving invalid results.
   if ~skipBB && (isnan(duObj) || isnan(mGap) || lamDiff < differentiableObjTol)
      skipBB = true;
      skipLS = true;
   end

   if ~stat  &&  duObj < 0
      stat = EXIT_PRIMAL_INFEASIBLE;
   end
   
   if ~stat  &&  duFeas > feaTol
      stat = EXIT_INFEASIBLE_ITERATE;
   end
  
   if ~stat  &&  testPrFeas && testOptiml
      stat = EXIT_OPTIMAL;
   end

   if ~stat  &&  testPrFeas  &&  StopOnFeasible
      stat = EXIT_PRIMAL_FEASIBLE;
   end
   
   if ~stat  &&  StopOnRelErr  &&  testPrimalRelErr  &&  testDualRelErr ...
      && testPrFeas
      stat = EXIT_RELATIVE_ERRORS_CONVERGED;
   end

   if ~stat  &&  StopOnRelErr  && testDualRelErr  &&  skipPFD   
      stat = EXIT_RELATIVE_ERRORS_CONVERGED;
   end

   if ~stat  &&  nit >= maxItns
      stat = EXIT_ITERATIONS;
   end

   if ~stat  &&  toc(tStart) >= maxTime
      stat = EXIT_TIME_OUT;
   end

   %------------------------------------------------------------------
   % Print log and act on exit conditions.
   %------------------------------------------------------------------
   if ~isempty(x0)
      iterData.x_error_rel(nit+1) = util.hermitianerror(x(:), x0(:), 'fro') / norm(x0(:))^2;
   end
   if verbosity
      if mod(nit, 25) == 0
         fprintf(fid,'\n');
         if printRelErr
            fprintf(fid,logH,...
               'itn','1/pr Obj','du Obj','mGap','prRelDiff','duRelDiff','l1-l2','xErr','xRelErr','step',...
               'prFeas','PFDpgNrm','dfp','spc','eigTol','numEigs','nls',...
               'PFDits','DFPits');
				 else
            fprintf(fid,logH,...
               'itn','1/pr Obj','du Obj','mGap','prRelDiff','duRelDiff','l1-l2','step',...
               'prFeas','PFDpgNrm','dfp','spc','eigTol','numEigs','nls',...
               'PFDits','DFPits');
				 end
      end
      if printRelErr
         if numel(x) == numel(x0)
            if strcmp(relErrNormType, 'Frobenius')
               x_error = util.hermitianerror(x(:), x0(:), 'fro');
               x_error_rel = x_error / norm(x0(:))^2;
            else
               x_error = min(norm(x(:) - x0(:), 1), norm(x(:) + x0(:), 1));
               x_error_rel = x_error / norm(x0(:), 1);
            end
         else
            x_error = nan;
            x_error_rel = nan;
         end
         % Note: norm(v, 2)^2 = norm(v*v', 'fro')        
         fprintf(fid, logB, ...
            nit, 1/prObj, duObj/scale, mGap, prRelErr, duRelErr, lamDiff, ...
            x_error, x_error_rel, step,   prFeas,...
            pfdstats.pgNorms(end), realsqrt(2*dfpstats.fOut), ...
            spaced, eigTol, eigsOpts.kPrev, LS.nls, pfdstats.iterations,...
            dfpstats.iterations);
      else
         fprintf(fid, logB, ...
            nit, 1/prObj, duObj/scale, mGap, prRelErr, duRelErr, lamDiff, step, prFeas,...
            pfdstats.pgNorms(end), realsqrt(2*dfpstats.fOut), ...
            spaced, eigTol, eigsOpts.kPrev, LS.nls, pfdstats.iterations,...
            dfpstats.iterations);
      end
   end
   
   if stat
      if skipPFD
         skipPFD = false; % Turns on primal recovery for final iterate to return solution signal ~WW
      else
         break
      end
   end
   
   %==================================================================
   % New iteration starts here.
   %==================================================================
   nit = nit + 1;
   primalResPrev = primalRes;
   duObjPrev = duObj;
   yPrev = y;

   %------------------------------------------------------------------
   % Linesearch.
   %------------------------------------------------------------------
   if skipLS
      % uses 'linesearch' to compute single-step projected gradient
      LS.maxit = 0;
   end
   eigTol = max(eigTolMin, min(eigTol, eigTolFact*min([prFeas, mGap, eigTol])));
   yk = y; gk = g;
   try
      tMain = tMain + toc(tMeas); tMeas = tic;
      [duObj, y, g, v, lamDiff, d, LS, stepLS, eigsOpts, iterData] = linesearch(@objective, @projection,...
         yk, gk, step, A, v, eigsOpts, eigTol, LS, iterData, logIterData, nit);
      tLinesearch = tLinesearch + toc(tMeas); tMeas = tic;
   catch
      stat = EXIT_OBJECTIVE_ERROR;
   end
   nfe = nfe + LS.nls;

   % Primal recovery.
   if skipPFD
      brhs = b-epsilon*y/normv(y, pr_norm);
      nx = sqrt(scale)*sqrt(max(0, rdot(g,brhs))) / normv(g); % scaling based on WFLOW method
      %nx = sqrt(scale/duObj); % scaling based on primal-gauge dual optimality
      x = nx*v;
      r = b - A.forward(x);
      pfdstats.pgNorms = nan;
      pfdstats.iterations = 0;
   else
      tMain = tMain + toc(tMeas); tMeas = tic;
      [x, r, pfdstats] = primal_recovery(A, b, epsilon, y, v, g);
      tPrimalRecovery = tPrimalRecovery + toc(tMeas); tMeas = tic;
   end
	iterData.eigsTol{nit+1} = eigTol;
   iterData.eigsTries{nit+1} = eTries;
   iterData.d{nit+1} = d;
   iterData.stepInit{nit+1} = step;
   iterData.stepLS{nit+1} = stepLS;
   iterData.skipBB{nit+1} = skipBB;
   if logIterData
      iterData.x{nit+1} = x;
      iterData.y{nit+1} = y;
   % Selects best signal among recent iterates, overwrites oldest signal
   else
      nitMod = mod(nit+1, numPrevItsForBestSignal);
      if nitMod == 0
         nitMod = numPrevItsForBestSignal;
      end
      if isequal(numel(x), A.n)
         xPrevIts(:,:,nitMod) = reshape(x, [A.n1, A.n2]);
      else
         xPrevIts(:,:,nitMod) = nan(A.n1, A.n2);
      end
   end
   
   % Spacer step.
   spaced = '';
   if ~skipDFP
      tMain = tMain + toc(tMeas); tMeas = tic;
      [yT, duObjT, gT, vT, lamDiffT, dfpstats, dT, eTol, eTries, infoEigsFFTs] ...
         = dual_recovery(A, b, x, y, v);
      tDualRecovery = tDualRecovery + toc(tMeas); tMeas = tic;
      if duObjT < LS.C
         LS.C    = (duObjT-duObj)/LS.Q+LS.C;
         y       = yT;
         d       = dT;
         duObj   = duObjT;
         g       = gT;
         v       = vT;
         lamDiff = lamDiffT;
         spaced  = 'Y';
         iterData.dRefined{nit+1} = d;
         iterData.eigsTolRefined{nit+1} = eTol;
         iterData.eigsTriesRefined{nit+1} = eTries;
         iterData.infoEigsFFTs_yRefined{nit+1} = infoEigsFFTs;
         if logIterData
            iterData.yRefined{nit+1} = y;
         end
      end
   end

   if skipBB
      maxStep = skipBBMaxStepScale / (nit + 1);
      minStep = skipBBMinStepScale / (nit + 1);
   end
   % Update Barzilai-Borwein step length.
   dy = y - yk;
   dg = g - gk;
   step = min(maxStep, max(minStep, rdot(dy,dy)/rdot(dy,dg)));

end % while true
time.solver = toc(tStart);
tMain = tMain + toc(tMeas);

% Selects optimal x based on minimal duality gap value
if epsilon == 0 || skipPFD
   % Returns last signal iterate if problem is noiseless or primal
   % refinement is off
elseif nit+1 <= numPrevItsForBestSignal || logIterData
   [~, xIdx] = min(iterData.res_mGap);
   if logIterData
      x = iterData.x{xIdx};
   else
      x = xPrevIts(:, :, xIdx);
   end
else
   nitMod = mod(nit+1, numPrevItsForBestSignal);
   mGapVals = iterData.res_mGap(nit+2-numPrevItsForBestSignal : nit+1);
   [~, minIdx] = min(mGapVals);
   xIdx = mod(nitMod + minIdx, numPrevItsForBestSignal);
   if xIdx == 0
      xIdx = numPrevItsForBestSignal;
   end
   x = xPrevIts(:, :, xIdx);
end

%----------------------------------------------------------------------
% Log footer.
%----------------------------------------------------------------------
nFFT = nfft2+nifft2;
nDWT = ndwt +nidwt ;

if verbosity
fprintf(fid,'\n EXIT -- %s\n\n',EXIT_MSG{stat+1});
fprintf(fid,' %-22s: %10i %5s'  ,'no. of measurements', A.nForward,'');
fprintf(fid,' %-22s: %10.2f %5s\n','time solve (sec)', time.solver,'');
fprintf(fid,' %-22s: %10i %5s'  ,'no. of adjoints', A.nAdjoint,'');
fprintf(fid,' %-22s: %10.2f %5s\n','time objective (sec)', time.objective,'');
fprintf(fid,' %-22s: %10i %5s'  ,'no. of FFTs', nFFT,'');
fprintf(fid,' %-22s: %10i %5s\n','no. of PFDObj. calls', A.nPFDObjective,'');
fprintf(fid,' %-22s: %10i %5s'  ,'no. of DWTs', nDWT,'');
fprintf(fid,' %-22s: %10i %5s\n','no. of DFPObj. calls', A.nDFPObjective,'');
fprintf(fid,'\n');
end

% Gather exit info.
info.prFeas = prFeas;
info.duObj = duObj;
info.prObj = prObj;
info.mGap = mGap;
info.duObjRelErr = duObjRelErr;
info.epsilon = epsilon;
info.scale = scale;
info.nit = nit;
info.nfe = nfe;
info.nfft = nFFT;
info.ndwt = nDWT;
info.nAdjoint = A.nAdjoint;
info.nMeasure = A.nForward;
info.main.tMain = tMain;
info.main.tObjective = tObjective;
info.main.tPrimalRecovery = tPrimalRecovery;
info.main.tDualRecovery = tDualRecovery;
info.main.tLinesearch = tLinesearch;
info.hop.nForward = A.nForward;
info.hop.tForward = A.tForward;
info.hop.nAdjoint = A.nAdjoint;
info.hop.tAdjoint = A.tAdjoint;
info.hop.nMeasure = A.nForward;
info.hop.nGDObjective = A.nGDObjective;
info.hop.tGDObjective = A.tGDObjective;
info.hop.nEigs = A.nEigs;
info.hop.tEigs = A.tEigs;
info.hop.nMultsInEigs = A.nMultsInEigs;
info.hop.nFFTsInEigs = A.nFFTsInEigs;
info.eigsOpts = eigsOpts;
info.hop.nPFDObjective = A.nPFDObjective;
info.hop.tPFDObjective = A.tPFDObjective;
info.hop.nDFPObjective = A.nDFPObjective;
info.hop.tDFPObjective = A.tDFPObjective;
info.time = time.solver;
info.stat = stat;
info.status = EXIT_MSG{stat+1};
info.y = y;
info.iterData = iterData;

% restore global FFT and DWT counters
[nfft2,nifft2,ndwt,nidwt] = deal(dftPrevious,idftPrevious,dwtPrevious,idwtPrevious);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NESTED FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   function [fd, g, v, lamDiff, d, eTries, eigsOpts, infoEigsFFTs] = objective(y, A, v, eigsOpts, eTol)
      %OBJECTIVE  The dual objective value lambda_1(A'*y).
      %
      % JUL 6 2015: EIGS routine can fail if the parameter `p` is too
      % small. The catch statement below will detect a failure and
      % increase `p` 3 times in a row. Otherwise, will throw the error.

      tObjStart    = tic;
      eigsOptsTemp = eigsOpts;
      eigsOptsTemp.tol = eTol;
      eTries = 0; % No. of times eigenvalue solver called.
      while true
         try
            [fd,g,V,d,infoEigsFFTs] = A.gdobjective(y,eigsOptsTemp.k,eigsOptsTemp,v);
            break
         catch ME
            if eTries <= 3 && strcmp(ME.identifier,'MATLAB:eigs:ARPACKroutineErrorMinus14')
               eigsOptsTemp.p = 2*eigsOptsTemp.p;
               eTries = eTries + 1;
            else
               causeException = MException('MATLAB:saga:objective','eigs failed');
               ME = addCause(ME,causeException);
               rethrow(ME)
            end
         end
      end
      v = V(:,1);
      if numel(d) > 1
         lamDiff = (d(1)-d(2))/d(1);
      else
         lamDiff = NaN;
      end

      % updates parameters for adaptive eigs method
      if eigsOpts.useAdaptiveEigs
               
%   eigsOpts.adaptiveEigsOpts.num_eigs_min = p.Results.num_eigs_min;
%   eigsOpts.adaptiveEigsOpts.num_eigs_max = p.Results.num_eigs_max;
%   eigsOpts.adaptiveEigsOpts.num_lin_interp_pts = p.Results.num_lin_interp_pts;
%   eigsOpts.adaptiveEigsData.num_eigs = p.Results.num_eigs_min;
%   eigsOpts.adaptiveEigsData.num_matvecs = [];
%   eigsOpts.adaptiveEigsData.EMEP_iter = 1;
%   eigsOpts.k = p.Results.num_eigs_min;
%   eigsOpts.kPrev = p.Results.num_eigs_min;   
   
         EMEP_iter = eigsOpts.adaptiveEigsData.EMEP_iter;
         eigsOpts.adaptiveEigsData.num_matvecs(EMEP_iter) = infoEigsFFTs.nMatVec;         
         
         [num_eigs_update] = adaptive_eigs_params(...
            eigsOpts.adaptiveEigsData.num_matvecs, ...
            eigsOpts.adaptiveEigsData.num_eigs, ...
            EMEP_iter, eigsOpts.adaptiveEigsOpts);
         
         eigsOpts.kPrev = eigsOpts.k;
         eigsOpts.k = num_eigs_update;
         eigsOpts.adaptiveEigsData.num_eigs(EMEP_iter + 1) = num_eigs_update;
         eigsOpts.adaptiveEigsData.EMEP_iter = eigsOpts.adaptiveEigsData.EMEP_iter + 1;
      end
      time.objective = time.objective + toc(tObjStart);
   end

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   function z = projection(z)
      % Projection operator onto the constraint
      % <b,y> - epsilon||y|| >= scale.
      z = project(z, b, epsilon, scale);
%     z = z + ((scale-rdot(b,z))/bNorm^2)*b; % projection onto <b,y>=1.
   end

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   function [x, r, pfdstats] = primal_recovery(A, b, epsilon, y, v, g)
      %PRIMAL_REOCOVERY  Recovery primal estimate from gradient.
      [x, pfdstats] = pfd(A, b, epsilon, y, v, g, pfdOpts);
%      [x, pfdstats] = pfd(A, b, 0, y, v, g, pfdOpts);
      r = b - A.forward(x);
   end

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   function [y, duObj, g, v, lamDiff, dfpstats, d, eTol, eTries, infoEigsFFTs] = dual_recovery(A, b, x, y, v)
      %DUAL_RECOVERY  Recover an dual estimate from primal.
      normr = util.normv(r);
      if epsilon > 0 && normr > 0
         y0 = r*(scale/(util.rdot(r,b)-epsilon*normr));
      else
         y0 = y;
      end
      [y, dfpstats] = dfp(A, b, epsilon, scale, x, y0, dfpOpts);
		eTol = eps;
      [duObj, g, v, lamDiff, d, eTries, eigsOpts, infoEigsFFTs] = objective(y, A, v, eigsOpts, eTol);
   end



end % function saga_sd

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LOCAL FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [duObj, y, g, v, lamDiff, d, LS, step, eigsOpts, iterData] = linesearch(objective, projection, yk, gk, step, A, v, eigsOpts, eigTol, LS, iterData, logIterData, nit)
   % Full nonmonote linesearch.
   import util.*;
   y = projection( yk - step * gk );
   [duObj, g, v, lamDiff, d, eTries, eigsOpts, infoEigsFFTs] = objective(y, A, v, eigsOpts, eigTol);
   iterData.dLinesearch{nit+1}{1} = d;
   iterData.eigsTolLinesearch{nit+1}{1} = eigTol;
   iterData.eigsTriesLinesearch{nit+1}{1} = eTries;
   iterData.infoEigsFFTs_yLS{nit+1}{1} = infoEigsFFTs;
   iterData.lam_diff{nit+1}{1} = lamDiff;
   if logIterData
      iterData.yLinesearch{nit+1}{1} = y;
   end
   dy = y - yk;
   done = false;
   if step <= LS.mu && duObj <= LS.C + LS.delta * rdot(gk,dy)
      % Armijo satisfied
      if rdot(g,dy) >= LS.sigma * rdot(gk,dy)
         % Wolfe satisfied
         done = true;
      else
         increaseStep = true;
      end
   else
      increaseStep = false;
   end
   nls = 1;
   while ~done && nls <= LS.maxit
      if increaseStep
         step = step * LS.rho;
         yT = projection( yk - step * gk );
         [duObjT, gT, vT, lamDiffT, dT, eTriesT, eigsOpts, infoEigsFFTsT] = objective(yT, A, v, eigsOpts, eigTol);
         iterData.dLinesearch{nit+1}{nls+1} = dT;
         iterData.eigsTolLinesearch{nit+1}{nls+1} = eigTol;
         iterData.eigsTriesLinesearch{nit+1}{nls+1} = eTriesT;
         iterData.infoEigsFFTs_yLS{nit+1}{nls+1} = infoEigsFFTsT;
         iterData.lam_diff{nit+1}{nls+1} = lamDiffT;
		   if logIterData
       	   iterData.yLinesearch{nit+1}{nls+1} = yT;
    		end
         dyT = yT - yk;
         armijoSatisfied = (step <= LS.mu) && (duObjT <= LS.C + LS.delta * rdot(gk,dyT));
         done = ~armijoSatisfied || (rdot(gT,dyT) >= LS.sigma * rdot(gk,dyT));
         if armijoSatisfied || nls == LS.maxit
            y = yT; d = dT; duObj = duObjT; g = gT; v = vT; lamDiff = lamDiffT;
         end
      else
         step = step / LS.rho;
         y = projection( yk - step * gk );
         [duObj, g, v, lamDiff, d, eTries, eigsOpts, infoEigsFFTs] = objective(y, A, v, eigsOpts, eigTol);
         iterData.dLinesearch{nit+1}{nls+1} = d;
         iterData.eigsTolLinesearch{nit+1}{nls+1} = eigTol;
         iterData.eigsTriesLinesearch{nit+1}{nls+1} = eTries;
         iterData.infoEigsFFTs_yLS{nit+1}{nls+1} = infoEigsFFTs;
         iterData.lam_diff{nit+1}{nls+1} = lamDiff;
         if logIterData
            iterData.yLinesearch{nit+1}{nls+1} = y;
         end
         dy = y - yk;
         done = (duObj <= LS.C + LS.delta * rdot(gk,dy));
      end
      nls = nls + 1;
   end
   % Update linesearch quantities.
   LS.nls = nls;
   LS.Q = LS.eta * LS.Q + 1;
   LS.C = ( duObj + LS.C * (LS.Q - 1) ) / LS.Q;
end

