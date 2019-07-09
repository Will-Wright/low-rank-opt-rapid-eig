function [results] = noisyimage_spectrumdist(varargin)
% experiments.figure.noisyimage_spectrumdist
%
% Generates plots of particular iterates for evolving matrices Ak = A'(yk)

ip = inputParser;

% Sets generated data parameters
ip.addParameter('results', []); 
ip.parse(varargin{:});
results = ip.Results.results;

L = 5;
noise_ratio = 0.15;

%iters = [1, 3, 15, 30, 45, 60];
iters = [1, 2, 3, 10, 50, 100];


% Computes default problem `results` struct if `results` not provided
if isempty(results)
   n = 4096;
   resizeImage = 1/1; 
   image_file = 'data/parrot_4k.jpg';
  
   % Loads cached (or solves and caches) results of test problem
   pathName = 'cache/figure_noisyimage_spectrumdist/';
   if exist(pathName) ~= 7
      mkdir(pathName);
   end
   matFile = fullfile(pathName, 'results.mat');
   if exist(matFile, 'file') ~= 2
      [results] ...
         = experiments.evolvingmatricespl('L', L, 'resizeImage', resizeImage, ...
                                          'noise_type', 'gaussian', ...
                                          'image_file', image_file, ...
                                          'noise_ratio', noise_ratio, 'iterations', 300, ...
                                          'verbosity', 1, 'logIterData', true, ...
                                          'printRelErr', true, ...
                                          'skipDFP', true, 'StopOnFeasible', false, ...
                                          'returnEigsIterData', true, ...
                                          'feaTol', 2e-4, 'optTol', 2e-4, ...
                                          'primalRelErrTol', 1e-4, ...
                                          'StopOnRelErr', true, ...
                                          'seed', 1);
      save(matFile, 'results');
   else
      load(matFile, 'results');
   end

   % Removes laters iters if solver ends before itn = 90
   totalIters = results.saga_sd.solInfo.nit;
   iters = iters(1:sum(iters<=totalIters));

   % Computes all eigenvalues of A_itn
   if ~isfield(results, 'spectrum')
      results.spectrum = cell(1);
   end
   
   fprintf('\n');
   for iter = iters
      if isempty(results.spectrum{iter})
         fprintf('Computing spectrum of A(%i)...\n', iter);
         results.spectrum{iter} = experiments.exputil.getspectrum(results, iter);
         save(matFile, 'results');
      end

   end
   
   % Removes unnecessary data to save smaller file
   if isfield(results.saga_sd.solInfo.iterData, 'x') ...
         || isfield(results.saga_sd.solInfo.iterData, 'y') ...
         || isfield(results.saga_sd.solInfo.iterData, 'yLinesearch')
      results.saga_sd.solInfo.iterData = rmfield(results.saga_sd.solInfo.iterData, 'x');
      results.saga_sd.solInfo.iterData = rmfield(results.saga_sd.solInfo.iterData, 'y');
      results.saga_sd.solInfo.iterData = rmfield(results.saga_sd.solInfo.iterData, 'yLinesearch');
      save(matFile, 'results');
   end

else
   results.spectrum = cell(1);
   for iter = iters
      fprintf('Computing spectrum of A(%i)...\n', iter);
      results = experiments.exputil.getspectrum(results, iter);
   end
   n = results.gendatapl.A.n;
end


% Plots eigenvalue distribution for A_itn
fprintf('\nAll matrix spectra have been computed.\n');

for i = 1:length(iters)
   iter = iters(i);
   subplot(2, length(iters)/2, i);
   plot(results.spectrum{iter}, 'b.');
   hold on
   plot(zeros(1, length(results.spectrum{iter})), 'r');
   hold off
%{
   if itn == 1
      max_eig_vals = length(results.spectrum{itn});
      plot(6:max_eig_vals, results.spectrum{itn}(6:max_eig_vals), 'b');
      hold on
      plot(1:6, results.spectrum{itn}(1:6), 'ro')
      hold off
   elseif itn == 3
      max_eig_vals = length(results.spectrum{itn});
      plot(2:max_eig_vals, results.spectrum{itn}(2:max_eig_vals), 'b');
      hold on
      plot(1:2, results.spectrum{itn}(1:2), 'ro')
      hold off
   else
      plot(results.spectrum{itn}, 'b');
   end
%}   
   %strTitle = strcat('A(', num2str(itn), ')');
   %strTitle = strcat('A_{', num2str(iter-1), '}');
   strTitle = strcat("EMEP iteration = ", num2str(iter));
   title(strTitle);
   xlabel('Eigenvalue index');
   ylabel('Eigenvalue');
   xMin = -50; xMax = n + 50;
   yMin = -2;%min(results.spectrum{itn});
   yMax = 2;%max(results.spectrum{itn});
   axis([xMin xMax yMin yMax]);
end




figure
for i = 1:length(iters)
   iter = iters(i);
   subplot(2, length(iters)/2, i);
   spec = vec(results.spectrum{iter});
   sub_spec_max = 20;
   plot(spec(1:sub_spec_max), 'b*');
%{
   if itn == 1
      max_eig_vals = length(results.spectrum{itn});
      plot(6:max_eig_vals, results.spectrum{itn}(6:max_eig_vals), 'b');
      hold on
      plot(1:6, results.spectrum{itn}(1:6), 'ro')
      hold off
   elseif itn == 3
      max_eig_vals = length(results.spectrum{itn});
      plot(2:max_eig_vals, results.spectrum{itn}(2:max_eig_vals), 'b');
      hold on
      plot(1:2, results.spectrum{itn}(1:2), 'ro')
      hold off
   else
      plot(results.spectrum{itn}, 'b');
   end
%}   
   %strTitle = strcat('A_{', num2str(iter-1), '}');
   strTitle = strcat("EMEP iteration = ", num2str(iter));
   title(strTitle);
   xlabel('Eigenvalue index');
   ylabel('Eigenvalue');
   xMin = 0; xMax = sub_spec_max;
   yMin = min(spec(1:sub_spec_max));
   yMax = max(spec(1:sub_spec_max));
   shift = 0.1*(yMax - yMin);
   yMin = yMin - shift;
   yMax = yMax + shift;
   axis([xMin xMax yMin yMax]);
end



end % end function recreatefigure_noisyimage_spectrumdist



