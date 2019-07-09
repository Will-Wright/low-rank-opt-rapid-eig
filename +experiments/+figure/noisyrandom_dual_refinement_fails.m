function [results] = noisyrandom_dual_refinement_fails(varargin)
% experiments.figure.noisyrandom_dual_refinement_fails
%
% Generates table with number of DFTs and runtimes for recovery of 
% noisy random gaussian signal with and without refinement subroutine.
% Also generates figure 

ip = inputParser;

ip.addParameter('noise_ratios', [0.05, 0.15, 0.3]);
ip.addParameter('iterations', 100);
ip.parse(varargin{:});
noise_ratios = ip.Results.noise_ratios;
iterations = ip.Results.iterations;

n = 128;
L = 10;


folder_name = 'cache/figure_noisyrandom_dual_refinement_fails';
if exist(folder_name) ~= 7
   mkdir(folder_name);
end
exists_default_experiment = exist(strcat(folder_name, '/results.mat'));
if exists_default_experiment
	load(strcat(folder_name, '/results.mat'));
else
	results = [];
end

if isempty(results)
   i = 1;
   for noise_ratio = noise_ratios
      results{i}.DFP = experiments.evolvingmatricespl('L', L, 'n', n, ...
                                             'signal', 'gaussian', 'scale', [], ...
                                             'noise_ratio', noise_ratio, 'iterations', iterations, ...
                                             'verbosity', 1, 'logIterData', false, ...
                                             'printRelErr', true, ...
                                             'skipDFP', false, 'StopOnFeasible', false, ...
                                             'returnEigsIterData', true, ...
                                             'optTol', 0, ...
                                             'StopOnRelErr', false);
      results{i}.noDFP = experiments.evolvingmatricespl('L', L, 'n', n, ...
                                                'signal', 'gaussian', 'scale', [], ...
                                                'noise_ratio', noise_ratio, 'iterations', iterations, ...
                                                'verbosity', 1, 'logIterData', false, ...
                                                'printRelErr', true, ...
                                                'skipDFP', true, 'StopOnFeasible', false, ...
                                                'returnEigsIterData', true, ...
                                                'optTol', 0, ...
                                                'StopOnRelErr', false);
      i = i + 1;
   end
   save(strcat(folder_name, '/results.mat'), 'results');
end



for i = 1:length(noise_ratios)
   x_re_DFP = results{i}.DFP.saga_sd.solInfo.iterData.x_error_rel;
   x_re_noDFP = results{i}.noDFP.saga_sd.solInfo.iterData.x_error_rel;
   
   subplot(1, length(noise_ratios), i);
   plot(x_re_DFP, 'b');
   hold on
   plot(x_re_noDFP, 'r--');
   strTitle = strcat("Noise ratio = ", num2str(noise_ratios(i)));
   title(strTitle);
   ylabel('Signal relative error');
   xlabel('Iteration');
   xMin = 1; xMax = iterations + 1;
   yMin = .9*min([x_re_DFP, x_re_noDFP]);
   yMax = 1.1*max([x_re_DFP, x_re_noDFP]);
   axis([xMin xMax yMin yMax]);
end

end