function [results] = noisyrandom_relative_errors_stagnate(varargin)
% experiments.figure.noisyrandom_relative_errors_stagnate
%
% Generates 10000 iterates from saga_sd of random noisy phase retrieval problem
% and creates figure with relative errors.

% noisy_random_signal_relative_errors_stagnate

ip = inputParser;
ip.addParameter('results', []);
ip.addParameter('iterations', 100);
ip.parse(varargin{:});

results = ip.Results.results;
iterations = ip.Results.iterations;

verbosity = 1;
n = 16;
L = 6;
noise_ratio = 0.3;
feaTol = 0;
optTol = 0;
feaTol_orig = 2e-4; % used in plot
dualRelErrTol = 0;

folder_name = 'cache/figure_noisyrandom_relative_errors_stagnate';
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
   [results, data] = experiments.exputil.solve_noisyimage_saga_vs_exact_SDP_solutions('iterations', iterations, ...
   'n', n, 'L', L, 'noise_ratio', noise_ratio, 'signal', 'gaussian', ...
   'feaTol', feaTol, 'optTol', optTol, 'dualRelErrTol', dualRelErrTol, ...
   'verbosity', verbosity);
   y_star = data.y_cvx_dual(:);
   data.y_rel_err = zeros(iterations, 1);
   for i = 1:iterations + 1
      data.y_rel_err(i) = norm(results.saga_sd.solInfo.iterData.y{i}(:) - y_star) / norm(y_star);
   end
   results.data = data;
   results.saga_sd.solInfo.iterData = rmfield(results.saga_sd.solInfo.iterData, 'x');
   results.saga_sd.solInfo.iterData = rmfield(results.saga_sd.solInfo.iterData, 'y');
   results.saga_sd.solInfo.iterData = rmfield(results.saga_sd.solInfo.iterData, 'yLinesearch');
   save(strcat(folder_name, '/results.mat'), 'results');
else
   data = results.data;
end




% Plots relative error of signal obseration A(x_k) against true observation
subplot(1, 3, 1);
norm_b = norm(results.gendatapl.b(:));
semilogx(results.saga_sd.solInfo.iterData.res_primal / norm_b, 'b');
hold on
%plot(noise_ratio*ones(results.saga_sd.solInfo.nit + 1, 1), 'r');
plot( (noise_ratio + feaTol_orig*(1 + 1/norm_b))*ones(results.saga_sd.solInfo.nit + 1, 1), 'r');
hold off
title('Primal relative error');
ylabel('Relative error');
xlabel('Iteration');
xMin = 0; xMax = iterations + 1;
yMin = min(results.saga_sd.solInfo.iterData.res_primal / norm_b) - 0.1;
yMax = max(results.saga_sd.solInfo.iterData.res_primal / norm_b) + 0.1;
axis([xMin xMax yMin yMax]);
set(gca,'XTick',[1, logspace(1,4,4)])


% Plots relative error of dual iterate y_k against y_star obtained via cvx
% (SDP solver with residuals around 1e-8)
subplot(1, 3, 2);
semilogx(data.y_rel_err, 'b');
title('Dual relative error');
ylabel('Relative error');
xlabel('Iteration');
xMin = 0; xMax = iterations + 1;
yMin = 0;
yMax = max(data.y_rel_err) + 100;
axis([xMin xMax yMin yMax]);
set(gca,'XTick',[1, logspace(1,4,4)])

% Plots strong duality gap
subplot(1, 3, 3);
semilogx(results.saga_sd.solInfo.iterData.res_mGap, 'b');
title('Duality gap');
ylabel('Error');
xlabel('Iteration');
xMin = 0; xMax = iterations + 1;
yMin = 0;
yMax = max(results.saga_sd.solInfo.iterData.res_mGap) + 100;
axis([xMin xMax yMin yMax]);
set(gca,'XTick',[1, logspace(1,4,4)])


%fprintf('Press enter to continue.\n');
%pause;
%experiments.exputil.view_space_of_dual_solutions('results', results, 'data', data);

end