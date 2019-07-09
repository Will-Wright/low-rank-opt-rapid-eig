function [results, x_rel_err] = noisyimage_signal_relerr_various_saga_iterates(varargin)
% experiments.figure.noisyimage_signal_relerr_various_saga_iterates
%
% Solves noisy image recovery problem and displays figures from 
% various iterates.  Also plots signal relative error.
% Demonstrates that iter 20 is essentially optimal.
%
% Note: solves RGB model, not grayscale model.


ip = inputParser;
ip.addParameter('results', []);
ip.addParameter('run_default_experiment', true); % must be set to `false` to load `results` struct
ip.addParameter('run_TeX_print_mode', true); % `false` prints single figure, `true` prints multiple
ip.addParameter('wflowOn', false);
ip.addParameter('image_file', 'data/parrot_16k.jpg');
ip.addParameter('noise_ratio', 0.30);
ip.addParameter('L', 8);
ip.addParameter('resizeImage', 1/1);
ip.addParameter('isGrayscale', false);
ip.addParameter('noise_type', 'gaussian');
ip.addParameter('iterations', 100);
ip.addParameter('StopOnRelErr', false); 
ip.addParameter('optTol', 0); 
ip.addParameter('feaTol', 0);
ip.addParameter('iterates_to_view', [1, 3, 21, 101]);
ip.parse(varargin{:});

results = ip.Results.results;
saga_iterations = ip.Results.iterations;
run_default_experiment = ip.Results.run_default_experiment;
run_TeX_print_mode = ip.Results.run_TeX_print_mode;
noise_ratio = ip.Results.noise_ratio;
L = ip.Results.L;
isGrayscale = ip.Results.isGrayscale;
iterates_to_view = ip.Results.iterates_to_view;
if run_default_experiment
   isGrayscale = false;
   noise_ratio = 0.3;
   L = 8;
   optsNoisyModel = struct('wflowOn', false, ... 
                           'L', L, 'noise_ratio', noise_ratio, ...
                           'noise_type', 'gaussian', ...
                           'resizeImage', 1/1, ...
                           'image_file', 'data/parrot_16k.jpg', ...
                           'iterations', saga_iterations, ...
                           'optTol', 0, ...
                           'feaTol', 0, ...
                           'StopOnRelErr', false, ...
                           'cutDim', 5, 'eigs_basis_size', 100);
   folder_name = 'cache/figure_noisyimage_signal_relerr_various_saga_iterates';
   if exist(folder_name) ~= 7
      mkdir(folder_name);
   end
   exists_default_experiment = exist(strcat(folder_name, '/results.mat'));
   if exists_default_experiment
      load(strcat(folder_name, '/results.mat'));
   else
      results = [];
   end
elseif isempty(results)
   optsNoisyModel = struct('wflowOn', ip.Results.wflowOn, ... 
                        'L', ip.Results.L, ...
                        'noise_type', ip.Results.noise_type, ...
                        'resizeImage', ip.Results.resizeImage, ...
                        'image_file', ip.Results.image_file, ...
                        'iterations', ip.Results.iterations, ...
                        'optTol', ip.Results.optTol, ...
                        'feaTol', ip.Results.feaTol, ...
                        'StopOnRelErr', ip.Results.StopOnRelErr, ...
                        'cutDim', 5, 'eigs_basis_size', 100);
end
num_iters_to_view = length(iterates_to_view);



if isempty(results)
   results = cell(0);
   optsNoisyModelTemp = optsNoisyModel;
   optsNoisyModelTemp.noise_ratio = noise_ratio;
   fprintf('Solving model with %f noise ratio.\n', noise_ratio);
   if isGrayscale
      results = experiments.evolvingmatricespl('channel', 'gray', optsNoisyModelTemp);
   else
      results = experiments.evolvingmatricespl('channel', 'RGB', optsNoisyModelTemp);
   end
end



if run_TeX_print_mode
   plot_multiple_figures(results, noise_ratio, isGrayscale, ...
      iterates_to_view, num_iters_to_view, saga_iterations)
else
   plot_single_figure(results, noise_ratio, isGrayscale, ...
      iterates_to_view, num_iters_to_view, saga_iterations)
end
 


% Saves default experiment if didn't previously exist.
if run_default_experiment && ~exists_default_experiment
   % Removes unnecessary data
   results.red.saga_sd.solInfo.iterData = rmfield(results.red.saga_sd.solInfo.iterData, 'y');
   results.green.saga_sd.solInfo.iterData = rmfield(results.green.saga_sd.solInfo.iterData, 'y');
   results.blue.saga_sd.solInfo.iterData = rmfield(results.blue.saga_sd.solInfo.iterData, 'y');
   results.red.saga_sd.solInfo.iterData = rmfield(results.red.saga_sd.solInfo.iterData, 'yLinesearch');
   results.green.saga_sd.solInfo.iterData = rmfield(results.green.saga_sd.solInfo.iterData, 'yLinesearch');
   results.blue.saga_sd.solInfo.iterData = rmfield(results.blue.saga_sd.solInfo.iterData, 'yLinesearch');
   results.red.saga_sd.solInfo.iterData = rmfield(results.red.saga_sd.solInfo.iterData, 'yRefined');
   results.green.saga_sd.solInfo.iterData = rmfield(results.green.saga_sd.solInfo.iterData, 'yRefined');
   results.blue.saga_sd.solInfo.iterData = rmfield(results.blue.saga_sd.solInfo.iterData, 'yRefined');
   for iter = 1:101
      save_x_iterate = sum(iter == iterates_to_view);
      if ~save_x_iterate
         results.red.saga_sd.solInfo.iterData.x{iter} = [];
         results.green.saga_sd.solInfo.iterData.x{iter} = [];
         results.blue.saga_sd.solInfo.iterData.x{iter} = [];         
      end
   end
   
   save(strcat(folder_name, '/results.mat'), 'results');
end
   

end





function [] = plot_single_figure(results, ...
   noise_ratio, isGrayscale, iterates_to_view, ...
   num_iters_to_view, saga_iterations)

figure;
hold on;
% Displays images of saga_sd results for various iterates
% Plots original image
if isGrayscale
   x_rel_err = results.saga_sd.solInfo.iterData.x_error_rel;
   n1 = results.gendatapl.A.n1;
   n2 = results.gendatapl.A.n2;
   imSize = [n1 n2];
   imOrig = real(reshape(results.gendatapl.x0, imSize));
   subplot(3, num_iters_to_view + 1, 1);
   image(imOrig);
   title('Original image');
   set(gca,'YTickLabel',[]); set(gca,'XTickLabel',[]);
   maxIt = results.saga_sd.solInfo.nit;
else
   x_rel_err = 1/3*(results.red.saga_sd.solInfo.iterData.x_error_rel ...
      + results.green.saga_sd.solInfo.iterData.x_error_rel ...
      + results.blue.saga_sd.solInfo.iterData.x_error_rel);
   n1 = results.red.gendatapl.A.n1;
   n2 = results.red.gendatapl.A.n2;
   imSize = [n1 n2];
   imOrig = cat(3, ...
      reshape(results.red.gendatapl.x0, imSize), ...
      reshape(results.green.gendatapl.x0, imSize), ...
      reshape(results.blue.gendatapl.x0, imSize));
   subplot(3, num_iters_to_view + 1, 1);
   image(imOrig);
   title('Original image');
   set(gca,'YTickLabel',[]); set(gca,'XTickLabel',[]);
   maxIt = min([results.red.saga_sd.solInfo.nit, ...
      results.green.saga_sd.solInfo.nit, ...
      results.blue.saga_sd.solInfo.nit]);
end

% Plots images of selected iterates
for iter_image = 1:num_iters_to_view
   iter_saga = iterates_to_view(iter_image);   
   if isGrayscale
      % Identifies correct sign of result 'x' and removes negative pixel values
      realX = real(reshape(results.saga_sd.solInfo.iterData.x{iter_saga}, imSize));
      idxMatPos = realX > 0;
      xSagaMat = zeros(imSize);
      if sum(sum(idxMatPos)) >= 0.5*results.gendatapl.A.n
         xSagaMat(idxMatPos) = realX(idxMatPos);
      else
         xSagaMat(~idxMatPos) = -realX(~idxMatPos);
      end
      imSaga = mat2gray(xSagaMat);

      subplot(3, num_iters_to_view + 1, 1 + iter_image);
      imshow(imSaga); title('Recovered with Algorithm 3');

   else
      % Identifies correct sign of result 'x' and removes negative pixel values
      realX = real(reshape(results.red.saga_sd.solInfo.iterData.x{iter_saga}, imSize));
      idxMatPos = realX > 0;
      xSagaMatRed = zeros(imSize);
      if sum(sum(idxMatPos)) >= 0.5*results.red.gendatapl.A.n
         xSagaMatRed(idxMatPos) = realX(idxMatPos);
      else
         xSagaMatRed(~idxMatPos) = -realX(~idxMatPos);
      end
      realX = real(reshape(results.green.saga_sd.solInfo.iterData.x{iter_saga}, imSize));
      idxMatPos = realX > 0;
      xSagaMatGreen = zeros(imSize);
      if sum(sum(idxMatPos)) >= 0.5*results.green.gendatapl.A.n
         xSagaMatGreen(idxMatPos) = realX(idxMatPos);
      else
         xSagaMatGreen(~idxMatPos) = -realX(~idxMatPos);
      end
      realX = real(reshape(results.blue.saga_sd.solInfo.iterData.x{iter_saga}, imSize));
      idxMatPos = realX > 0;
      xSagaMatBlue = zeros(imSize);
      if sum(sum(idxMatPos)) >= 0.5*results.blue.gendatapl.A.n
         xSagaMatBlue(idxMatPos) = realX(idxMatPos);
      else
         xSagaMatBlue(~idxMatPos) = -realX(~idxMatPos);
      end

      imSaga = cat(3, xSagaMatRed, xSagaMatGreen, xSagaMatBlue); 

      norm_b_red = norm(results.red.gendatapl.b(:));
      norm_b_green = norm(results.green.gendatapl.b(:));
      norm_b_blue = norm(results.blue.gendatapl.b(:));
      pr_rel_err = 1/3*(results.red.saga_sd.solInfo.iterData.res_primal / norm_b_red ...
         + results.green.saga_sd.solInfo.iterData.res_primal / norm_b_green ...
         + results.blue.saga_sd.solInfo.iterData.res_primal / norm_b_blue);

      subplot(3, num_iters_to_view + 1, iter_image + 1);
      image(imSaga); title(['Iterate ' num2str(iter_saga-1)]);
      set(gca,'YTickLabel',[]); set(gca,'XTickLabel',[]);
      xlabel(['obs rel err ' num2str(pr_rel_err(iter_saga))]);
   end
end

% Plots observation rel err |Ax-b| / |b|
subplot(3, num_iters_to_view + 1, (num_iters_to_view+1) + 1 ...
                                    : 2*(num_iters_to_view+1) );
plot(0:maxIt, pr_rel_err, 'b');
hold on
plot(noise_ratio*ones(saga_iterations + 1, 1) , 'r');
title('Primal relative error');
ylabel('Relative error');
xlabel('Iteration');
xMin = 0; xMax = maxIt;
yMin = 0.9*min(pr_rel_err);
yMax = 1.1*max(pr_rel_err);
axis([xMin xMax yMin yMax]);
for iter_image = 1:num_iters_to_view
   iter_saga = iterates_to_view(iter_image);
   plot(iter_saga - 1, pr_rel_err(iter_saga), 'or');
end

% Plots signal rel err |x - x_star| / |x_star|
subplot(3, num_iters_to_view+1, 2*(num_iters_to_view+1) + 1 ...
                                  : 3*(num_iters_to_view+1) );
plot(0:maxIt, x_rel_err, 'b');
hold on
title('Signal relative error');
ylabel('Relative error');
xlabel('Iteration');
xMin = 0; xMax = maxIt;
yMin = 0.9*min(x_rel_err);
yMax = 1.1*max(x_rel_err);
axis([xMin xMax yMin yMax]);
for iter_image = 1:num_iters_to_view
   iter_saga = iterates_to_view(iter_image);
   plot(iter_saga - 1, x_rel_err(iter_saga), 'or');
end     

end




function [] = plot_multiple_figures(results, ...
   noise_ratio, isGrayscale, iterates_to_view, ...
   num_iters_to_view, saga_iterations)

% Displays figure 1 of 3
figure;
hold on;
% Displays images of saga_sd results for various iterates
% Plots original image
if isGrayscale
   x_rel_err = results.saga_sd.solInfo.iterData.x_error_rel;
   n1 = results.gendatapl.A.n1;
   n2 = results.gendatapl.A.n2;
   imSize = [n1 n2];
   imOrig = real(reshape(results.gendatapl.x0, imSize));
   subplot(1, num_iters_to_view + 1, 1);
   image(imOrig);
   title('Original image');
   set(gca,'YTickLabel',[]); set(gca,'XTickLabel',[]);
   maxIt = results.saga_sd.solInfo.nit;
else
   x_rel_err = 1/3*(results.red.saga_sd.solInfo.iterData.x_error_rel ...
      + results.green.saga_sd.solInfo.iterData.x_error_rel ...
      + results.blue.saga_sd.solInfo.iterData.x_error_rel);
   n1 = results.red.gendatapl.A.n1;
   n2 = results.red.gendatapl.A.n2;
   imSize = [n1 n2];
   imOrig = cat(3, ...
      reshape(results.red.gendatapl.x0, imSize), ...
      reshape(results.green.gendatapl.x0, imSize), ...
      reshape(results.blue.gendatapl.x0, imSize));
   subplot(1, num_iters_to_view + 1, 1);
   image(imOrig);
   title('Original image');
   set(gca,'YTickLabel',[]); set(gca,'XTickLabel',[]);
   maxIt = min([results.red.saga_sd.solInfo.nit, ...
      results.green.saga_sd.solInfo.nit, ...
      results.blue.saga_sd.solInfo.nit]);
end

% Plots images of selected iterates
for iter_image = 1:num_iters_to_view
   iter_saga = iterates_to_view(iter_image);   
   if isGrayscale
      % Identifies correct sign of result 'x' and removes negative pixel values
      realX = real(reshape(results.saga_sd.solInfo.iterData.x{iter_saga}, imSize));
      idxMatPos = realX > 0;
      xSagaMat = zeros(imSize);
      if sum(sum(idxMatPos)) >= 0.5*results.gendatapl.A.n
         xSagaMat(idxMatPos) = realX(idxMatPos);
      else
         xSagaMat(~idxMatPos) = -realX(~idxMatPos);
      end
      imSaga = mat2gray(xSagaMat);

      subplot(1, num_iters_to_view + 1, 1 + iter_image);
      imshow(imSaga); title('Recovered with Algorithm 3');

   else
      % Identifies correct sign of result 'x' and removes negative pixel values
      realX = real(reshape(results.red.saga_sd.solInfo.iterData.x{iter_saga}, imSize));
      idxMatPos = realX > 0;
      xSagaMatRed = zeros(imSize);
      if sum(sum(idxMatPos)) >= 0.5*results.red.gendatapl.A.n
         xSagaMatRed(idxMatPos) = realX(idxMatPos);
      else
         xSagaMatRed(~idxMatPos) = -realX(~idxMatPos);
      end
      realX = real(reshape(results.green.saga_sd.solInfo.iterData.x{iter_saga}, imSize));
      idxMatPos = realX > 0;
      xSagaMatGreen = zeros(imSize);
      if sum(sum(idxMatPos)) >= 0.5*results.green.gendatapl.A.n
         xSagaMatGreen(idxMatPos) = realX(idxMatPos);
      else
         xSagaMatGreen(~idxMatPos) = -realX(~idxMatPos);
      end
      realX = real(reshape(results.blue.saga_sd.solInfo.iterData.x{iter_saga}, imSize));
      idxMatPos = realX > 0;
      xSagaMatBlue = zeros(imSize);
      if sum(sum(idxMatPos)) >= 0.5*results.blue.gendatapl.A.n
         xSagaMatBlue(idxMatPos) = realX(idxMatPos);
      else
         xSagaMatBlue(~idxMatPos) = -realX(~idxMatPos);
      end

      imSaga = cat(3, xSagaMatRed, xSagaMatGreen, xSagaMatBlue); 

      norm_b_red = norm(results.red.gendatapl.b(:));
      norm_b_green = norm(results.green.gendatapl.b(:));
      norm_b_blue = norm(results.blue.gendatapl.b(:));

      subplot(1, num_iters_to_view + 1, iter_image + 1);
      image(imSaga); title(['Iterate ' num2str(iter_saga-1)]);
      set(gca,'YTickLabel',[]); set(gca,'XTickLabel',[]);
   end
end





% Displays figure 2 of 3
figure;
hold on;
% Displays images of saga_sd results for various iterates
% Plots original image
if isGrayscale
   x_rel_err = results.saga_sd.solInfo.iterData.x_error_rel;
   n1 = results.gendatapl.A.n1;
   n2 = results.gendatapl.A.n2;
   imSize = [n1 n2];
   imOrig = real(reshape(results.gendatapl.x0, imSize));
   subplot(2, num_iters_to_view + 1, 1);
   image(imOrig);
   title('Original image');
   set(gca,'YTickLabel',[]); set(gca,'XTickLabel',[]);
   maxIt = results.saga_sd.solInfo.nit;
else
   x_rel_err = 1/3*(results.red.saga_sd.solInfo.iterData.x_error_rel ...
      + results.green.saga_sd.solInfo.iterData.x_error_rel ...
      + results.blue.saga_sd.solInfo.iterData.x_error_rel);
   norm_b_red = norm(results.red.gendatapl.b(:));
   norm_b_green = norm(results.green.gendatapl.b(:));
   norm_b_blue = norm(results.blue.gendatapl.b(:));
   pr_rel_err = 1/3*(results.red.saga_sd.solInfo.iterData.res_primal / norm_b_red ...
      + results.green.saga_sd.solInfo.iterData.res_primal / norm_b_green ...
      + results.blue.saga_sd.solInfo.iterData.res_primal / norm_b_blue);
   n1 = results.red.gendatapl.A.n1;
   n2 = results.red.gendatapl.A.n2;
   imSize = [n1 n2];
   imOrig = cat(3, ...
      reshape(results.red.gendatapl.x0, imSize), ...
      reshape(results.green.gendatapl.x0, imSize), ...
      reshape(results.blue.gendatapl.x0, imSize));
   subplot(2, num_iters_to_view + 1, 1);
   image(imOrig);
   title('Original image');
   set(gca,'YTickLabel',[]); set(gca,'XTickLabel',[]);
   maxIt = min([results.red.saga_sd.solInfo.nit, ...
      results.green.saga_sd.solInfo.nit, ...
      results.blue.saga_sd.solInfo.nit]);
end

% Plots images of selected iterates
for iter_image = 1:num_iters_to_view
   iter_saga = iterates_to_view(iter_image);   
   if isGrayscale
      % Identifies correct sign of result 'x' and removes negative pixel values
      realX = real(reshape(results.saga_sd.solInfo.iterData.x{iter_saga}, imSize));
      idxMatPos = realX > 0;
      xSagaMat = zeros(imSize);
      if sum(sum(idxMatPos)) >= 0.5*results.gendatapl.A.n
         xSagaMat(idxMatPos) = realX(idxMatPos);
      else
         xSagaMat(~idxMatPos) = -realX(~idxMatPos);
      end
      imSaga = mat2gray(xSagaMat);

      subplot(2, num_iters_to_view + 1, 1 + iter_image);
      imshow(imSaga); title('Recovered with Algorithm 3');

   else
      % Identifies correct sign of result 'x' and removes negative pixel values
      realX = real(reshape(results.red.saga_sd.solInfo.iterData.x{iter_saga}, imSize));
      idxMatPos = realX > 0;
      xSagaMatRed = zeros(imSize);
      if sum(sum(idxMatPos)) >= 0.5*results.red.gendatapl.A.n
         xSagaMatRed(idxMatPos) = realX(idxMatPos);
      else
         xSagaMatRed(~idxMatPos) = -realX(~idxMatPos);
      end
      realX = real(reshape(results.green.saga_sd.solInfo.iterData.x{iter_saga}, imSize));
      idxMatPos = realX > 0;
      xSagaMatGreen = zeros(imSize);
      if sum(sum(idxMatPos)) >= 0.5*results.green.gendatapl.A.n
         xSagaMatGreen(idxMatPos) = realX(idxMatPos);
      else
         xSagaMatGreen(~idxMatPos) = -realX(~idxMatPos);
      end
      realX = real(reshape(results.blue.saga_sd.solInfo.iterData.x{iter_saga}, imSize));
      idxMatPos = realX > 0;
      xSagaMatBlue = zeros(imSize);
      if sum(sum(idxMatPos)) >= 0.5*results.blue.gendatapl.A.n
         xSagaMatBlue(idxMatPos) = realX(idxMatPos);
      else
         xSagaMatBlue(~idxMatPos) = -realX(~idxMatPos);
      end

      imSaga = cat(3, xSagaMatRed, xSagaMatGreen, xSagaMatBlue); 

      subplot(2, num_iters_to_view + 1, iter_image + 1);
      image(imSaga); title(['Iterate ' num2str(iter_saga-1)]);
      set(gca,'YTickLabel',[]); set(gca,'XTickLabel',[]);
%      xlabel(['SRE = ' num2str(x_rel_err(iter_saga))]);
   end
end

% Plots signal rel err |x - x_star| / |x_star|
subplot(2, num_iters_to_view+1, 1*(num_iters_to_view+1) + 1 ...
                                  : 2*(num_iters_to_view+1) );
plot(0:maxIt, x_rel_err, 'b');
hold on
title('Signal relative error');
ylabel('Relative error');
xlabel('Iteration');
xMin = 0; xMax = maxIt;
yMin = 0.95*min(x_rel_err);
yMax = 1.05*max(x_rel_err);
axis([xMin xMax yMin yMax]);
for iter_image = 1:num_iters_to_view
   iter_saga = iterates_to_view(iter_image);
   plot(iter_saga - 1, x_rel_err(iter_saga), 'or');
end 








% Displays figure 3 of 3
figure;
hold on;
% Plots signal and observation rel errs

% First plots signal rel err |x - x_star| / |x_star|
%subplot(2, num_iters_to_view+1, 0*(num_iters_to_view+1) + 1 ...
%                                  : 1*(num_iters_to_view+1) );
subplot(1, 2, 1);
plot(0:maxIt, x_rel_err, 'b');
hold on
title('Signal relative error');
ylabel('Relative error');
xlabel('Iteration');
xMin = 0; xMax = maxIt;
yMin = 0.95*min(x_rel_err);
yMax = 1.05*max(x_rel_err);
axis([xMin xMax yMin yMax]);
for iter_image = 1:num_iters_to_view
   iter_saga = iterates_to_view(iter_image);
   plot(iter_saga - 1, x_rel_err(iter_saga), 'or');
end
%NumTicks = 5;
%L = get(gca,'YLim');
%set(gca,'YTick',linspace(0.2, 0.4, NumTicks))

% Next plots observation rel err |Ax-b| / |b|
%subplot(2, num_iters_to_view + 1, 1*(num_iters_to_view+1) + 1 ...
%                                    : 2*(num_iters_to_view+1) );
subplot(1, 2, 2);
plot(0:maxIt, pr_rel_err, 'b');
hold on
plot(noise_ratio*ones(saga_iterations + 1, 1) , 'r');
title('Primal relative error');
ylabel('Relative error');
xlabel('Iteration');
xMin = 0; xMax = maxIt;
yMin = 0.95*min(pr_rel_err);
yMax = 1.05*max(pr_rel_err);
axis([xMin xMax yMin yMax]);
for iter_image = 1:num_iters_to_view
   iter_saga = iterates_to_view(iter_image);
   plot(iter_saga - 1, pr_rel_err(iter_saga), 'or');
end
%NumTicks = 5;
%L = get(gca,'YLim');
%set(gca,'YTick',linspace(0.25, 0.45, NumTicks))
hold off;


end

