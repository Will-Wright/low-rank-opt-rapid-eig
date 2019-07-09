function [results, data] = noisyimage_new_termination_criteria(varargin)
% experiments.figure.noisyimage_new_termination_criteria
%
% Solves noisy image recovery problem and displays figures and table
% indicating best residuals are y_diff and pr_diff, and du_gap to pick
% best iterate of previous 20.
%
% Note: Solves grayscale model.

iter_max = 1000;
L_range = [5, 10];
noise_ratio_range = [0.05, 0.15, 0.3];
print_figs_and_tables = true;
figure_1_on = true;
figure_2_on = true;
figure_3_on = true;
figure_4_on = true;
figure_5_on = true;
figure_6_on = true;

% Sets parameters for experiment
optsNoisyModel = struct('wflowOn', false, ... 
                        'image_file', 'data/parrot_4k.jpg', ...
                        'iterations', iter_max, 'maxTime', 10*60*60, 'optTol', 0, ...
                        'StopOnRelErr', false, 'resizeImage', 1/1, ...
                        'eigs_basis_size', 200, 'cutDim', 50, 'skipLS', false, ...
                        'skipBB', false, 'logIterData', false, 'eigIts', 10000);
                     
folder_name = 'cache/figure_noisyimage_new_termination_criteria';
if exist(folder_name) ~= 7
   mkdir(folder_name);
end
exists_default_experiment = exist(strcat(folder_name, '/results.mat'));
if exists_default_experiment
      load(strcat(folder_name, '/results.mat'));
else
      results = [];
end

% Runs sequence of experiments for various oversampling (L) and noise ratio
if isempty(results)
   L_idx = 1;
   noise_ratio_idx = 1;
   results = cell(1,1);
   for L = L_range
      for noise_ratio = noise_ratio_range
         results{L_idx, noise_ratio_idx} = experiments.evolvingmatricespl(...
            optsNoisyModel, 'L', L, 'noise_ratio', noise_ratio);
         save(strcat(folder_name, '/results.mat'), 'results');
         noise_ratio_idx = noise_ratio_idx + 1;
      end
      noise_ratio_idx = 1;
      L_idx = L_idx + 1;
   end
end


if print_figs_and_tables

   % Prints results tables and figures
   sig_rel_err = []; sig_diff = []; du_gap = []; du_gap_diff = []; 
   du_obj = []; du_obj_diff = []; y_diff = []; stepLS = [];
   pr_rel_err = []; pr_diff = [];
   iter_start = 2;
   iter_max = 1001;

   for L_idx = 1:2
      for noise_ratio_idx = 1:3
         for iter = iter_start:iter_max
            norm_b = norm(results{L_idx,noise_ratio_idx}.gendatapl.b(:));
            
            sig_rel_err(L_idx, noise_ratio_idx, iter) = abs(results{L_idx,noise_ratio_idx}.saga_sd.solInfo.iterData.x_error_rel(iter));
            sig_diff(L_idx, noise_ratio_idx, iter) = abs(results{L_idx,noise_ratio_idx}.saga_sd.solInfo.iterData.x_error_rel(iter) - results{L_idx,noise_ratio_idx}.saga_sd.solInfo.iterData.x_error_rel(iter -1)) ...
               / abs(results{L_idx,noise_ratio_idx}.saga_sd.solInfo.iterData.x_error_rel(iter));

            stepLS(L_idx, noise_ratio_idx, iter) = cell2mat(results{L_idx,noise_ratio_idx}.saga_sd.solInfo.iterData.stepLS(iter));

            du_gap(L_idx, noise_ratio_idx, iter) = results{L_idx,noise_ratio_idx}.saga_sd.solInfo.iterData.res_mGap(iter);
            du_gap_diff(L_idx, noise_ratio_idx, iter) = abs(results{L_idx,noise_ratio_idx}.saga_sd.solInfo.iterData.res_mGap(iter) - results{L_idx,noise_ratio_idx}.saga_sd.solInfo.iterData.res_mGap(iter-1)) ...
               / abs(results{L_idx,noise_ratio_idx}.saga_sd.solInfo.iterData.res_mGap(iter));

            pr_rel_err(L_idx, noise_ratio_idx, iter) = abs(results{L_idx,noise_ratio_idx}.saga_sd.solInfo.iterData.res_primal(iter)) / norm_b;
            pr_diff(L_idx, noise_ratio_idx, iter) = abs(results{L_idx,noise_ratio_idx}.saga_sd.solInfo.iterData.res_primal(iter) - results{L_idx,noise_ratio_idx}.saga_sd.solInfo.iterData.res_primal(iter-1)) ...
               / abs(results{L_idx,noise_ratio_idx}.saga_sd.solInfo.iterData.res_primal(iter));

            du_obj(L_idx, noise_ratio_idx, iter) = results{L_idx,noise_ratio_idx}.saga_sd.solInfo.iterData.d{1,iter}(1);
            du_obj_diff(L_idx, noise_ratio_idx, iter) = abs(results{L_idx,noise_ratio_idx}.saga_sd.solInfo.iterData.d{1,iter}(1) - results{L_idx,noise_ratio_idx}.saga_sd.solInfo.iterData.d{1,iter-1}(1)) ...
               / abs(results{L_idx,noise_ratio_idx}.saga_sd.solInfo.iterData.d{1,iter}(1));

            % basically same as A_k
            y_diff(L_idx, noise_ratio_idx, iter) = results{L_idx,noise_ratio_idx}.saga_sd.solInfo.iterData.res_duRelErr(iter);
         end
      end
   end

   % Identifies iterate at which linesearch cut off due to algo
   % determining obj fun was nondifferentiable
   modelNondiff = NaN(2,3);
   for L_idx = 1:2
      for noise_ratio_idx = 1:3
         for iter = 2:1001
            if cell2mat(results{L_idx,noise_ratio_idx}.saga_sd.solInfo.iterData.skipBB(iter)) == true
               modelNondiff(L_idx, noise_ratio_idx) = iter - 1;
               break;
            end
         end
      end
   end
   
   
   
   if figure_1_on
      % Plots figure of iterates vs signal relative error across models
      figure
      iter_start = 2;
      iter_max = 1001;
      i = 1;
      for noise_ratio_idx = 1:3
         for L_idx = 1:2
            L = L_range(L_idx);
            noise_ratio = noise_ratio_range(noise_ratio_idx);
            subplot(3, 2, i);
            seq_plot = vec(sig_rel_err(L_idx, noise_ratio_idx, iter_start:iter_max));
            plot(iter_start-1:iter_max-1, seq_plot, 'b')
            title_str = strcat("L = ", num2str(L), ", noise ratio = ", num2str(noise_ratio));
            title(title_str)
            xlabel('Iteration')
            ylabel('Relative Error')
            %xlim([iter_start - 2, iter_max-1]);
            plot_buffer = (max(seq_plot) - min(seq_plot)) / 20;
            ylim([min(seq_plot)- plot_buffer, max(seq_plot) + plot_buffer])
            % Proto code to control vertical axis ticks
%            format long
%            [min(seq_plot)- plot_buffer, max(seq_plot) + plot_buffer]
%            set(gca,'XLim',[103.6 104],'XTick',[103.6:0.10:104])
%            set(gca,'XLim',[1.25 1.50],'YTick',[1.25:0.05:1.50])
            % Plots point at which iterate becomes nondifferentiable
            hold on
            iter_nondiff = modelNondiff(L_idx, noise_ratio_idx);
            if ~isnan(iter_nondiff)
               plot(iter_nondiff - 1, seq_plot(iter_nondiff - 1), 'or');
            end
            hold off

%{            
            subplot(6, 2, 2*i);
            iter_start_temp = 101;
            seq_plot_temp = vec(sig_rel_err(L_idx, noise_ratio_idx, iter_start_temp:iter_max));
            plot(iter_start_temp-1:iter_max-1, seq_plot_temp, 'b')
            title_str = strcat('latter iterates: L = ', num2str(L), ', noise ratio =  ', num2str(noise_ratio));
            title(title_str)
            xlabel('iteration')
            ylabel('error')
            %xlim([iter_start_temp - 2, iter_max-1])

            plot_buffer = (max(seq_plot_temp) - min(seq_plot_temp)) / 20;
            xlim([iter_start_temp-1, iter_max-1]);
            ylim([min(seq_plot_temp)- plot_buffer, max(seq_plot_temp) + plot_buffer])
            % Plots point at which iterate becomes nondifferentiable
            hold on
            iter_nondiff = modelNondiff(L_idx, noise_ratio_idx);
            if ~isnan(iter_nondiff)
               plot(iter_nondiff - 1, seq_plot(iter_nondiff - 1), 'or');
            end
            hold off
%}
            i = i + 1;
         end
      end

      % Prints table of final duality gap values for models
      fprintf('\n')
      fprintf('Duality gap values for various models at 1000 iterations\n');
      fprintf('             |  L = %1i | L = %1i \n', L_range(1), L_range(2)');
      noise_ratio_idx = 1;
      fprintf('noise = %1.2f |  %5.2f |  %5.2f \n', noise_ratio_range(noise_ratio_idx), ...
         du_gap(1, noise_ratio_idx, iter_max), ...
         du_gap(2, noise_ratio_idx, iter_max));
      noise_ratio_idx = 2;
      fprintf('noise = %1.2f |  %5.2f |  %5.2f \n', noise_ratio_range(noise_ratio_idx), ...
         du_gap(1, noise_ratio_idx, iter_max), ...
         du_gap(2, noise_ratio_idx, iter_max));   
      noise_ratio_idx = 3;
      fprintf('noise = %1.2f |  %5.2f |  %5.2f \n', noise_ratio_range(noise_ratio_idx), ...
         du_gap(1, noise_ratio_idx, iter_max), ...
         du_gap(2, noise_ratio_idx, iter_max));

      % Prints table of iterate at which each model became primal feasible
      prFeasZero = NaN(L_idx, noise_ratio_idx);
      for L_idx = 1:2
         for noise_ratio_idx = 1:3
            noise_ratio = noise_ratio_range(noise_ratio_idx);
            val = find(vec(pr_rel_err(L_idx, noise_ratio_idx, iter_start:iter_max)) <= noise_ratio, 1, 'first');
            if ~isempty(val)
               prFeasZero(L_idx, noise_ratio_idx) = val(1,1);
            end
         end
      end
      fprintf('\n')
      fprintf('Iterate at which model became primal feasible for various models\n');
      fprintf('             |  L = %1i | L = %1i \n', L_range(1), L_range(2)');
      noise_ratio_idx = 1;
      fprintf('noise = %1.2f |   %3i  |   %3i  \n', noise_ratio_range(noise_ratio_idx), ...
         prFeasZero(1, noise_ratio_idx), ...
         prFeasZero(2, noise_ratio_idx));
      noise_ratio_idx = 2;
      fprintf('noise = %1.2f |   %3i  |   %3i  \n', noise_ratio_range(noise_ratio_idx), ...
         prFeasZero(1, noise_ratio_idx), ...
         prFeasZero(2, noise_ratio_idx));   
      noise_ratio_idx = 3;
      fprintf('noise = %1.2f |   %3i  |   %3i  \n', noise_ratio_range(noise_ratio_idx), ...
         prFeasZero(1, noise_ratio_idx), ...
         prFeasZero(2, noise_ratio_idx));
   end
   
   
   
   if figure_2_on
      % Plots primal residual vs iterates for model which never becomes
      % feasible
      figure;
      L_idx = 1; noise_ratio_idx = 2;
      L = L_range(L_idx);
      noise_ratio = noise_ratio_range(noise_ratio_idx);
      seq_plot = vec(pr_rel_err(L_idx, noise_ratio_idx, iter_start:iter_max));
      plot(iter_start-1:iter_max-1, seq_plot, 'b')
      hold on
      %seq_plot_prFea1 = (noise_ratio + 2e-4*(1+norm_b)/norm_b)*ones(iter_max - iter_start + 1, 1);
      seq_plot_prFea = (noise_ratio)*ones(iter_max - iter_start + 1, 1);
      plot(seq_plot_prFea, 'r')
      %title_str = strcat('L = ', num2str(L), ', noise ratio =  ', num2str(noise_ratio));
      title('Primal relative error and feasibility threshold')
      xlabel('Iteration')
      ylabel('Relative error')
      %xlim([iter_start-2, iter_max-1]);
      plot_buffer = (max(seq_plot) - min(seq_plot)) / 20;
      ylim([min([seq_plot; noise_ratio])- plot_buffer, max([seq_plot; noise_ratio]) + plot_buffer])
      hold off

      % Prints 
      norm_b = norm(results{L_idx,noise_ratio_idx}.gendatapl.b(:));
      prTol = noise_ratio + 2e-4*(1+norm_b)/norm_b;

      optsNoisyModel.wflowOn = true;
      optsNoisyModel.wflowUseOptimalStep = false;
      optsNoisyModel.saga_sdOn = false;
      optsNoisyModel.wflowVerbosity = false;
      results_wf = experiments.evolvingmatricespl(optsNoisyModel, 'L', 5, 'noise_ratio', 0.15);
      b = results_wf.gendatapl.b; 
      wflow_pr_rel_err = norm(vec(results_wf.gendatapl.A.forward(results_wf.wflow.x) - b)) / norm(b(:));

      fprintf('\n')
      fprintf('The following results were recovered for model with L = 5, noise_ratio = 0.15\n')
      fprintf('saga_sd primal error: %f\n', pr_rel_err(1, 2, end))
      fprintf('wflow primal error:   %f\n', wflow_pr_rel_err)
      fprintf('primal feas tol:      %f\n', prTol)
      fprintf('\n')
   end
    
   
   
   if figure_3_on
      
      % min and max iter vals for proper convergence window
      % for each model
      model_iter_min = [200, 200, 100; 50, 50, 50];
      model_iter_max = [400, 400, 200; 200, 100, 100];

      figure;
      % Plots iterate cutoffs for various primal rel err tolerances
      pr_range = 1e-3:-1e-7:1e-7;
      plot_idx = 0;
      for L_idx = 1:2
         for noise_ratio_idx = 1:3
            idx = 1;
            for tol = pr_range
               val = find(pr_diff(L_idx, noise_ratio_idx, 2:end) <= tol, 1);
               if ~isempty(val)
                  pr_tol_iter(L_idx, noise_ratio_idx, idx) = val;
               else
                  pr_tol_iter(L_idx, noise_ratio_idx, idx) = NaN;
               end
               idx = idx + 1;
            end
            subplot(6, 2, 1 + 2*plot_idx)
            semilogy(vec(pr_tol_iter(L_idx, noise_ratio_idx, :)), pr_range, 'b')
            hold on
            
            % plots tol vals
            for tol_temp = [1e-3, 1e-4, 1e-5, 1e-6, 1e-7]         
               pr_idx = find(vec(pr_diff(L_idx, noise_ratio_idx, 2:end)) <= tol_temp, 1);
               if ~isempty(pr_idx)
                  plot(pr_idx, tol_temp, 'or');
               end
            end
                        
            % plots window of termination
            iter_temp_min = model_iter_min(L_idx, noise_ratio_idx);
            iter_temp_max = model_iter_max(L_idx, noise_ratio_idx);
            x_points = [iter_temp_min, iter_temp_min, iter_temp_max, iter_temp_max];  
            y_points = [1e-7, 1e-3, 1e-3, 1e-7];
            color = [0, 0, 1];
            a = fill(x_points, y_points, color);
            a.FaceAlpha = 0.1;
            
            xlim([1, 1001])
            ylim([min(pr_range), max(pr_range)])
            set(gca,'YTick',flip(logspace(-3,-7,3)))
            L = L_range(L_idx);
            noise_ratio = noise_ratio_range(noise_ratio_idx);
            title_str = strcat("Primal: L = ", num2str(L), ", noise ratio = ", num2str(noise_ratio));
            title(title_str)
            hold off
            plot_idx = plot_idx + 1;
         end
      end   

      % Plots iterate cutoffs for various dual variable rel err tolerances
      y_range = 1e-2:-1e-6:1e-6;
      plot_idx = 0;
      for L_idx = 1:2
         for noise_ratio_idx = 1:3
            idx = 1;
            for tol = y_range
               val = find(y_diff(L_idx, noise_ratio_idx, 2:end) <= tol, 1);
               if ~isempty(val)
                  y_tol_iter(L_idx, noise_ratio_idx, idx) = val;
               else
                  y_tol_iter(L_idx, noise_ratio_idx, idx) = NaN;
               end
               idx = idx + 1;
            end
            subplot(6, 2, 2 + 2*plot_idx)
            semilogy(vec(y_tol_iter(L_idx, noise_ratio_idx, :)), y_range, 'b')
            hold on
            
            % plots tol vals
            for tol_temp = [1e-2, 1e-3, 1e-4, 1e-5, 1e-6]
               y_idx = find(vec(y_diff(L_idx, noise_ratio_idx, 2:end)) <= tol_temp, 1);
               if ~isempty(y_idx)
                  plot(y_idx, tol_temp, 'or');
               end
            end
            
            % plots window of termination
            iter_temp_min = model_iter_min(L_idx, noise_ratio_idx);
            iter_temp_max = model_iter_max(L_idx, noise_ratio_idx);
            x_points = [iter_temp_min, iter_temp_min, iter_temp_max, iter_temp_max];  
            y_points = [1e-6, 1e-2, 1e-2, 1e-6];
            color = [0, 0, 1];
            a = fill(x_points, y_points, color);
            a.FaceAlpha = 0.1;
            
            xlim([1, 1001])
            ylim([min(y_range), max(y_range)])
            set(gca,'YTick',flip(logspace(-2,-6,3)))
            L = L_range(L_idx);
            noise_ratio = noise_ratio_range(noise_ratio_idx);
            title_str = strcat("Dual variable: L = ", num2str(L), ", noise ratio = ", num2str(noise_ratio));
            title(title_str)
            hold off
            plot_idx = plot_idx + 1;
         end
      end




      % Prints table of iteration at which saga would terminate for prescribed
      % tolerances
      yTol = 1e-4;
      prTol = 1e-5;
      model_converged_iter = [];
      for L_idx = 1:2
         for noise_ratio_idx = 1:3
            for i = 2:1001
               val1 = y_diff(L_idx, noise_ratio_idx, i) <= yTol;
               val2 = pr_diff(L_idx, noise_ratio_idx, i) <= prTol;
               if val1 && val2
                  model_converged_iter(L_idx, noise_ratio_idx) = i;
                  break
               end
            end
         end
      end
      fprintf('\n')
      fprintf('Iterate at which model would terminate for primal tol %1.2e and dual iterate tol %1.2e.\n', prTol, yTol);
      fprintf('             |  L = %1i | L = %1i \n', L_range(1), L_range(2)');
      noise_ratio_idx = 1;
      fprintf('noise = %1.2f |   %3i  |   %3i  \n', noise_ratio_range(noise_ratio_idx), ...
         model_converged_iter(1, noise_ratio_idx), ...
         model_converged_iter(2, noise_ratio_idx));
      noise_ratio_idx = 2;
      fprintf('noise = %1.2f |   %3i  |   %3i  \n', noise_ratio_range(noise_ratio_idx), ...
         model_converged_iter(1, noise_ratio_idx), ...
         model_converged_iter(2, noise_ratio_idx));   
      noise_ratio_idx = 3;
      fprintf('noise = %1.2f |   %3i  |   %3i  \n', noise_ratio_range(noise_ratio_idx), ...
         model_converged_iter(1, noise_ratio_idx), ...
         model_converged_iter(2, noise_ratio_idx));
   end

   
   
   if figure_4_on
      
      % min and max iter vals for proper convergence window
      % for each model
      model_iter_min = [200, 200, 100; 50, 50, 50];
      model_iter_max = [400, 400, 200; 200, 100, 100];

      figure;
      
      % Plots iterate cutoffs for various dual objective rel err tolerances
      du_obj_range = 1e-2:-1e-6:1e-6;
      plot_idx = 0;
      for noise_ratio_idx = 1:3
         for L_idx = 1:2
            idx = 1;
            for tol = du_obj_range
               val = find(du_obj_diff(L_idx, noise_ratio_idx, 2:end) <= tol, 1);
               if ~isempty(val)
                  du_obj_tol_iter(L_idx, noise_ratio_idx, idx) = val;
               else
                  du_obj_tol_iter(L_idx, noise_ratio_idx, idx) = NaN;
               end
               idx = idx + 1;
            end

            subplot(3, 2, 1+plot_idx)
            semilogy(vec(du_obj_tol_iter(L_idx, noise_ratio_idx, :)), du_obj_range, 'b')
            hold on
            % plots tol vals
            for tol_temp = [1e-2, 1e-3, 1e-4, 1e-5, 1e-6]
               du_idx = find(vec(du_obj_diff(L_idx, noise_ratio_idx, 2:end)) <= tol_temp, 1);
               if ~isempty(du_idx)
                  plot(du_idx, tol_temp, 'or');
               end
            end
            
            % plots window of termination
            iter_temp_min = model_iter_min(L_idx, noise_ratio_idx);
            iter_temp_max = model_iter_max(L_idx, noise_ratio_idx);
            x_points = [iter_temp_min, iter_temp_min, iter_temp_max, iter_temp_max];  
            y_points = [1e-6, 1e-2, 1e-2, 1e-6];
            color = [0, 0, 1];
            a = fill(x_points, y_points, color);
            a.FaceAlpha = 0.1;
            
            xlim([1, 1001])
            ylim([min(du_obj_range), max(du_obj_range)])
            set(gca,'YTick',flip(logspace(-2,-6,3)))         
            L = L_range(L_idx);
            noise_ratio = noise_ratio_range(noise_ratio_idx);
            title_str = strcat("Dual objective: L = ", num2str(L), ", noise ratio = ", num2str(noise_ratio));
            title(title_str)
            hold off
            plot_idx = plot_idx + 1;
         end
      end
      
   end
   
   
   
   
   if figure_5_on
      % Plots signal relative error and du_gap to demonstrate that signal outliers
      % can be identified reliably using du_gap outliers
      iter_start = 50;
      iter_end = 200; 
   %iter_start = 145; iter_end = 149;
      iter_terminated = 148;
      L_idx = 2;
      noise_ratio_idx = 1;

      % Finds iterate with lowest du_gap value and identifies index
      du_gap_vec = vec(du_gap(L_idx, noise_ratio_idx, 2:end));
      [min_val, min_idx] = min(du_gap_vec(iter_terminated-19:iter_terminated));
      min_idx = min_idx - 1 + (iter_terminated-19);

      figure
      subplot(2,1,1)
      pr_vec = vec(sig_rel_err(L_idx, noise_ratio_idx, 2:end));
      plot(pr_vec)
      hold on
      outlier_idx = [];
      for iter = iter_start:iter_end
         subseq_start = max(iter-5, iter_start);
         subseq_end = min(iter+5, iter_end);
         scaleSignal = mean(pr_vec(subseq_start:subseq_end,1));
         if pr_vec(iter) >= 1.0005*scaleSignal
            outlier_idx = [outlier_idx; iter];
         end
      end
      xlim([iter_start, iter_end])
      ylim([.999*min(pr_vec(iter_start:iter_end)), 1.001*max(pr_vec(iter_start:iter_end))])

      plot(outlier_idx, pr_vec(outlier_idx), 'o') 
      plot(min_idx, pr_vec(min_idx), 'k*')
      title('Signal relative error')
      ylabel('Relative error')
      xlabel('Iteration')
      hold off

      %scaledu_gap = mean(vec(data.mGap(1,1,iter_start:iter_end)));
      %du_gapRescaled = (scaleSignal/scaledu_gap)*vec(data.du_gap(1,1,iter_start:iter_end));
      subplot(2,1,2)
      plot(du_gap_vec)
      hold on
      plot(outlier_idx, du_gap_vec(outlier_idx), 'o')
      %plot(iter_terminated, du_gap_vec(iter_terminated), 'x')

      plot(min_idx, min_val, 'k*')
      xlim([iter_start, iter_end])
      ylim([.95*min(du_gap_vec(iter_start:iter_end)), 1.05*max(du_gap_vec(iter_start:iter_end))])
   %   title_str = strcat('Dual variable: L = ', num2str(L_range(L_idx)), ...
   %      ', noise ratio =  ', num2str(noise_ratio_range(noise_ratio_idx)));
      title('Duality gap')
      ylabel('Error')
      xlabel('Iteration')

      hold off
   end

   
   
   if figure_6_on
   % Plots figure of iterates vs signal relative error across models
   % with red circle for each termination iterate
      figure
      iter_start = 2;
      iter_max = 501;
      i = 1;
      for noise_ratio_idx = 1:3
         for L_idx = 1:2
            L = L_range(L_idx);
            noise_ratio = noise_ratio_range(noise_ratio_idx);
            subplot(3, 2, i);
            seq_plot = vec(sig_rel_err(L_idx, noise_ratio_idx, iter_start:iter_max));
            plot(iter_start-1:iter_max-1, seq_plot, 'b')
            title_str = strcat("L = ", num2str(L), ", noise ratio = ", num2str(noise_ratio));
            title(title_str)
            xlabel('Iteration')
            ylabel('Relative error')
            %xlim([iter_start - 2, iter_max-1]);
            plot_buffer = (max(seq_plot) - min(seq_plot)) / 20;
            ylim([min(seq_plot)- plot_buffer, max(seq_plot) + plot_buffer])
            % Plots point at which iterate converged
            hold on
            iter_converged = model_converged_iter(L_idx, noise_ratio_idx);
            if ~isnan(iter_converged)
               plot(iter_converged, seq_plot(iter_converged), 'or');
            end
            legend('Signal rel err', strcat("Final iter = ", num2str(iter_converged)))
            hold off
%{
            subplot(6, 2, 2*i);
            iter_start_temp = 101;
            seq_plot_temp = vec(sig_rel_err(L_idx, noise_ratio_idx, iter_start_temp:iter_max));
            plot(iter_start_temp-1:iter_max-1, seq_plot_temp, 'b')
            title_str = strcat('latter iterates: L = ', num2str(L), ', noise ratio =  ', num2str(noise_ratio));
            title(title_str)
            xlabel('iteration')
            ylabel('error')
            %xlim([iter_start_temp - 2, iter_max-1])

            plot_buffer = (max(seq_plot_temp) - min(seq_plot_temp)) / 20;
            xlim([iter_start_temp-1, iter_max-1]);
            ylim([min(seq_plot_temp)- plot_buffer, max(seq_plot_temp) + plot_buffer])
            % Plots point at which iterate converged
            hold on
            iter_converged = model_converged_iter(L_idx, noise_ratio_idx);
            if ~isnan(iter_converged)
               plot(iter_converged, seq_plot(iter_converged), 'or');
            end
            hold off
%}
            i = i + 1;
         end
      end    
   end
   
end % print_figs_and_tables



data = struct('sig_rel_err', sig_rel_err, 'sig_diff', sig_diff, ...
   'du_gap', du_gap, 'du_gap_diff', du_gap_diff, 'du_obj', du_obj, ...
   'du_obj_diff', du_obj_diff, 'y_diff', y_diff, 'stepLS', stepLS, ...
   'pr_rel_err', pr_rel_err, 'pr_diff', pr_diff);

end

