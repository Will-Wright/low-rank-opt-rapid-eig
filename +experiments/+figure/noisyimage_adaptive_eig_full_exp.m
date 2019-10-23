function [data] = noisyimage_adaptive_eig_full_exp()
% function [data] = noisyimage_adaptive_eig_full_exp()
%
% Plots figures showing performance of `eigs`
% solving EMEP with various parameters (e.g., blocksize).


signal = 'image';
image_file = 'data/parrot_4k.jpg';
file_test = 'eigs.mat';
folder_exp = 'parrot_4k/';
solver = 'eigs';
n = 4096;
folder_root = 'cache/figure_noisyimage_adaptive_eig_full_exp/';


adaptive_shift_rate = 1; % default 1

use_lin_interp_in_ada = true; % default true
num_lin_interp_pts = 4; % default 4
lin_interp_weights = ones(num_lin_interp_pts,1);


print_num_mat_vecs = true;
print_num_mat_vecs_in_TeX_mode = false;

plot_2_figs_num_mvs_adaptive_various_blocksize = false;

print_table_opt_vs_ada_block_size_40 = true;
print_opt_vs_ada_in_TeX_mode = false;

plot_num_eigs = false;

plot_mat_vec_orig_vs_opt_surf = false;

plot_num_eigs_ada_vs_orig = true;

plot_mat_vec_and_eig_diff_surf = false;

plot_mat_vec_for_blocksize_vs_num_eigs_surf = false;



iter_start_min_mat_vec_2 = 1;
iter_end_min_mat_vec_2 = 50;

rel_err_fail_tol = 5e-11;



% Loads solver performance data for generating figures
if exist(strcat(folder_root, folder_exp, file_test)) ~= 2 % no test file found
   % Solves EMEPs with requested solver and varying parameters.
   % Saves solver info used in figures below.
   experiments.exputil.solve_EMEPs_various_solvers_and_params(...
      'signal', signal, 'n', n, 'solver', solver, ...
      'folder_root', folder_root, 'image_file', image_file);
   load(strcat(folder_root, folder_exp, file_test), 'results_EMEP_solvers')
else
   load(strcat(folder_root, folder_exp, file_test), 'results_EMEP_solvers')
   if ~results_EMEP_solvers.exp_is_complete
      experiments.exputil.solve_EMEPs_various_solvers_and_params(...
         'signal', signal, 'n', n, 'solver', solver, ...
         'folder_root', folder_root, 'image_file', image_file);
      load(strcat(folder_root, folder_exp, file_test), 'results_EMEP_solvers')
   end
end 
data.eigs.results = results_EMEP_solvers.eigs;
[data.eigs.L_range, data.eigs.noise_ratio_range, ...
   data.eigs.block_size_range, data.eigs.num_eigs_range] ...
   = size(data.eigs.results);



% Collects eigs results for num mat-vecs across EMEP tests
for L_idx = 1:data.eigs.L_range
   for noise_ratio_idx = 1:data.eigs.noise_ratio_range
      num_EMEP_iters = length(data.eigs.results{L_idx, noise_ratio_idx, 1, 1}.d);
      eigs_mv_array = zeros(num_EMEP_iters, data.eigs.block_size_range, data.eigs.num_eigs_range);
      eigs_failed_array = zeros(num_EMEP_iters, data.eigs.block_size_range, data.eigs.num_eigs_range);
      % Collects data for fixed L, noise_ratio, 
      % across eigs parameters and EMEP iterates
      % Also determines if solver failed
      for EMEP_iter = 1:num_EMEP_iters
         for block_size_idx = 1:data.eigs.block_size_range
            for num_eigs_idx = 1:data.eigs.num_eigs_range
               
             
               % Stops early on this loop for smallest block size due to
               % experimental structure (num eigs limited by small block size)
               if block_size_idx > 1 || num_eigs_idx < 15
                  eigs_mv_array(EMEP_iter, block_size_idx, num_eigs_idx) = data.eigs.results{...
                     L_idx, noise_ratio_idx, block_size_idx, num_eigs_idx}.num_mat_vec{1, EMEP_iter};

                  tol_temp = max(rel_err_fail_tol, (10^(-(EMEP_iter-4))) );
                  % Determines if this call to eigs failed
                  rel1 = data.eigs.results{L_idx, noise_ratio_idx, block_size_idx, num_eigs_idx}.rel_err1{1,EMEP_iter};
                  rel2 = data.eigs.results{L_idx, noise_ratio_idx, block_size_idx, num_eigs_idx}.rel_err2{1,EMEP_iter};
                  if ~isempty(rel1) && ~isempty(rel2) && ~isnan(rel1) && ~isnan(rel2) && (rel1 <= tol_temp) && (rel2 <= tol_temp)
                     eigs_failed_array(EMEP_iter, block_size_idx, num_eigs_idx) = 0;
                  else
                     eigs_failed_array(EMEP_iter, block_size_idx, num_eigs_idx) = 1;
                  end
               end
               
            end
         end
      end
      data.eigs.mv_array{L_idx, noise_ratio_idx} = eigs_mv_array;
      data.eigs.failed_array{L_idx, noise_ratio_idx} = eigs_failed_array;
   end
end





%%%%%%%%%%%%%
%           %
%  ANALYZE  %
%           %
%%%%%%%%%%%%%



% Determines optimal eigs parameter for each EMEP
for block_size_idx = 1:data.eigs.block_size_range
   for L_idx = 1:data.eigs.L_range
      for noise_ratio_idx = 1:data.eigs.noise_ratio_range
         for EMEP_iter = 1:length(data.eigs.results{L_idx, noise_ratio_idx, 1, 1}.d)
            if block_size_idx == 1
               % smallest block size has smaller max num eigs
               [opt_num_mvs, opt_num_eigs_idx] ...
                  = min(data.eigs.mv_array{L_idx, noise_ratio_idx}(EMEP_iter,block_size_idx,1:14));
            else
               [opt_num_mvs, opt_num_eigs_idx] ...
                  = min(data.eigs.mv_array{L_idx, noise_ratio_idx}(EMEP_iter,block_size_idx,:));
            end
            data.eigs.opt_params{L_idx, noise_ratio_idx, block_size_idx}.num_mvs(EMEP_iter) = opt_num_mvs;
            data.eigs.opt_params{L_idx, noise_ratio_idx, block_size_idx}.num_eigs(EMEP_iter) = 1+opt_num_eigs_idx;
            data.eigs.opt_params{L_idx, noise_ratio_idx, block_size_idx}.fails(EMEP_iter) ...
               = data.eigs.failed_array{L_idx, noise_ratio_idx}(EMEP_iter, block_size_idx, opt_num_eigs_idx);
         end   
      end
   end
end



% Generates results for original eigs parameters
block_size_idx = 1;
num_eigs_idx = 1;
for L_idx = 1:data.eigs.L_range
   for noise_ratio_idx = 1:data.eigs.noise_ratio_range      
      data.eigs.orig_params{L_idx, noise_ratio_idx}.num_mvs ...
         = data.eigs.mv_array{L_idx, noise_ratio_idx}(:,block_size_idx,num_eigs_idx);
      data.eigs.orig_params{L_idx, noise_ratio_idx}.fails ...
         = data.eigs.failed_array{L_idx, noise_ratio_idx}(:, block_size_idx, num_eigs_idx);
   end
end



% Simulates adaptive strategy for eigs parameter
for block_size_idx = 1:data.eigs.block_size_range
   if block_size_idx > 1
      num_eigs_range = data.eigs.num_eigs_range;
   else
      num_eigs_range = 14;
   end
   for L_idx = 1:data.eigs.L_range
      for noise_ratio_idx = 1:data.eigs.noise_ratio_range
         % Initializes num_eigs = 2
         
         opts.num_eigs_min = 2;
         opts.num_eigs_max = num_eigs_range+1;
         opts.num_lin_interp_pts = num_lin_interp_pts;
         
         shift_rate = adaptive_shift_rate;
         EMEP_iter = 1;
         num_eigs_idx(EMEP_iter) = 1;
         num_eigs(EMEP_iter) = num_eigs_idx(EMEP_iter) + 1;
         ada_num_mvs(EMEP_iter) = data.eigs.mv_array{L_idx, noise_ratio_idx}(EMEP_iter,block_size_idx,num_eigs_idx(EMEP_iter));
         data.eigs.ada_params{L_idx, noise_ratio_idx, block_size_idx}.num_mvs(EMEP_iter) = ada_num_mvs(EMEP_iter);
         data.eigs.ada_params{L_idx, noise_ratio_idx, block_size_idx}.num_eigs(EMEP_iter) = 1+num_eigs_idx(EMEP_iter);
         data.eigs.ada_params{L_idx, noise_ratio_idx, block_size_idx}.fails(EMEP_iter) ...
            = data.eigs.failed_array{L_idx, noise_ratio_idx}(EMEP_iter, block_size_idx, num_eigs_idx(EMEP_iter));

         for EMEP_iter = 2:length(data.eigs.results{L_idx, noise_ratio_idx, 1, 1}.d)         
            % Determine new num_eigs_idx
            
            [num_eigs_update] = adaptive_eigs_params(ada_num_mvs, num_eigs, EMEP_iter-1, opts);
            num_eigs(EMEP_iter) = num_eigs_update;
            num_eigs_idx(EMEP_iter) = num_eigs_update - 1;
            

%{
            % PROTOTYPING CODE
            if EMEP_iter == 2 || (num_eigs_idx(EMEP_iter-1) > num_eigs_idx(EMEP_iter-2))
               shift_dir = 1;
            else
               shift_dir = -1;
            end

            % Computes linear interpolation of previous results
            if (EMEP_iter > num_lin_interp_pts)
               row_idx = EMEP_iter-num_lin_interp_pts:EMEP_iter-1; 
               
               Y = ada_num_mvs(1,row_idx)';               
               X = [ones(num_lin_interp_pts,1), num_eigs_idx(1,row_idx)'];               
               B = lscov(X,Y,lin_interp_weights); 
               data.proto{L_idx, noise_ratio_idx}.B(EMEP_iter, :) = B;
               data.proto{L_idx, noise_ratio_idx}.Bn(EMEP_iter, :) = B/mean(ada_num_mvs(1,row_idx));               
               data.proto{L_idx, noise_ratio_idx}.std(EMEP_iter, 1) = std(ada_num_mvs(1,row_idx));
               data.proto{L_idx, noise_ratio_idx}.stdn(EMEP_iter, 1) = std(ada_num_mvs(1,row_idx)/mean(ada_num_mvs(1,row_idx)));               
               
               A = [ones(num_lin_interp_pts, 1), num_eigs_idx(1,row_idx)'];
               b = ada_num_mvs(1,row_idx)';
               x = A\b;
               delta_4 = -sign(x(2,1));
               
            end
                        
            % Forces shift away from min or max inner_iter value         
            if num_eigs_idx(EMEP_iter-1) == 1
               num_eigs_idx(EMEP_iter) = num_eigs_idx(EMEP_iter-1) + shift_rate;
            elseif num_eigs_idx(EMEP_iter-1) == num_eigs_range
               num_eigs_idx(EMEP_iter) = num_eigs_idx(EMEP_iter-1) - shift_rate;
               
            % Shifts based on previous results  
            elseif use_lin_interp_in_ada && (EMEP_iter > num_lin_interp_pts) ...
                  && ( delta_4 == sign((ada_num_mvs(EMEP_iter-2) / ada_num_mvs(EMEP_iter-1) - 1)*shift_dir) )
               %fprintf('L = %i, nr = %i, EMEP iter = %i num eigs = %i \n', L_idx, noise_ratio_idx, EMEP_iter, num_eigs_idx(EMEP_iter-1));

               shift_dir = 2*delta_4;
               num_eigs_idx(EMEP_iter) = num_eigs_idx(EMEP_iter-1) + shift_dir;
               
               num_eigs_idx(EMEP_iter) = min(max(num_eigs_idx(EMEP_iter), 1), num_eigs_range);   

            else
               
               if ada_num_mvs(EMEP_iter-1) <= ada_num_mvs(EMEP_iter-2)
                  num_eigs_idx(EMEP_iter) = num_eigs_idx(EMEP_iter-1) + shift_rate*shift_dir;
               else
                  num_eigs_idx(EMEP_iter) = num_eigs_idx(EMEP_iter-1) - shift_rate*shift_dir;
               end
            end
%}            
            
            
            
            % Determine new ada_num_mvs
            ada_num_mvs(EMEP_iter) = data.eigs.mv_array{L_idx, noise_ratio_idx}(EMEP_iter, block_size_idx, num_eigs_idx(EMEP_iter));
            data.eigs.ada_params{L_idx, noise_ratio_idx, block_size_idx}.num_mvs(EMEP_iter) = ada_num_mvs(EMEP_iter);
            data.eigs.ada_params{L_idx, noise_ratio_idx, block_size_idx}.num_eigs(EMEP_iter) = 1+num_eigs_idx(EMEP_iter);    
            data.eigs.ada_params{L_idx, noise_ratio_idx, block_size_idx}.fails(EMEP_iter) ...
               = data.eigs.failed_array{L_idx, noise_ratio_idx}(EMEP_iter, block_size_idx, num_eigs_idx(EMEP_iter));
         end   
      end
   end
end






%%%%%%%%%%%%%
%           %
%  FIGURES  %
%           %
%%%%%%%%%%%%%


% Prints maximum number of mat-vecs and number of fails
if print_num_mat_vecs
   fprintf('\nTable of num mat-vecs for original EMEP setting vs adaptive inner iteration size of IRAM\n')
   fprintf('                 |  original  |       adaptive inner iter size (m), (second term is num fails) \n')
   fprintf(' L    nr    nIts |  j=2, m=20 |   m = 20   |   m = 40   |   m = 60   |   m = 80   |   m = 100 \n')
   for L_idx = 1:data.eigs.L_range      
      for noise_ratio_idx = 1:data.eigs.noise_ratio_range
         L = data.eigs.results{L_idx, noise_ratio_idx}.exp_params.L;
         noise_ratio = data.eigs.results{L_idx, noise_ratio_idx}.exp_params.noise_ratio;
         num_EMEP_iters = length(data.eigs.results{L_idx, noise_ratio_idx,1,1}.d);
         for block_size_idx = 1:data.eigs.block_size_range
            mvs_ada(block_size_idx) = sum(data.eigs.ada_params{L_idx, noise_ratio_idx, block_size_idx}.num_mvs);
            fails_ada(block_size_idx) = sum(data.eigs.ada_params{L_idx, noise_ratio_idx, block_size_idx}.fails);
         end
         
         mvs_orig = sum(data.eigs.orig_params{L_idx, noise_ratio_idx}.num_mvs);
         fails_orig = sum(data.eigs.orig_params{L_idx, noise_ratio_idx}.fails);
         
         if print_num_mat_vecs_in_TeX_mode
            fprintf(' %2i &  %3.2f & %3i & %7i  ', L, noise_ratio, num_EMEP_iters, mvs_orig);
            for block_size_idx = 1:data.eigs.block_size_range
               fprintf('& %7i  ', mvs_ada(block_size_idx));
            end
            fprintf('\\\\ \n')
         else
            fprintf(' %2i   %3.2f   %3i | %7i %2i ', L, noise_ratio, num_EMEP_iters, mvs_orig, fails_orig);
            for block_size_idx = 1:data.eigs.block_size_range
               fprintf('| %7i %2i ', mvs_ada(block_size_idx), fails_ada(block_size_idx));
            end
            fprintf('\n')
         end
      end
   end
end








% Plots number of mat-vecs across EMEP iters for two experiments
if plot_2_figs_num_mvs_adaptive_various_blocksize
   figure;
   plot_num = 1;
   for L_idx = 1
      for noise_ratio_idx = 2:data.eigs.noise_ratio_range
         subplot(2,1,plot_num)
         hold on;
         L = data.eigs.results{L_idx, noise_ratio_idx}.exp_params.L;
         noise_ratio = data.eigs.results{L_idx, noise_ratio_idx}.exp_params.noise_ratio;
         num_EMEP_iters = length(data.eigs.results{L_idx, noise_ratio_idx,1,1}.d);
         
         %mv_orig = data.eigs.orig_params{L_idx, noise_ratio_idx}.num_mvs;
         mv_ada20 = data.eigs.ada_params{L_idx, noise_ratio_idx, 1}.num_mvs(1:end);
         mv_ada40 = data.eigs.ada_params{L_idx, noise_ratio_idx, 2}.num_mvs(1:end);
         mv_ada80 = data.eigs.ada_params{L_idx, noise_ratio_idx, 4}.num_mvs(1:end);
         
         %plot(1:num_EMEP_iters, mv_orig, 'x');
         plot(1:num_EMEP_iters, mv_ada20, '*');
         plot(1:num_EMEP_iters, mv_ada40, 'o');
         plot(1:num_EMEP_iters, mv_ada80, '+');
         legend('m = 20', 'm = 40', 'm = 80')
         xlim([1, num_EMEP_iters])
         
         title_str = strcat("L = ", num2str(L), ", noise ratio = ", num2str(noise_ratio));
         title(title_str)
         ylabel('Number mat-vecs')
         xlabel('EMEP iteration')
         
         
         hold off;
         plot_num = plot_num + 1;
      end
   end   
end



% Prints table comparing total mat-vecs for optimal inner iter vs adaptive
block_size_idx = 2;   % fixes block size = 40
if print_table_opt_vs_ada_block_size_40
   fprintf('\nTable of num mat-vecs for optimal EMEP setting vs adaptive inner iteration size of IRAM\n')
   fprintf('                 |  original |   optimal   | adaptive inner iter size (m), second value is %% decr from original \n')
   fprintf(' L    nr    nIts | m=20, j=2 |    m = 40   |  m = 40   \n')
   for L_idx = 1:data.eigs.L_range      
      for noise_ratio_idx = 1:data.eigs.noise_ratio_range
         L = data.eigs.results{L_idx, noise_ratio_idx}.exp_params.L;
         noise_ratio = data.eigs.results{L_idx, noise_ratio_idx}.exp_params.noise_ratio;
         num_EMEP_iters = length(data.eigs.results{L_idx, noise_ratio_idx,1,1}.d);
         mvs_orig = sum(data.eigs.orig_params{L_idx, noise_ratio_idx}.num_mvs);
         mvs_ada = sum(data.eigs.ada_params{L_idx, noise_ratio_idx, block_size_idx}.num_mvs); 
         mvs_ada_perc_decr = 1 - mvs_ada/mvs_orig;
         mvs_opt = sum(data.eigs.opt_params{L_idx, noise_ratio_idx, block_size_idx}.num_mvs);
         mvs_opt_perc_decr = 1 - mvs_opt/mvs_orig;
         if print_opt_vs_ada_in_TeX_mode
            fprintf(' %2i &  %3.2f & %3i & %7i  ', L, noise_ratio, num_EMEP_iters, mvs_orig);
            fprintf('& %7i & %2.0f\\%% ', mvs_opt, 100*mvs_opt_perc_decr);
            fprintf('& %7i & %2.0f\\%% ', mvs_ada, 100*mvs_ada_perc_decr);
            fprintf('\\\\ \n')
         else
            fprintf(' %2i   %3.2f   %3i |  %7i  ', L, noise_ratio, num_EMEP_iters, mvs_orig);
            fprintf('|%7i %3.2f ', mvs_opt, mvs_opt_perc_decr);
            fprintf('|%7i %3.2f ', mvs_ada, mvs_ada_perc_decr);
            fprintf('\n')
         end
      end
   end
end



if plot_num_eigs_ada_vs_orig
   
   block_size_idx = 2; % fixes block size at 40
   iter_start_array = [1, 1, 1; 1, 1, 1];
   iter_end_array = [length(data.eigs.results{1, 1,1,1}.d),...
      length(data.eigs.results{1, 2,1,1}.d),...
      length(data.eigs.results{1, 3,1,1}.d);...
      length(data.eigs.results{2, 1,1,1}.d),...
      length(data.eigs.results{2, 2,1,1}.d),...
      length(data.eigs.results{2, 3,1,1}.d)];
   num_eigs_start_array = [2, 2, 2; 2, 2, 2];
   num_eigs_end_array = [25, 25, 25; 25, 25, 25];

   
   figure;
   L_idx = 1;
   noise_ratio_idx = 2;
   hold on

   L = data.eigs.results{L_idx, noise_ratio_idx}.exp_params.L;
   noise_ratio = data.eigs.results{L_idx, noise_ratio_idx}.exp_params.noise_ratio;
   iter_start = iter_start_array(L_idx, noise_ratio_idx);
   iter_end = iter_end_array(L_idx, noise_ratio_idx);

   plot(data.eigs.ada_params{L_idx, noise_ratio_idx,block_size_idx}.num_eigs(iter_start:iter_end), 'r*')
   plot(data.eigs.opt_params{L_idx, noise_ratio_idx,block_size_idx}.num_eigs(iter_start:iter_end), 'bo')
   legend('IGDD algorithm', 'Optimal params')

   title_str = strcat("L = ", num2str(L), ", noise ratio = ", num2str(noise_ratio));
   title(title_str)
   xlim([1 , iter_end])
   ylabel('Number eigenvalues requested')
   xlabel('EMEP iteration')
      
   hold off
   
   
   L_idx = 2;
   noise_ratio_idx = 1;
   figure;
   hold on

   L = data.eigs.results{L_idx, noise_ratio_idx}.exp_params.L;
   noise_ratio = data.eigs.results{L_idx, noise_ratio_idx}.exp_params.noise_ratio;
   iter_start = iter_start_array(L_idx, noise_ratio_idx);
   iter_end = iter_end_array(L_idx, noise_ratio_idx);

   plot(data.eigs.ada_params{L_idx, noise_ratio_idx,block_size_idx}.num_eigs(iter_start:iter_end), 'r*')
   plot(data.eigs.opt_params{L_idx, noise_ratio_idx,block_size_idx}.num_eigs(iter_start:iter_end), 'bo')
   legend('IGDD algorithm', 'Optimal params')

   title_str = strcat("L = ", num2str(L), ", noise ratio = ", num2str(noise_ratio));
   title(title_str)
   xlim([1 , iter_end])
   ylabel('Number eigenvalues requested')
   xlabel('EMEP iteration')
      
   hold off

end




if plot_num_eigs
   figure;
   
   block_size_idx = 2; % fixes block size at 40
   iter_start_array = [1, 1, 1; 1, 1, 1];
   iter_end_array = [length(data.eigs.results{1, 1,1,1}.d),...
      length(data.eigs.results{1, 2,1,1}.d),...
      length(data.eigs.results{1, 3,1,1}.d);...
      length(data.eigs.results{2, 1,1,1}.d),...
      length(data.eigs.results{2, 2,1,1}.d),...
      length(data.eigs.results{2, 3,1,1}.d)];
   num_eigs_start_array = [2, 2, 2; 2, 2, 2];
   num_eigs_end_array = [25, 25, 25; 25, 25, 25];
   L_range = 1;
   noise_ratio_range = 2:3;
   
   num_rows = 1;
   num_cols = length(L_range)*length(noise_ratio_range);
   plot_num = 1;
   
   for L_idx = L_range
      for noise_ratio_idx = noise_ratio_range
                 
         subplot(num_rows, num_cols, plot_num)
         hold on
         
         
         L = data.eigs.results{L_idx, noise_ratio_idx}.exp_params.L;
         noise_ratio = data.eigs.results{L_idx, noise_ratio_idx}.exp_params.noise_ratio;
         iter_start = iter_start_array(L_idx, noise_ratio_idx);
         iter_end = iter_end_array(L_idx, noise_ratio_idx);
         num_eigs_start = num_eigs_start_array(L_idx, noise_ratio_idx);
         num_eigs_end = num_eigs_end_array(L_idx, noise_ratio_idx);
      
         
         plot(data.eigs.ada_params{L_idx, noise_ratio_idx,block_size_idx}.num_eigs(iter_start:iter_end), 'r*')
         
         title_str = strcat("L = ", num2str(L), ", noise ratio = ", num2str(noise_ratio));
         title(title_str)
         xlim([1 , iter_end])
         ylabel('Number eigenvalues requested')
         xlabel('EMEP iteration')

         
         hold off
         plot_num = plot_num + 1;
      end
   end
   
   
end






if plot_mat_vec_orig_vs_opt_surf
     
   block_size_idx = 2; % fixes block size at 40
   iter_start_array = [1, 50, 1; 1, 1, 1];
   iter_end_array = [length(data.eigs.results{1, 1,1,1}.d),...
      200,... %length(data.eigs.results{1, 2,1,1}.d),...
      length(data.eigs.results{1, 3,1,1}.d);...
      length(data.eigs.results{2, 1,1,1}.d),...
      length(data.eigs.results{2, 2,1,1}.d),...
      length(data.eigs.results{2, 3,1,1}.d)];
   num_eigs_start_array = [2, 2, 2; 2, 2, 2];
   num_eigs_end_array = [25, 25, 25; 25, 25, 25];
%   num_eigs_end_array = [15, 15, 15; 15, 15, 15];
   
   mv_shift_for_dots = [20, 20, 20; 20, 20, 20];
      
   normalize_mat_vecs = false;
   trim_max_vals = true;
   trim_mv_val = 1500;
   print_dots_on_adaptive_vals = true;
   
   
   L_idx = 1;
   noise_ratio_idx = 2;
   figure;

   hold on

   L = data.eigs.results{L_idx, noise_ratio_idx}.exp_params.L;
   noise_ratio = data.eigs.results{L_idx, noise_ratio_idx}.exp_params.noise_ratio;
   iter_start = iter_start_array(L_idx, noise_ratio_idx);
   iter_end = iter_end_array(L_idx, noise_ratio_idx);
   num_eigs_start = num_eigs_start_array(L_idx, noise_ratio_idx);
   num_eigs_end = num_eigs_end_array(L_idx, noise_ratio_idx);

   [X, Y] = meshgrid(iter_start:iter_end, num_eigs_start:num_eigs_end);

   mat_vec_temp = data.eigs.mv_array{L_idx, noise_ratio_idx}(iter_start:iter_end, block_size_idx, num_eigs_start-1:num_eigs_end-1);
   mat_vec_temp = permute(mat_vec_temp, [3,1,2]);

   if normalize_mat_vecs
      for i = 1 : iter_end - iter_start + 1
         mat_vec_temp(:,i) = mat_vec_temp(:,i) / max(mat_vec_temp(:,i));
      end
   elseif trim_max_vals
      % Eliminates outliers in plot
      mat_vec_temp = min(mat_vec_temp, trim_mv_val*ones(size(mat_vec_temp)));
   end

   Z_mv = reshape(mat_vec_temp, size(X));
   surf(X, Y, Z_mv)
   xlim([iter_start_array(L_idx, noise_ratio_idx), iter_end_array(L_idx, noise_ratio_idx)]);
   title_str = strcat("L = ", num2str(L), ", noise ratio = ", num2str(noise_ratio));
   title(title_str)
   zlabel('Number mat-vecs')
   ylabel('Number eigenvalues requested')
   xlabel('EMEP iteration')



   % Plots dots on inner iter val for optimal IRAM (eigs) parameters
   if print_dots_on_adaptive_vals      
      for EMEP_iter = iter_start:iter_end
         eig_num_idx = data.eigs.opt_params{L_idx, noise_ratio_idx,block_size_idx}.num_eigs(EMEP_iter);
         mv_val = data.eigs.opt_params{L_idx, noise_ratio_idx,block_size_idx}.num_mvs(EMEP_iter);
         if trim_max_vals
            mv_val = min(mv_val, trim_mv_val);
         end
         scatter3(EMEP_iter, eig_num_idx, mv_shift_for_dots(L_idx, noise_ratio_idx) + mv_val, 50, 'filled', 'k')
      end
   end

   hold off         

   figure;
   hold on

   L = data.eigs.results{L_idx, noise_ratio_idx}.exp_params.L;
   noise_ratio = data.eigs.results{L_idx, noise_ratio_idx}.exp_params.noise_ratio;
   num_EMEP_iters = length(data.eigs.results{L_idx, noise_ratio_idx,1,1}.d);

   mv_opt = data.eigs.opt_params{L_idx, noise_ratio_idx, block_size_idx}.num_mvs(1:end);
   mv_orig = data.eigs.orig_params{L_idx, noise_ratio_idx}.num_mvs(1:end);

   plot(1:num_EMEP_iters, mv_orig, '*');
   plot(1:num_EMEP_iters, mv_opt, '+');
   legend('Original params (GDD algorithm)', 'Optimal params')
   xlim([1, num_EMEP_iters])

   title_str = strcat("L = ", num2str(L), ", noise ratio = ", num2str(noise_ratio));
   title(title_str)
   ylabel('Number mat-vecs')
   xlabel('EMEP iteration')
   hold off
   
end








% Plot surfaces of mat-vec counts
if plot_mat_vec_and_eig_diff_surf
   
   block_size_idx = 2; % fixes block size at 40
   iter_start_array = [10, 50, 40; 1, 1, 1];
   iter_end_array = [length(data.eigs.results{1, 1,1,1}.d),...
      200,...%length(data.eigs.results{1, 2,1,1}.d),...
      length(data.eigs.results{1, 3,1,1}.d);...
      length(data.eigs.results{2, 1,1,1}.d),...
      length(data.eigs.results{2, 2,1,1}.d),...
      length(data.eigs.results{2, 3,1,1}.d)];
   num_eigs_start_array = [2, 2, 2; 2, 2, 2];
   num_eigs_end_array = [25, 25, 25; 25, 25, 25];
   
   mv_shift_for_dots = [20, 20, 20; 20, 20, 20];
   eig_diff_shift_for_dots = [2e-5, 1e-5, 1e-5; 2e-5, 2e-5, 2e-5];
   
   L_range = 1;
   noise_ratio_range = 2:3;
   
   
   normalize_mat_vecs = false;
   trim_max_vals = true;
   trim_mv_val = 2000;
   trim_eig_diff_val = 1e-3;
   print_dots_on_adaptive_vals = true;

   num_rows = length(L_range)*length(noise_ratio_range);
   num_cols = 2;
   plot_num = 1;
   
   for L_idx = L_range
      for noise_ratio_idx = noise_ratio_range
         figure;
         plot_num = 1;
         
         subplot(2, 1, plot_num)
         hold on
         
         L = data.eigs.results{L_idx, noise_ratio_idx}.exp_params.L;
         noise_ratio = data.eigs.results{L_idx, noise_ratio_idx}.exp_params.noise_ratio;
         iter_start = iter_start_array(L_idx, noise_ratio_idx);
         iter_end = iter_end_array(L_idx, noise_ratio_idx);
         num_eigs_start = num_eigs_start_array(L_idx, noise_ratio_idx);
         num_eigs_end = num_eigs_end_array(L_idx, noise_ratio_idx);

         [X, Y] = meshgrid(iter_start:iter_end, num_eigs_start:num_eigs_end);
         
         mat_vec_temp = data.eigs.mv_array{L_idx, noise_ratio_idx}(iter_start:iter_end, block_size_idx, num_eigs_start-1:num_eigs_end-1);
         mat_vec_temp = permute(mat_vec_temp, [3,1,2]);
         
         if normalize_mat_vecs
            for i = 1 : iter_end - iter_start + 1
               mat_vec_temp(:,i) = mat_vec_temp(:,i) / max(mat_vec_temp(:,i));
            end
         elseif trim_max_vals
            % Eliminates outliers in plot
            mat_vec_temp = min(mat_vec_temp, trim_mv_val*ones(size(mat_vec_temp)));
         end

         Z_mv = reshape(mat_vec_temp, size(X));
         surf(X, Y, Z_mv)
         xlim([iter_start_array(L_idx, noise_ratio_idx), iter_end_array(L_idx, noise_ratio_idx)]);
         title_str = strcat("L = ", num2str(L), ", noise ratio = ", num2str(noise_ratio));
         title(title_str)
         zlabel('Number mat-vecs')
         ylabel('Number eigenvalues requested')
         xlabel('EMEP iteration')
         
         

         % Plots dots on inner iter val for adaptive method
         if print_dots_on_adaptive_vals      
            for EMEP_iter = iter_start:iter_end
               eig_num_idx = data.eigs.ada_params{L_idx, noise_ratio_idx,block_size_idx}.num_eigs(EMEP_iter);
               mv_val = data.eigs.ada_params{L_idx, noise_ratio_idx,block_size_idx}.num_mvs(EMEP_iter);
               if trim_max_vals
                  mv_val = min(mv_val, trim_mv_val);
               end
               scatter3(EMEP_iter, eig_num_idx, mv_shift_for_dots(L_idx, noise_ratio_idx) + mv_val, 50, 'filled', 'k')
            end
         end
   
         hold off
         plot_num = plot_num + 1;
         
         
         
         % Plot surfaces of differences between eigenvalues
         
         subplot(2, 1, plot_num)
         hold on
         
         [X, Y] = meshgrid(iter_start:iter_end, num_eigs_start:num_eigs_end);
         
         %mat_vec_temp = data.eigs.mv_array{L_idx, noise_ratio_idx}(iter_start:iter_end, block_size_idx, num_eigs_start-1:num_eigs_end-1);
         %mat_vec_temp = permute(mat_vec_temp, [3,1,2]);
         %Z_mv = reshape(mat_vec_temp, size(X));
         
         eig_diffs = nan(1+num_eigs_end-num_eigs_start, 1+iter_end-iter_start);
         
         i = 1;
         for EMEP_iter = iter_start:iter_end
            d = data.eigs.results{L_idx, noise_ratio_idx,5,29}.d{1,EMEP_iter};
            eig_diffs(:, i) = d(num_eigs_start:num_eigs_end) - d(num_eigs_start+1:num_eigs_end+1);
            i = i + 1;
         end

         if trim_max_vals
            eig_diffs = min(eig_diffs, trim_eig_diff_val*ones(size(eig_diffs)));            
         end
                  
         surf(X, Y, eig_diffs)
         xlim([iter_start_array(L_idx, noise_ratio_idx), iter_end_array(L_idx, noise_ratio_idx)]);
         title_str = strcat("L = ", num2str(L), ", noise ratio = ", num2str(noise_ratio));
         title(title_str)
         ylabel('Number eigenvalues requested')
         zlabel('Eigenvalue difference')
         xlabel('EMEP iteration')

         % Plots dots on inner iter val for adaptive method
         if print_dots_on_adaptive_vals      
            i = 1;
            for EMEP_iter = iter_start:iter_end
               eig_num_val = data.eigs.ada_params{L_idx, noise_ratio_idx,block_size_idx}.num_eigs(EMEP_iter);
               eig_diff_val = eig_diffs(eig_num_val-num_eigs_start+1, i);
               if trim_max_vals
                  eig_diff_val = min(eig_diff_val, trim_eig_diff_val);
               end               
               scatter3(EMEP_iter, eig_num_val, eig_diff_shift_for_dots(L_idx, noise_ratio_idx) + eig_diff_val, 50, 'filled', 'k')
               i = i + 1;
            end
         end
         
   
         hold off
         plot_num = plot_num + 1;

         
      end
   end
end






% xxx here

if plot_mat_vec_for_blocksize_vs_num_eigs_surf
       
   block_size_idx = 2; % fixes block size at 40
   iter_start_array = [1, 1, 1; 1, 1, 1];
   iter_end_array = [length(data.eigs.results{1, 1,1,1}.d),...
      length(data.eigs.results{1, 2,1,1}.d),...
      length(data.eigs.results{1, 3,1,1}.d);...
      length(data.eigs.results{2, 1,1,1}.d),...
      length(data.eigs.results{2, 2,1,1}.d),...
      length(data.eigs.results{2, 3,1,1}.d)];
   num_eigs_start_array = [2, 2, 2; 2, 2, 2];
   num_eigs_end_array = [25, 25, 25; 25, 25, 25];
%   num_eigs_end_array = [15, 15, 15; 15, 15, 15];
   
   mv_shift_for_dots = [20, 20, 20; 20, 20, 20];
      
   trim_max_vals = true;
   trim_mv_val = 3000;
   print_dots_on_adaptive_vals = true;
   
   m_range = 20:20:100;
   j_range = 2:1:15;
   
   
   figure;
   hold on
   
   iter = 180;
   L_idx = 1;
   noise_ratio_idx = 2;

   L = data.eigs.results{L_idx, noise_ratio_idx}.exp_params.L;
   noise_ratio = data.eigs.results{L_idx, noise_ratio_idx}.exp_params.noise_ratio;

   [X, Y] = meshgrid(m_range, j_range);

   mat_vec_temp = data.eigs.mv_array{L_idx, noise_ratio_idx}(iter, 1:5, 1:14);
   mat_vec_temp = permute(mat_vec_temp, [3,2,1]);

   if trim_max_vals
      % Eliminates outliers in plot
      mat_vec_temp = min(mat_vec_temp, trim_mv_val*ones(size(mat_vec_temp)));
   end
   
   Z_mv = reshape(mat_vec_temp, size(X));

   surf(X, Y, Z_mv)
   ylim([j_range(1), j_range(end)]);
   xticks([20, 40, 60, 80, 100]);
   title_str = strcat("L = ", num2str(L), ...
      ", noise ratio = ", num2str(noise_ratio), ...
      ", EMEP iter = ", num2str(iter));
   title(title_str)
   zlabel('Number mat-vecs')
   ylabel('Number eigenvalues requested')
   xlabel('Arnoldi decomp sz')

   hold off 

   
   figure;
   hold on
   
   iter = 85;
   L_idx = 1;
   noise_ratio_idx = 3;

   L = data.eigs.results{L_idx, noise_ratio_idx}.exp_params.L;
   noise_ratio = data.eigs.results{L_idx, noise_ratio_idx}.exp_params.noise_ratio;

   [X, Y] = meshgrid(m_range, j_range);

   mat_vec_temp = data.eigs.mv_array{L_idx, noise_ratio_idx}(iter, 1:5, 1:14);
   mat_vec_temp = permute(mat_vec_temp, [3,2,1]);

   if trim_max_vals
      % Eliminates outliers in plot
      mat_vec_temp = min(mat_vec_temp, trim_mv_val*ones(size(mat_vec_temp)));
   end
   
   Z_mv = reshape(mat_vec_temp, size(X));

   surf(X, Y, Z_mv)
   ylim([j_range(1), j_range(end)]);
   xticks([20, 40, 60, 80, 100]);
   title_str = strcat("L = ", num2str(L), ...
      ", noise ratio = ", num2str(noise_ratio), ...
      ", EMEP iter = ", num2str(iter));
   title(title_str)
   zlabel('Number mat-vecs')
   ylabel('Number eigenvalues requested')
   xlabel('Arnoldi decomp sz')

   hold off   
   
end












%{
L_idx = 2;
noise_ratio_idx = 3;


% Collects and organizes specific solver data for figures and tables
iter_range = length(data.eigs.results{L_idx, noise_ratio_idx, 1, 1}.num_mat_vec);

% Note: 3rd dim of mat_vec_array begins at 1, but num_eigs begins at 2
mat_vec_array = nan(iter_range, data.eigs.block_size_range, data.eigs.num_eigs_range);

% Gets number of mat-vec products for each experiment
for iter = 1:iter_range
   for block_size_idx = 1:data.eigs.block_size_range
      for num_eigs_idx = 1:data.eigs.num_eigs_range
         if ~isempty(data.eigs.results{L_idx, noise_ratio_idx, block_size_idx, num_eigs_idx})
            mat_vec_array(iter, block_size_idx, num_eigs_idx) ...
               = data.eigs.results{L_idx, noise_ratio_idx, block_size_idx, num_eigs_idx}.num_mat_vec{iter};
         end

      end      
   end
end

% Finds minimum number of mat-vecs with fixed blocksize for all iterates
for iter = 1:iter_range
   mat_vec_temp = reshape(mat_vec_array(iter, :, :), [data.eigs.block_size_range, data.eigs.num_eigs_range]);
   [~, min_num_eigs_idx] = min(mat_vec_temp, [], 2, 'omitnan');
   for block_size_idx = 1:data.eigs.block_size_range
      info.block_size(iter, block_size_idx) = data.eigs.results{L_idx, noise_ratio_idx, block_size_idx, min_num_eigs_idx(block_size_idx)}.exp_params.block_size;
      info.min_num_eigs(iter, block_size_idx) = data.eigs.results{L_idx, noise_ratio_idx, block_size_idx, min_num_eigs_idx(block_size_idx)}.exp_params.num_eigs;
      info.min_num_mat_vec(iter, block_size_idx) = data.eigs.results{L_idx, noise_ratio_idx, block_size_idx, min_num_eigs_idx(block_size_idx)}.num_mat_vec{iter};
   end
end
[opt_mat_vec, opt_block_size_idx] = min(info.min_num_mat_vec, [], 2);





% Plots surfaces of eigenvalue diff avgs
if plot_eig_val_diff
   figure;
   normalize_mat_vecs = false;
   normalize_eig_diffs = false;
   print_dots_on_mins = true;
   eig_start = 2;
   eig_end = 24;
   num_prev_diffs_to_avg = 4; % 0 corresponds to no averaging
   shift_start = 2; % 0 corresponds to no averaging
   shift_end = 2;
   
   iter_start = iter_start_min_mat_vec_2; 
   iter_end = iter_end_min_mat_vec_2; 

   for iter = 1:iter_end
      d = results_EMEP_solvers.eigs{L_idx, noise_ratio_idx, end, end}.d{iter};
      eig_diffs(:, iter) = d(1:end-1) - d(2:end);
   end

   plot_idx = 1;
   for shift_val = shift_start:shift_end

   %      subplot( (shift_end - shift_start + 1)/2, 2, plot_idx)

      % Averages each eig diff with neighboring diffs (smooths plot)
      Z_row = 1;
      Z_col = 1;
      for i = iter_start:iter_end
         for j = eig_start:eig_end
            idx_start = max(eig_start, j-shift_val);
            Z_eig_diff_avgs(Z_row, Z_col) = mean(eig_diffs(idx_start:j, i));
            Z_row = Z_row + 1;
         end
         Z_col = Z_col + 1;
         Z_row = 1;
      end
      if normalize_eig_diffs
         for i = 1 : iter_end - iter_start + 1
            Z_eig_diff_avgs(:, i) = Z_eig_diff_avgs(:, i) / max(Z_eig_diff_avgs(:, i));
         end
      end

      Z_eig_diff_avgs_temp = Z_eig_diff_avgs; Z_eig_diff_avgs = [];
      % Averages current eig diff with previous ones (smooths plot)
      [num_rows, num_cols] = size(Z_eig_diff_avgs_temp);
      for i = 1:num_rows
         for j = 1:num_cols
            idx_start = max(1, j-num_prev_diffs_to_avg);
            Z_eig_diff_avgs(i, j) = mean(Z_eig_diff_avgs_temp(i, idx_start:j));
         end
      end

      [X, Y] = meshgrid(iter_start:iter_end, eig_start+1:eig_end+1);
      surf(X, Y, Z_eig_diff_avgs)
      xlabel('EMEP iterate')
      ylabel('num eigs')   

      title(strcat('avg dist btwn current & ', num2str(shift_val+1), ' prev eigvals'))

      
      % Plots dots for min num mat-vecs for each iterate
      if print_dots_on_mins      
         [X, Y] = meshgrid(iter_start:iter_end, eig_start+1:eig_end+1);
         mat_vec_temp = [];
         mat_vec_temp = squeeze(permute(mat_vec_array(iter_start:iter_end, 2, eig_start:eig_end), [3,2,1]));
         if normalize_mat_vecs
         for i = 1 : iter_end - iter_start + 1
         mat_vec_temp(:,i) = mat_vec_temp(:,i) / max(mat_vec_temp(:,i));
         end
         end

         % Eliminates outliers in plot
         mat_vec_temp  = min(mat_vec_temp, 2000*ones(size(mat_vec_temp)));
         Z_mv = reshape(mat_vec_temp, size(X));
         
         hold on
         for i = 1 : iter_end - iter_start + 1
            [~, idx] = min(Z_mv(:,i));
            val = Z_eig_diff_avgs(idx,i);
            scatter3(iter_start+i-1, eig_start+idx, val, 'filled', 'k')
         end
         hold off
      end

      plot_idx = plot_idx + 1;
   end
end
%}



   


end









