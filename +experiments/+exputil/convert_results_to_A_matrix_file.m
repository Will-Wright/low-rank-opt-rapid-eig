function [A_data] = convert_results_to_A_matrix_file(results, varargin)

if nargin >= 2
   save_file_name = varargin{1};
   pathName = 'cache/';
   str_end = save_file_name(end-3:end);
   if ~strcmp(str_end, '.mat')
      save_file_name = strcat(save_file_name, '.mat');
   end
   matFile = fullfile(pathName, save_file_name);
end

A = results.gendatapl.A;
iter_range = 1:results.saga_sd.solInfo.nit + 1;
[y_cell, eigsTolCell, ~, d_cell] = experiments.exputil.getmatrixiteratesequence(results);
y_matrix = [];
eig_vals = [];
eigs_tols = [];
y_is_backtrack_step = [];
i = 1;
row = 1;
for col = iter_range
   for sub_iter = 1:length(y_cell{row, col})
      eig_vals(:, i) = d_cell{row, col}{row, sub_iter};
      y_matrix(:, i) = y_cell{row, col}{row, sub_iter}(:);
      eigs_tols(i) = eigsTolCell{row, col}{row, sub_iter};
      if sub_iter == 1
         y_is_backtrack_step(i) = 0;
      else
         y_is_backtrack_step(i) = 1;
      end
      i = i + 1;
   end
end

A_data = struct('isreal', A.isrealin, ...
   'm1', A.m1, 'm2', A.m2, 'm3', A.m3, 'n1', A.n1, 'n2', A.n2, ...
   'n', A.n, 'm', A.m, 'masks', A.masks, ...
   'y_matrix', y_matrix, 'eigs_tols', eigs_tols, ...
   'y_is_backtrack_step', y_is_backtrack_step, ...
   'eig_vals', eig_vals);


%A_data.eig_vals = results.saga_sd.solInfo.iterData.d;
if exist('save_file_name')
   save(matFile, 'A_data');
end
   
end