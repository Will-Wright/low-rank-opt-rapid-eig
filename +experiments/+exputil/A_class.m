classdef A_class < handle
   properties
      A_data
      num_mat_vec_mults
      num_calls
   end
   methods
      function self = A_class(A_data)
         self.A_data = A_data;
         self.num_mat_vec_mults = [];
         self.num_calls = [];
      end
      
      function [Ax] = A_fun(self, iter, X)
         % function [Ax] = A_fun(A_data, iter, X)
         %
         % This function computes matrix-vector products based on the
         % adjoint of the PhaseLift sensing operator A.
         %
         % This operator has forward and adjoint maps.  The forward map
         % is from the lifted object domain to the Fourier domain.  
         % A : C^(n x n) ---> R^m
         %             x |--> b
         %
         % The adjoint map is from the Fourier domain to the object domain.
         % A' : R^m ---> C^(n x n)
         %        y |--> A'(y)
         % where A'(y) = sum_{k = 1:L}  C_k' F' Diag(y_k) F C_k,
         % with C_k are sensing masks, F is the discrete Fourier transform,
         % and y_k is a subsegment of the dual iterate y.
         %
         % Since C_k and Diag(y_k) are diagonal, and F is computed with the fast 
         % Fourier transform, matrix-vector products [A'(y)]*x are computed 
         % implicitly, without forming A'(y).

         [X_rows, X_cols] = size(X);
         y = self.A_data.y_matrix(:, iter);

         X   = reshape(X,[self.A_data.n1, self.A_data.n2, 1, numel(X)/self.A_data.n, 1]);
         y   = reshape(y,[self.A_data.m1, self.A_data.m2, self.A_data.m3, 1, numel(y)/self.A_data.m]);
         Ax = sum(bsxfun(@times, conj(self.A_data.masks), ...
                         ifft2(bsxfun(@times, y, ...
                                      fft2(bsxfun(@times, self.A_data.masks, X))))), 3);

         if self.A_data.isreal
            Ax = reshape(real(Ax(:)), [X_rows, X_cols]);
         else
            Ax = reshape(Ax(:), [X_rows, X_cols]);
         end
         if length(self.num_mat_vec_mults) < iter || length(self.num_calls) < iter
            self.num_mat_vec_mults(iter) = 0;
            self.num_calls(iter) = 0;
         end
         self.num_mat_vec_mults(iter) = self.num_mat_vec_mults(iter) + numel(X)/self.A_data.n;
         self.num_calls(iter) = self.num_calls(iter) + 1;
      end
   end
end

