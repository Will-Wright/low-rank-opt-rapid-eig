function [spectrum] = getspectrum(results, itn)
import util.*;

A = results.gendatapl.A;

% Caches true solution to eigenvalue problems 
eigsopts.maxit = 5000;
eigsopts.issym = true;
eigsopts.isreal = A.isrealin;
eigsopts.returnIterData = true;

if itn == 1
   y = results.saga_sd.solInfo.iterData.y{1, itn};
else
   % iterates through linesearch vecs for `y`
   LS_idx = 2;
   for i = 2:results.saga_sd.solInfo.nit
      for j = 1:length(results.saga_sd.solInfo.iterData.yLinesearch{1, i})
         if LS_idx == itn            
            y = results.saga_sd.solInfo.iterData.yLinesearch{1, i}{j};
         end
         LS_idx = LS_idx + 1;
      end      
   end
end


% Computes spectrum of Aty with eigs

Afun_temp = @(x) real(vec(A.adjoint(y,x)));
k1 = ceil(A.n/2);
k2 = A.n - k1;

sigma = 'LA';
[V1, d1, ~, iterData] = eigs(Afun_temp, A.n, k1, sigma, eigsopts);
[d1, permutation] ...
 = sort(real(diag(d1)),'descend');
V1 = V1(:, permutation);

sigma = 'SA';
[V2, d2, ~, iterData] = eigs(Afun_temp, A.n, k2, sigma, eigsopts);
[d2, permutation] ...
 = sort(real(diag(d2)),'descend');
V2 = V2(:, permutation);

spectrum = [d1; d2];



end % end function getspectrum