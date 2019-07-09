function [Aty] = getfulldualmatrix(A, y)
% Computes the dual product as full matrix: At(y)
%
% TODO: Write test for this method

if A.m2 == 1
   y_array = reshape(y, [A.n, 1, A.m3]);
   Aty = sum(bsxfun(@times,conj(A.masks), ...
                  ifft2(bsxfun(@times,y_array, ...
                  fft2(bsxfun(@times,A.masks, ...
                  eye(A.n, A.n)))))),3);
else
   y_array = y(:);
   Aty = zeros(A.n, A.n);
   for i = 1:length(y_array)
      submask_idx = mod(i, A.n);
      if submask_idx == 0
      submask_idx = A.n;
      end
      mask_num = (i - submask_idx)/A.n + 1;
      ei = zeros(A.n, 1);
      ei(submask_idx, 1) = 1;
      ei = reshape(ei, [A.m1, A.m2]);
      ai = sqrt(A.n)*bsxfun(@times, conj(A.masks(:, :, mask_num)), ifft2(ei));
      Ai = ai(:)*ai(:)';
      Aty = Aty + y(i)*Ai;
   end
end

%{
Ak_test = zeros(A.n, A.n);
I = eye(A.n, A.n);
yStart = 1;
yEnd = A.n;
for l = 1:10
   FC = sqrt(A.n)*diag(vec(conj(A.masks(:, :, l))))*ifft2(I);
   %FC = sqrt(A.n)*bsxfun(@times, conj(A.masks(:, :, l)), ifft2(I));
   Ak_test = Ak_test + FC*diag(y_array(yStart:yEnd))*FC';
   yStart = yStart + A.n;
   yEnd = yEnd + A.n;
end

norm(Ak_test - Ak)
%}
end