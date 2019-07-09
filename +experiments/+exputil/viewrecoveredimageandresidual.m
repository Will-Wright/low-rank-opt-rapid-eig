function [imSagaRes, imWflowRes, imOrig, imSaga] = viewrecoveredimageandresidual(results, varargin)
% EXPERIMENT.viewrecoveredimageandresidual
%
% Plots figures of original image (x0), recovered image (saga_sd.x), and
% residual image (abs(x0 - saga_sd.x))

ip = inputParser;
ip.addParameter('iterate', []);
ip.parse(varargin{:});

iterate = ip.Results.iterate;

imSagaRes = [];
imWflowRes = [];
imOrig = [];
imSaga = [];

if isfield(results, 'red') && isfield(results, 'green') ...
         && isfield(results, 'blue')
   saga_is_color = true;
else
   saga_is_color = false;
end


nPlotRows = 0;
if isfield(results, 'saga_sd') ...
      || (isfield(results, 'red') && isfield(results, 'green') && isfield(results, 'blue'))
   saga_sdIsField = true;
   wflowDisplayPosition = 3;
   if saga_is_color
      imageIsGrayscale = false;
      nPlotRows = nPlotRows + 1;
   elseif isfield(results.saga_sd, 'x')
      imageIsGrayscale = true;
      nPlotRows = nPlotRows + 1;
   else
      fprintf('Invalid results struct.  No color or grayscale signal provided.\n');
      return;
   end
else
   saga_sdIsField = false;
   wflowDisplayPosition = 0;
end

if (~saga_is_color && isfield(results, 'wflow')) ...
      || (saga_is_color && isfield(results.red, 'wflow') && isfield(results.green, 'wflow') && isfield(results.blue, 'wflow'))
   wflowIsField = true;
   if (saga_is_color && isfield(results.red, 'wflow') && isfield(results.green, 'wflow') && isfield(results.blue, 'wflow'))
      imageIsGrayscale = false;
      nPlotRows = nPlotRows + 1;
   elseif isfield(results.wflow, 'x')
      imageIsGrayscale = true;
      nPlotRows = nPlotRows + 1;
   else
      fprintf('Invalid results struct.  No color or grayscale signal provided.\n');
      return;
   end
else
   wflowIsField = false;
end

figure;
hold on;

% Displays image of original image, saga_sd result, and residual
if saga_sdIsField
   if imageIsGrayscale
      [n1, n2] = size(results.saga_sd.x);
      imSize = [n1 n2];
      n = n1*n2;
%      figure;
      x0Mat = real(reshape(results.gendatapl.x0, imSize));
      imOrig = mat2gray(x0Mat);
      subplot(nPlotRows,3,1), imshow(imOrig); title('Original image');
%      hold on;

      x_saga = results.saga_sd.x;
      % Identifies correct sign of result 'x' and removes negative pixel values
      realX = real(reshape(x_saga, imSize));
      idxMatPos = realX > 0;
      xSagaMat = zeros(imSize);
      if sum(sum(idxMatPos)) >= 0.5*n
         xSagaMat(idxMatPos) = realX(idxMatPos);
      else
         xSagaMat(~idxMatPos) = -realX(~idxMatPos);
      end
      imSaga = mat2gray(xSagaMat);
      
      subplot(nPlotRows,3,2), imshow(imSaga); title('Recovered with saga\_sd');
      imSagaRes = mat2gray(abs(xSagaMat - x0Mat));
      subplot(nPlotRows,3,3), imshow(imSagaRes); title('saga\_sd residual');

   else
      [n1, n2] = size(results.red.saga_sd.x);
      imSize = [n1 n2];
      n = n1*n2;
%      figure;
      imOrig = cat(3, ...
         reshape(results.red.gendatapl.x0, imSize), ...
         reshape(results.green.gendatapl.x0, imSize), ...
         reshape(results.blue.gendatapl.x0, imSize)); 
      subplot(nPlotRows,3,1), image(imOrig); title('Original image');
%      hold on;

      % Identifies correct sign of result 'x' and removes negative pixel values
      realX = real(reshape(results.red.saga_sd.x, imSize));
      idxMatPos = realX > 0;
      xSagaMatRed = zeros(imSize);
      if sum(sum(idxMatPos)) >= 0.5*n
         xSagaMatRed(idxMatPos) = realX(idxMatPos);
      else
         xSagaMatRed(~idxMatPos) = -realX(~idxMatPos);
      end
      realX = real(reshape(results.green.saga_sd.x, imSize));
      idxMatPos = realX > 0;
      xSagaMatGreen = zeros(imSize);
      if sum(sum(idxMatPos)) >= 0.5*n
         xSagaMatGreen(idxMatPos) = realX(idxMatPos);
      else
         xSagaMatGreen(~idxMatPos) = -realX(~idxMatPos);
      end
      realX = real(reshape(results.blue.saga_sd.x, imSize));
      idxMatPos = realX > 0;
      xSagaMatBlue = zeros(imSize);
      if sum(sum(idxMatPos)) >= 0.5*n
         xSagaMatBlue(idxMatPos) = realX(idxMatPos);
      else
         xSagaMatBlue(~idxMatPos) = -realX(~idxMatPos);
      end

      imSaga = cat(3, xSagaMatRed, xSagaMatGreen, xSagaMatBlue); 

      subplot(nPlotRows,3,2), image(imSaga); title('Recovered with saga\_sd');
      imSagaRes = abs(imOrig - imSaga);
      subplot(nPlotRows,3,3), image(imSagaRes); title('saga\_sd residual');
   end
end


% Displays image of original image, wflow result, and residual
if wflowIsField
   if imageIsGrayscale
      n1 = results.gendatapl.A.n1;
      n2 = results.gendatapl.A.n2;
      n = n1*n2;
      imSize = [n1 n2];
%      figure;
      x0Mat = real(reshape(results.gendatapl.x0, imSize));
      imOrig = mat2gray(x0Mat);
      subplot(nPlotRows,3,1+wflowDisplayPosition), imshow(imOrig); title('Original image');
%      hold on;

      xWflowMat = abs(real(reshape(results.wflow.x, imSize)));
      
      imWflow = mat2gray(xWflowMat);
      subplot(nPlotRows,3,2+wflowDisplayPosition), imshow(imWflow); title('Recovered with wflow');
      imWflowRes = mat2gray(abs(xWflowMat - x0Mat));
      subplot(nPlotRows,3,3+wflowDisplayPosition), imshow(imWflowRes); title('wflow residual');

   else
      n1 = results.red.gendatapl.A.n1;
      n2 = results.red.gendatapl.A.n2;
      imSize = [n1 n2];
%      figure;
      imOrig = cat(3, ...
         reshape(results.red.gendatapl.x0, imSize), ...
         reshape(results.green.gendatapl.x0, imSize), ...
         reshape(results.blue.gendatapl.x0, imSize)); 
      subplot(nPlotRows,3,1+wflowDisplayPosition), image(imOrig); title('Original image');
%      hold on;

      imWflow = cat(3, ...
         reshape(abs(results.red.wflow.x), imSize), ...
         reshape(abs(results.green.wflow.x), imSize), ...
         reshape(abs(results.blue.wflow.x), imSize));

      subplot(nPlotRows,3,2+wflowDisplayPosition), image(imWflow); title('Recovered with wflow');
      imWflowRes = abs(imOrig - imWflow);
      subplot(nPlotRows,3,3+wflowDisplayPosition), image(imWflowRes); title('wflow residual');
   end
end

end