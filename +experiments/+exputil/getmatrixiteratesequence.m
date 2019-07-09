function [yCell, eigsTolCell, infoEigsFFTsCell, dCell] = getmatrixiteratesequence(results)
% Concatenates y values from linesearch, linesearch result, and refinement
%
% Output:
%     yCell{i}{j}    cell of dual variables y, with i = main saga_sd iteration
%                    and j = subroutine iterates in iterate i
%
% Note: saga_sd.solInfo.iterData.y{k} = saga_sd.solInfo.iterData.yLinesearch{1, k}{1, 1}
% for NOISY PROBLEMS, but NOT EQUAL FOR NOISELESS
%

dCell = cell(1);
yCell = cell(1);
eigsTolCell = cell(1);
infoEigsFFTsCell = cell(1);

nitn = results.saga_sd.solInfo.nit;
iterData = results.saga_sd.solInfo.iterData;
dCell{1}{1} = iterData.d{1};
yCell{1}{1} = iterData.y{1};
eigsTolCell{1}{1} = iterData.eigsTol{1};
if ~isempty(iterData.yRefined{1})
   yCell{1}{2} = iterData.yRefined{1};
   dCell{1}{2} = iterData.dRefined{1};
   eigsTolCell{1}{2} = iterData.eigsTolRefined{1};
   infoEigsFFTsCell{1}{2} = iterData.infoEigsFFTs_yRefined{1}{1};
end
% saga_sd iterations
for itnSaga = 1:nitn
   if ~isempty(iterData.yLinesearch{itnSaga+1})
      for itnLS = 1:length(iterData.yLinesearch{itnSaga+1})
         yCell{itnSaga+1}{itnLS} = iterData.yLinesearch{itnSaga+1}{itnLS};
         dCell{itnSaga+1}{itnLS} = iterData.dLinesearch{itnSaga+1}{itnLS};
         eigsTolCell{itnSaga+1}{itnLS} = iterData.eigsTolLinesearch{itnSaga+1}{itnLS};
         infoEigsFFTsCell{itnSaga+1}{itnLS} = iterData.infoEigsFFTs_yLS{itnSaga+1}{itnLS};
      end
   else
      itnLS = 1;
      yCell{itnSaga+1}{itnLS} = iterData.y{itnSaga+1};
      dCell{itnSaga+1}{itnLS} = iterData.d{itnSaga+1};
      eigsTolCell{itnSaga+1}{itnLS} = iterData.eigsTol{itnSaga+1};
      infoEigsFFTsCell{itnSaga+1}{itnLS} = iterData.infoEigsFFTs_y{itnSaga+1}{1,1};
   end
   if ~isempty(iterData.yRefined{1})
      yCell{itnSaga+1}{itnLS+1} = iterData.yRefined{itnSaga+1};
      dCell{itnSaga+1}{itnLS+1} = iterData.dRefined{itnSaga+1};
      eigsTolCell{itnSaga+1}{itnLS+1} = iterData.eigsTolRefined{itnSaga+1};
      infoEigsFFTsCell{itnSaga+1}{itnLS+1} = iterData.infoEigsFFTs_yRefined{itnSaga+1}{1,1};
   end
end


end

