classdef hop < handle
% hop class definition (Michael Friedlander, 2016)
% Modified to track runtimes (Will Wright, 2017)


   properties
      isrealin
      nForward
      tForward
      nAdjoint
      tAdjoint
      nGDObjective
      tGDObjective
      nEigs
      tEigs
      nMultsInEigs
      nFFTsInEigs
      nPFDObjective
      tPFDObjective
      nDFPObjective
      tDFPObjective
   end % properties
   
   methods ( Abstract )
      out       = forward(self,x1,x2,scatter)
      out       = adjoint(self,y,x)
      [f,g,V,d] = gdobjective(self,y,k,opts,v)
      [f,g,r]   = pfdobjective(self,rhs,x)
      [f,g,r]   = dfpobjective(self,xlhs,xrhs,y)
      sz        = size(self,dim)
   end % methods ( Abstract )
   
   methods

      function self = hop(isrealin)
         if (nargin == 0) || isempty(isrealin) || ~islogical(isrealin)
            isrealin = false;
         end
         self.isrealin      = isrealin;
         self.nForward      = 0;
         self.tForward      = 0.0;
         self.nAdjoint      = 0;
         self.tAdjoint      = 0.0;
         self.nGDObjective  = 0;
         self.tGDObjective  = 0.0;
         self.nEigs         = 0;
         self.tEigs         = 0.0;
         self.nMultsInEigs  = 0;
         self.nFFTsInEigs   = 0;
         self.nPFDObjective = 0;
         self.tPFDObjective = 0.0;
         self.nDFPObjective = 0;
         self.tDFPObjective = 0.0;
      end % hop (constructor)
      
   end % methods
   
end % classdef hop < handle
