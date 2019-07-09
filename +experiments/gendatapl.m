function [A, b, x0, info] = gendatapl(varargin)
% EXPERIMENT.gendatapl  Generate data for a PhaseLift experiment.
%
% Original code modified to read general images.  ~WW
import util.*
ip = inputParser;
ip.addParameter('signal', 'gaussian');    % can be 'gaussian' or 'image'
% Default image is the "Wusterland2" image from the Hubble Telescope
ip.addParameter('image_file', 'data/nebula.tif');
ip.addParameter('noise_type', 'gaussian');   % can be 'gaussian', 'dual', or 'exponential'
ip.addParameter('mask', 'gaussian');
ip.addParameter('mu', 1); % parameter for exponential distribution of noise term
ip.addParameter('noise_ratio', 0);     % noise in data, factor of norm(b).
ip.addParameter('n', 128);
ip.addParameter('m', 1);
ip.addParameter('L', []);
ip.addParameter('normalize_signal', false); % set true to scale ||x0|| = 1.
ip.addParameter('normalize_mask', false); % set true to scale ||x0|| = 1.
ip.addParameter('regularize_mask', false);
ip.addParameter('seed', 0);
ip.addParameter('resizeImage', 1);
ip.addParameter('channel', 'gray');
ip.parse(varargin{:});
signal = ip.Results.signal;
image_file = ip.Results.image_file;
noise_type = ip.Results.noise_type;
if strcmp(signal, 'dual')
   signal = 'gaussian';
   noise_type = 'dual';
end
n = ip.Results.n; % signal length
m = ip.Results.m; % signal width
L = ip.Results.L; % number of measurements
noise_ratio = ip.Results.noise_ratio;
mask = ip.Results.mask;
normalize_signal = ip.Results.normalize_signal;
normalize_mask = ip.Results.normalize_mask;
regularize_mask = ip.Results.regularize_mask;
mu = ip.Results.mu;

rng(ip.Results.seed); % Set the RNG

% ---------------------------------------------------------------------
% Generate signal
% ---------------------------------------------------------------------
switch signal
   case 'image'
      % Default image is the "Wusterland2" image from the Hubble Telescope:
      % http://hubblesite.org/newscenter/archive/releases/nebula/2015/12/image/a/
      % Use an octanary mask.
      x0 = imread(image_file);
      if ip.Results.resizeImage < 1
         x0 = imresize(x0, ip.Results.resizeImage);
      end
      x0 = mat2gray(x0);
      switch lower(ip.Results.channel)
         case 'red'
            x0 = x0(:,:,1);
         case 'green'
            x0 = x0(:,:,2);
         case 'blue'
            x0 = x0(:,:,3);
         case 'gray'
            x0 = rgb2gray(x0);
      end
      [n, m] = size(x0);
      mask = 'octanary';
      normalize_signal = false;
      
   case 'smooth'
      % Generate a smooth solution. This bit of code is taken from
      %    http://web.stanford.edu/~mahdisol/PRcode.html
      % It generates a random low-pass signal.
      m = 1;
      DFT_matrix = fftshift(fft(eye(n)))./sqrt(n);
      M = ceil(n/8);
      measurements_lowpass = (n/2-round(M/2-1)):(n/2+round(M/2));
      DFT_lowpass = DFT_matrix(measurements_lowpass,:);
      x0 = DFT_lowpass'*(randn(M,1) + 1i*randn(M,1));
      
   case 'gaussian'       % {'gaussian','dual'}
      x0 = randn(n,m) + 1i*randn(n,m);

end

% Ensure signal is a vector (possibly normalized).
x0 = x0(:);
if normalize_signal
   x0 = x0 / norm(x0);
end

% ---------------------------------------------------------------------
% Generate mask.
% ---------------------------------------------------------------------
switch mask

   case 'binary'
      if isempty(L)
         L = 10;
      end
      M = cdp.make('binary', n, m, L, ...
         'forceregular',regularize_mask, 'normalize', normalize_mask);
      
   case 'gaussian'
      if isempty(L)
         L = 10;
      end
      M = cdp.make('cgaussian', n, m, L, ...
         'forceregular',regularize_mask, 'normalize', normalize_mask);

   case 'octanary'
      if isempty(L)
         L = 10;
      end
      M = cdp.make('octanary', n, m, L, ...
         'forceregular',regularize_mask, 'normalize', normalize_mask);
end

% ---------------------------------------------------------------------
% Create operator, measurements (RHS), etc.
% ---------------------------------------------------------------------
A = hop.pl(M,isreal(x0));  % phase-lift operator
b = A.forward(x0);         % measurements

% noise_type = 'dual' condition generates a sythetic model with known dual 
% solution to the image recovery problem.
% This model generation method creates problems that solve much more
% quickly than if observation b = A.forward(x0) has gaussian noise.
if strcmp(noise_type,'dual') && strcmp(signal, 'gaussian') ...
      && (noise_ratio > 0)
   y = randn(size(b));
   A_temp = hop.pl(M,isreal(x0));
   [~,~,V,d] = A_temp.gdobjective(y,1,struct('tol',eps,'p',5,'maxit',512));
   if d(1) < 0
      y = -y;
      [~,~,V,d] = A_temp.gdobjective(y,1,struct('tol',eps,'p',5,'maxit',512));
   end
   if normalize_signal
      x0  = reshape(V(:,1),n,m);
   else
      x0  = reshape(V(:,1)/realsqrt(d(1)),n,m);
   end
   b   = A_temp.forward(x0);
   eta = y/normv(y);
   info.y0 = y;
elseif strcmp(noise_type, 'dual') && ~strcmp(signal, 'gaussian') ...
      && (noise_ratio > 0)
   A_temp = hop.pl(M,isreal(x0));
   fprintf('Generating dual noise parameter.\n');
   fprintf('Note that this can take a while for models with low oversampling (small L).\n')
   % Must solve noiseless PLGD model and recover optimal `y` for synthetic `dual` noise
   [~, ~, solInfo] = saga_sd(A_temp, b, 'verbosity', 0, 'printRelErr', true, 'trueSignalx0', x0);
   fprintf('  Parameter found.\n');
   y = solInfo.y;
   eta = y/normv(y);
   info.y0 = y;
elseif strcmp(noise_type, 'gaussian') && (noise_ratio > 0)
   eta = randn(size(b));
   eta = eta/normv(eta);
elseif strcmp(noise_type, 'exponential') && (noise_ratio > 0)
   eta = exprnd(mu, size(b));
   eta = eta/normv(eta,1);
end

% Add noise.
b0 = b;
if (noise_ratio > 0) && (strcmp(noise_type, 'gaussian') || strcmp(noise_type, 'dual'))
   e       = noise_ratio;
%   c       = [rdot(eta*(1+e),eta*(1-e)),...
%                    -2*rdot(e*b0,e*eta),...
%                       -rdot(e*b0,e*b0)];
   c      = [1-e^2, -2*e^2*rdot(b0,eta),...
                       -e^2*rdot(b0,b0)];
   epsilon = max(real(roots(c)));
%    noise_ratio = e*normv(b0);
   b = b0+epsilon*eta;
elseif (noise_ratio > 0) && strcmp(noise_type, 'exponential')
   % This computation assumes both b0 >= 0 and eta >= 0
   e       = noise_ratio;
   epsilon = e*normv(b0,1) / (1-e);
   b = b0+epsilon*eta;
else
   epsilon = 0;
end


% Save stats to info
info.epsilon = epsilon;
info.b0 = b0;

end
