%{
    Tests nuclear norm minimization,

    min_X ||X||_*
s.t.
    P(X) == b, or || P(X) - b || <= eps

where X is a M x N matrix, and P is an observation operator

The solvers solve a regularized version, using
||X||_* + mu/2*||X-X_0||_F^2

THIS VERSION TESTS THE PACKSVD FORM

%}
nuclear_norm = @(x) sum(svd(x,'econ'));

% Try to load the problem from disk
fileName = 'nuclearNorm_problem1_noiseless';
if exist([fileName,'.mat'],'file')
    load(fileName);
    fprintf('Loaded problem from %s\n', fileName );
else
    
    % Generate a new problem

    randn('state',sum('NuclearNorm'));
    rand('state',sum('NuclearNorm2'));
    
    M = 30; N = 30; R = 2;

    df = R*(M+N-R);
    oversample = 5;
    Left  = randn(M,R);
    Right = randn(N,R);
    k = round(oversample*df); 
    k = min( k, round(.8*M*N) );
    omega = randperm(M*N);
    omega = sort(omega(1:k)).';

    X_original = Left*Right';       % the "original" signal -- may not be optimal value
    b_original = X_original(omega); 
    EPS = 0;        % noise level
    b = b_original;
    obj_original = nuclear_norm(X_original);

    % get reference via CVX
    cvx_begin
        cvx_precision best
        variable Xcvx(M,N)
        minimize norm_nuc(Xcvx)
        subject to
            Xcvx(omega) == b
    cvx_end
    X_reference = Xcvx;         % the nuclear norm minimizer
    obj_reference = nuclear_norm(X_reference);
    
    save(fileName,'X_original','X_reference','omega','b','obj_original',...
        'Left','Right','EPS','b_original','R','obj_reference');
    fprintf('Saved data to file %s\n', fileName);
    
end

[M,N]           = size(X_reference);
norm_X_reference = norm(X_reference,'fro');
er_reference    = @(x) norm(x-X_reference,'fro')/norm_X_reference;
resid           = @(x) norm(x(omega)-b)/norm(b);  % change if b is noisy

[omegaI,omegaJ] = ind2sub([M,N],omega);
mat = @(x) reshape(x,M,N);
vec = @(x) x(:);
    
k  = length(omega);
p  = k/(M*N);
df = R*(M+N-R);
fprintf('%d x %d rank %d matrix, observe %d = %.1f x df = %.1f%% entries\n',...
    M,N,R,k,k/df,p*100);
fprintf(' Nuclear norm solution and original matrix differ by %.2e\n',...
    norm(X_reference-X_original,'fro')/norm_X_reference );

%% set some parameters

testOptions = {};
opts = [];
opts.alg    = 'AT';
% opts.alg    = 'GRA';
% opts.errFcn = @(f,dual,x) norm( x -
% X_reference,'fro')/norm(X_reference,'fro');
nrm = @(X) sqrt( tfocs_normsq(X) );
opts.errFcn = @(f,dual,x) nrm( x - X_reference)/nrm(X_reference);

mu      = .001;
% X0      = 0; 
X0      = packSVD(M,N);

opts.largeScale = true;
% opts.largeScale = false;

% [ x, out, opts ] = solver_sNuclearBP_packSVD( {M,N,omega}, b, mu, X0, [], opts ); 
[ x, out, opts ] = solver_sNuclearBP( {M,N,omega}, b, mu, X0, [], opts ); 

%% for comparison, try it without packSVD
[ x, out, opts ] = solver_sNuclearBP( {M,N,omega}, b, mu, X0, [], opts ); 

