% Clear workspace, command window, and add required paths
clear;
clc;

addpath('./CDCS-master/packages/+cdcs_utils');
addpath('./utilities/BundleSDP/');
addpath('./utilities/General/')


% Set problem parameters
opts.numOfConstraints = 100;
opts.seed = 123;
opts.dynamicAlpha = true;
opts.alpha = 0.5;
opts.maxIter = 25;
opts.eps = 10^-5;
opts.beta = 0.1;
opts.mu = 0.5;
opts.ml = 0.01;
opts.r_p = 0;
opts.r_c = 4;

% Set data parameters for the SDP problem
nCones = 1;               % Number of cones with block-arrow sparsity pattern
nBlk = [1];              % Number of diagonal blocks for each cone
BlkSize = [40];            % Block size for each cone
ArrowHead = [2];          % Arrow head size for each PSD cone
iternconnection = [0];    % Number of overlapping elements for neighboring cones

fprintf('\nSetting up random block-arrow SDP with n = %i, m = %i  ... \n',BlkSize*nBlk+ArrowHead, opts.numOfConstraints);
tsetup = tic;

% Generate data for the SDP problem
[At, b, c, K, sparsityPat, X] = GenStrictData(opts.numOfConstraints, nCones, ...
    nBlk, BlkSize, ArrowHead);
opts.sparsityPat = sparsityPat;
opts.dim = length(opts.sparsityPat);
A = reshape(full(At), [opts.dim, opts.dim, opts.numOfConstraints]);
B = full(b);
C = reshape(full(c), [opts.dim, opts.dim]);

% Initialize X with the generated data
X = mat(full(X{1}));
X(end - ArrowHead + 1:end, end - ArrowHead + 1:end) = zeros(ArrowHead);

opts.At = At;
opts.b = b;
opts.c = c;

% Calculate P_star and rho
P_star = trace(C' * X);
opts.rho = trace(C);

% Solve the SDP problem using Bundle Decomposition
[X_star, count, a] = BundleDecomposedSDP(A, B, C, opts);

% Plot convergence
semilogy(1:opts.maxIter-1, abs((a.Obj(2:end) - P_star)) / abs(P_star))

% Store optimization results
Out.gap = a.gap;
Out.costOpt = (a.Obj(end) - P_star) / abs(P_star);
Out.affineOpt = a.dualAffine;
Out.P_star = P_star;
Out.obj = a.Obj;

