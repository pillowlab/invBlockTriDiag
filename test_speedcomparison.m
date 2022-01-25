% test_speedcomparison.m
%
% Compare non-symmetric and symmetric versions of code

N = 50; % size of blocks
T = 1e4; % number of blocks

% =====  Build sparse block-tridiagonal matrix =============
fprintf('Building sparse block-tridiagonal matrix...\n');
tic;
% Build a random bi-diagonal matrix
[jj,ii] = meshgrid(1:N);
ii = ii(:); % convert to column vector
jj = jj(:); % convert to column vector

% First, compute indices for block bi-diagonal matrix
nelts = N^2*(T*2-1);  % number of elements
Iinds = zeros(nelts,1); % row indices
Jinds = zeros(nelts,1); % column indices
irow = 0;
for k = 1:T-1
    Iinds(irow+1:irow+(N^2*2)) = [ii+(k-1)*N; ii+k*N];
    Jinds(irow+1:irow+(N^2*2)) = [jj+(k-1)*N; jj+(k-1)*N];
    irow = irow+N^2*2;
end
Iinds(irow+1:irow+N^2) = ii+(T-1)*N; % last block
Jinds(irow+1:irow+N^2) = jj+(T-1)*N; % last block
D = sparse(Iinds,Jinds,randn(length(Iinds),1),T*N,T*N); % fill block bi-diag w/ random entries

% Now make symmetric sparse tri-diagonal matrix
B = D'*D+speye(N*T,N*T)*.1; 

% Compute total build time
buildTime = toc

% Make plot showing block structure of L
ni = min(1000,nelts); % # of rows+cols to plot
clf; spy(B(1:ni,1:ni));
title('symmetric block-tridiagonal matrix');
drawnow;

%% ====  Compare timing of full and symmetric versions ============

tic;
[vdiag1,Dblocks1,OffDblocks1] = invblktridiag(B,N);
computeTime1 = toc;

tic;
[vdiag2,Dblocks2,OffDblocks2] = invblktridiag_sym(B,N);
computeTime2 = toc;

fprintf('Compute time using invblktridiag:     %.2f\n',computeTime1);
fprintf('Compute time using invblktridiag_sym: %.2f\n',computeTime2);
