% test_unittest
%
% Unit tests for accuracy of functions for computing central and
% off-diagonal blocks of the inverse of a symmetric or  non-symmetric
% block-tridiagonal matrix, as implemented in:  
%
% 1. invblktridiag
% 2. invblktridiag_sym

N = 5; % size of blocks
T = 20; % number of blocks
TOL = 1e-10;  % tolerance for numerical error

%% =================================================
% Test for non-symmetric block-tridiagonal matrices
% =================================================

% 1. Build non-symmetric block-tridiagonal matrix
% --------------------------------------------

% Make diagonal blocks
BdiagMats = mat2cell(randn(N,N*T),N,N*ones(1,T)); % diagonal blocks
Bmat = blkdiag(BdiagMats{:});

% Make below-diagonal and above-diagonal blocks
Bdiag1Mats = mat2cell(randn(N,N*(T-1)),N,N*ones(1,T-1)); % diagonal blocks
Bdiag2Mats = mat2cell(randn(N,N*(T-1)),N,N*ones(1,T-1)); % diagonal blocks
Bdiag1 = blkdiag(Bdiag1Mats{:});
Bdiag2 = blkdiag(Bdiag2Mats{:});

% Assemble block-tridiagonal matrix
Bmat(N+1:end,1:end-N) = Bmat(N+1:end,1:end-N) + Bdiag1; % below-diagonal block
Bmat(1:end-N,N+1:end) = Bmat(1:end-N,N+1:end) + Bdiag2; % below-diagonal block

%% 2. compute explicit inverse, non-symmetric matrix

Binv = inv(Bmat);

% Use our function
[vdiag1,B1,B2] = invblktridiag(Bmat,N);

%% 3. Test for agreement

% Check diagonal elements
vdiag0 = diag(Binv); % diagonal elements of explicit inverse
diagdiff = sum(abs(vdiag1-vdiag0));
if diagdiff>TOL
    warning('invblktridiag unit test FAILED: diagonals don''t match');
else
    fprintf('invblktridiag unit test PASSED: diagonals match\n');
end

% Check diagonal blocks
MatchFAIL = 0;
for j_ind = 1:T
    ii = (j_ind-1)*N+1:j_ind*N;  % row/column indices
    aa = Binv(ii,ii);   % block from direct inverse 
    bb = B1(:,:,j_ind); % computed block
    errs = (abs((aa(:)-bb(:))));
    if any(errs>TOL)
        MatchFAIL = 1;
    end
end
if MatchFAIL
    warning('invblktridiag unit test FAILED: diagonal blocks don''t match');
else
    fprintf('invblktridiag unit test PASSED: diagonal blocks match\n');
end

% Check above-diagonal blocks
MatchFAIL = 0;
for j_ind = 1:T-1
    ii = (j_ind-1)*N+1:j_ind*N; % rows
    jj = j_ind*N+1:(j_ind+1)*N; % columns
    aa = Binv(ii,jj);  % block from direct inverse 
    bb = B2(:,:,j_ind);
    errs = (abs((aa(:)-bb(:))));
    if any(errs>TOL)
        MatchFAIL = 1;
    end
end
if MatchFAIL
    warning('invblktridiag unit test FAILED: above-diagonal blocks don''t match');
else
    fprintf('invblktridiag unit test PASSED: above-diagonal blocks match\n');
end
    
%% =================================================
% Now test for symmetric block-tridiagonal matrices
% ==================================================

% 1. Make symettric block-tridiagonal
Bmat = Bmat + Bmat';  % make symmetric matrix

% 2.  Compute explicit inverse
Binv = inv(Bmat);

% 3. Compute blocks of inverse using our function
[vdiag1,B1,B2] = invblktridiag_sym(Bmat,N);

% 4. Check for agreement

% Check diagonal elements
vdiag0 = diag(Binv); % diagonal elements of explicit inverse
diagdiff = sum(abs(vdiag1-vdiag0));
if diagdiff>TOL
    warning('invblktridiag_sym unit test FAILED: diagonals don''t match');
else
    fprintf('invblktridiag_sym unit test PASSED: diagonals match\n');
end

% 4a. Check diagonal blocks
MatchFAIL = 0;
for j_ind = 1:T
    ii = (j_ind-1)*N+1:j_ind*N;  % row/column indices
    aa = Binv(ii,ii);   % block from direct inverse 
    bb = B1(:,:,j_ind); % computed block
    errs = (abs((aa(:)-bb(:))));
    if any(errs>TOL)
        MatchFAIL = 1;
    end
end
if MatchFAIL
    warning('invblktridiag_sym unit test FAILED: diagonal blocks don''t match');
else
    fprintf('invblktridiag_sym unit test PASSED: diagonal blocks match\n');
end

% 4b. Check above-diagonal blocks
MatchFAIL = 0;
for j_ind = 1:T-1
    ii = (j_ind-1)*N+1:j_ind*N; % rows
    jj = j_ind*N+1:(j_ind+1)*N; % columns
    aa = Binv(ii,jj);  % block from direct inverse 
    bb = B2(:,:,j_ind);
    errs = (abs((aa(:)-bb(:))));
    if any(errs>TOL)
        MatchFAIL = 1;
    end
end
if MatchFAIL
    warning('invblktridiag_sym unit test FAILED: above-diagonal blocks don''t match');
else
    fprintf('invblktridiag_sym unit test PASSED: above-diagonal blocks match\n');
end
