function [MinvDiag,MinvDiagBlocks,MinvAboveDiagBlocks] = invblktridiag(M,nn)
% [MinvDiag,MinvDiagBlocks,MinvAboveDiagBlocks] = invblktridiag(M,nn)
%
% Find main and above-diagonal blocks of the inverse of a block tri-diagonal 
% matrix using O(N) recursive method of Rybicki & Hummer (1991). 
%
% Inputs:  
% ------
%   M = block-tridiagonal matrix
%  nn = size of blocks
%
% Output: 
% -------
%  MinvDiag - vector with diagonal elements of inverse
%  MinvDiagBlocks - central blocks of inverse (nn x nn x T tensor)
%  MinvAboveDiagBlocks - above-diagonal blocks of inverse (nn x nn x T-1 tensor)
%
% NOTE: doesn't yet compute the below-diagonal blocks!

nblocks = size(M,1)/nn; % number of total blocks

% Matrices to store during recursions
A = zeros(nn,nn,nblocks); % for below-diagonal blocks
B = zeros(nn,nn,nblocks); % for diagonal blocks
C = zeros(nn,nn,nblocks); % for above-diagonal blocks
D = zeros(nn,nn,nblocks); % quantity to compute
E = zeros(nn,nn,nblocks); % quantity to compute

% Initialize first D block
inds = 1:nn; % indices for 1st block
B(:,:,1) = M(inds,inds);
C(:,:,1) = M(inds,inds+nn); 
D(:,:,1) = B(:,:,1)\C(:,:,1);

% Initialize last E block
inds = (nblocks-1)*nn+(1:nn);  % indices for last block
A(:,:,nblocks) = M(inds,inds-nn);
B(:,:,nblocks) = M(inds,inds);
E(:,:,nblocks) = B(:,:,nblocks)\A(:,:,nblocks);

% Extract blocks A, B, and C
for ii = 2:nblocks-1
    inds = (ii-1)*nn+1:ii*nn; % indices for center block
    A(:,:,ii) = M(inds,inds-nn); % below-diagonal block
    B(:,:,ii) = M(inds,inds); % middle diagonal block
    C(:,:,ii) = M(inds,inds+nn); % above diagonal block
end
    
% Make a pass through data to compute D and E
for ii = 2:nblocks-1
    % Forward recursion
    D(:,:,ii) = (B(:,:,ii)-A(:,:,ii)*D(:,:,ii-1))\C(:,:,ii); 
    
    % backward recursion
    jj = nblocks-ii+1;
    E(:,:,jj) = (B(:,:,jj)-C(:,:,jj)*E(:,:,jj+1))\A(:,:,jj); 
end
    
% Now form blocks of inverse covariance
I = eye(nn);
MinvDiagBlocks = zeros(nn,nn,nblocks);
MinvAboveDiagBlocks = zeros(nn,nn,nblocks-1);
MinvDiagBlocks(:,:,1) = (B(:,:,1)*(I-D(:,:,1)*E(:,:,2)))\I;
MinvDiagBlocks(:,:,nblocks) = (B(:,:,nblocks)-A(:,:,nblocks)*D(:,:,nblocks-1))\I;
for ii = 2:nblocks-1
    % compute diagonal blocks of inverse
    MinvDiagBlocks(:,:,ii) = ((B(:,:,ii)-A(:,:,ii)*D(:,:,ii-1))*(I-D(:,:,ii)*E(:,:,ii+1)))\I;
    % compute above-diagonal blocks
    MinvAboveDiagBlocks(:,:,ii-1) = -(D(:,:,ii-1)*MinvDiagBlocks(:,:,ii));
end
MinvAboveDiagBlocks(:,:,nblocks-1) = -(D(:,:,nblocks-1)*MinvDiagBlocks(:,:,nblocks)); 

% Extract just the diagonal elements
MinvDiag = zeros(nn*nblocks,1);
for ii = 1:nblocks
    MinvDiag((ii-1)*nn+1:ii*nn,1) = diag(MinvDiagBlocks(:,:,ii));
end



