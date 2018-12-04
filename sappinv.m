%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Aleksandar Haber, Last update: December 2014
% Description:
% * MATLAB function that computes an approximate (sparse) inverse of a sparse
%   matrix.
% * Notation: 
%   M is an approximate inverse
%   A is the original matrix whose inverse is computed, M\approx A^{-1}
%   Apriori is the matrix whose non-zero entries describe the a priori sparisty pattern of M
% * The matrix M is computed by solving the global optimization problem  
%   by splitting it into a series of smaller least squares problems. For more
%   details see for example:
%
%  * Approximate sparsity patterns for the inverse of a matrix and
%  preconditioning, T. Huckle, Applied Numerical Mathematics 30 (1999) 291-303
%
%  * A comparative study of sparse approximate inverse preconditioners, M.
%  Benzi, M. Tuma, Applied Numerical Mathematics 30 (1999) 305-340.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function M=sappinv(A,Apriori)
n=length(A);
I=speye(n);
M=sparse(n,n);

for j=1:n %use parfor
    j % comment this line if you do not want to see the progress
    sj=I(:,j);
    
    % find the non-zero rows of a priori pattern
    col_el=find(Apriori(j,:)); 
    A1=A(:,col_el);
    [rows_el,~]=find(A1);
    rows_el=unique(rows_el);
    A1=A1(rows_el,:);
    sj=sj(rows_el,:);
 
    % the QR implementation
    bj=A1\sj;
    % this is to test
    % bj=inv(A1'*A1)*A1'*sj; 
    mj=sparse(n,1);
    mj(col_el)=bj;
    M(:,j)=mj;
end

end