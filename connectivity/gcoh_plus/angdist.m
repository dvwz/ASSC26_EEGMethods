function angmat = angdist(eigvecs_a,eigvecs_b,modes)
% ANGDIST Computes pairwise angular distances between eigenvectors.
%
% Usage:
%    angmat = angdist(eigvecs_a, eigvecs_b)
%    angmat = angdist(eigvecs_a, eigvecs_b, modes)
%
% Inputs:
%    eigvecs_a, eigevecs_b 
%            - A 3-D matrix with eigenvectors as columns in the first two dimensions.
%              The third dimension represents different measurements (Nchans x Nchans x Nmw...).
%    modes   - (Optional) Specifies the modes (columns) to use for computation.
%              Defaults to the first column if not provided.
%
% Outputs:
%    angmat  - A symmetric matrix of angular distances, where the entry (i, j) is the angular 
%              distance between the i-th and j-th eigenvectors.
%
% This function computes the pairwise angular distances between the columns of the 'eigvecs'
% matrix. The angular distance is computed as the principal angle between two subspaces spanned
% by the eigenvectors, using the MATLAB 'subspace' function.
%
% If 'modes' is not provided, the function uses the first column of 'eigvecs' by default.
%
% Example:
%    eigenvectors = gcohout.eigenvectors_l;  % Generate some random eigenvectors
%    angMatrix = angdist(eigenvectors, [1, 2]);  % Compute angular distances for the first two modes
%
% See also: SUBSPACE

if ~exist('modes','var') || isempty(modes)
    modes = 1;
end

assert(length(size(eigvecs_a)) == length(size(eigvecs_b)))

n_eigvecs_a = size(eigvecs_a,3);
n_eigvecs_b = size(eigvecs_b,3);

eigvecs_a = eigvecs_a(:,modes,:);
eigvecs_b = eigvecs_b(:,modes,:);

angmat = zeros(n_eigvecs_a, n_eigvecs_b);

for i = 1:n_eigvecs_a
    for j = 1:n_eigvecs_b
        
        if i < j
            angmat(i,j) = subspace(eigvecs_a(:,:,i),eigvecs_b(:,:,j));
        end
    
    end
end

angmat = angmat' + angmat;

end

