function angcube = angdist_boot(eigvecs_a,eigvecs_b,modes)
% ANGDIST_BOOt Computes pairwise angular distances between eigenvectors for bootstrap iterations.
%
% Usage:
%    angcube = angdist_boot(eigvecs_a, eigvecs_b)
%    angcube = angdist_boot(eigvecs_a, eigvecs_b, modes)
%
% Inputs:
%    eigvecs_a, eigevecs_b 
%            - A 4-D matrix with eigenvectors as columns in the first two dimensions.
%              The third dimension represents different measurements and the fourth dimension
%              represents bootstrap iterations (Nchans x Nchans x Nmw x Niter).
%    modes   - (Optional) Specifies the modes (columns) to use for computation.
%              Defaults to the first column if not provided.
%
% Outputs:
%    angcube - A cube of angular distance matrices, where each matrix corresponds to an
%              iteration of the bootstrap.
%
% This function computes the pairwise angular distances between the columns of the 'eigvecs'
% matrix for each bootstrap iteration. The angular distance is computed as the principal angle 
% between two subspaces spanned by the eigenvectors, using the MATLAB 'subspace' function.
%
% If 'modes' is not provided, the function uses the first column of 'eigvecs' by default.
%
% Example:
%    eigenvectors = rand(5, 5, 10, 100);  % Generate some random eigenvectors
%    angcube = angdist_boot(eigenvectors, eigenvectors, [1, 2]);  % Compute angular distances for the first two modes
%
% See also: SUBSPACE

if ~exist('modes','var') || isempty(modes)
    modes = 1;
end

assert(length(size(eigvecs_a)) == length(size(eigvecs_b)))

n_eigvecs_a = size(eigvecs_a,3);
n_eigvecs_b = size(eigvecs_b,3);
n_iter = size(eigvecs_a,4); % Number of bootstrap iterations

eigvecs_a = eigvecs_a(:,modes,:,:);
eigvecs_b = eigvecs_b(:,modes,:,:);

angcube = zeros(n_eigvecs_a, n_eigvecs_b, n_iter);

for k = 1:n_iter
    for i = 1:n_eigvecs_a
        for j = 1:n_eigvecs_b
            if i < j
                angcube(i,j,k) = subspace(eigvecs_a(:,:,i,k),eigvecs_b(:,:,j,k));
            end
        end
    end
end

angcube = angcube + permute(angcube, [2 1 3]);

end
