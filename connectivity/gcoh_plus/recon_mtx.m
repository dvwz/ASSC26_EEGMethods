%RECON_MTX Reconstruct cross-spectrum from selected eigenvalues and eigenvectors
%
% This function is designed to reconstruct the cross-spectrum of a given 
% data set based on selected eigenvalues and eigenvectors. It includes an optional
% 'modes' parameter, allowing users to specify which modes should be included in the 
% reconstruction. In case of not specifying the 'modes', it defaults to 1. It also
% handles input matrices of different dimensions.
%
% Usage:  
%    Pxy_recon = recon_mtx(eigvecs_l, eigvals, eigvecs_r)
%    Pxy_recon = recon_mtx(eigvecs_l, eigvals, eigvecs_r, modes)
%
% Inputs:
%    eigvecs_l - A matrix containing left eigenvectors.
%    eigvals - A matrix containing eigenvalues.
%    eigvecs_r - A matrix containing right eigenvectors.
%    modes - (optional) A scalar specifying which modes to be used in reconstruction.
%
% Outputs:
%    Pxy_recon - The reconstructed cross-spectrum matrix.
%
% Example: 
%    output = gcoh_plus(data, epochs, mtparams, gcohparams);
%    Pxy = recon_mtx(output.eigenvectors_l, output.eigenvalues, output.eigenvectors_r);
%    Pxy_selected_modes = recon_mtx(output.eigenvectors_l, output.eigenvalues, output.eigenvectors_r, [1:3]);
%
% Other m-files required: none
% Other requirements: none
%
% See also: GCOH_PLUS
%
% Copyright Apr-2020, David Zhou, dwzhou@mit.edu
% Last revision 06-Apr-2020
%------------------------------------------------

function Pxy_recon = recon_mtx(eigvecs_l,eigvals,eigvecs_r,modes)

%---------------------------------------
% PROCESS INPUTS
%---------------------------------------

if ~exist('modes','var') || isempty(modes)
    modes = 1;
end
 
%VARARGIN
switch length(size(eigvals))
    case 2
        disp('gcohout dims not implemented')
        return
    case 3 % Nch x Nfreqs x (Niter/Nmw)
        Nch = size(eigvals,1);
        Nfreqs = size(eigvals,2);
        
        Pxy_recon = nan(size(eigvecs_l));
    case 4 % Nch x Nfreqs x Nmw x Niter
        disp('gcohout dims not implemented')
        Pxy_recon = nan(size(eigvecs_l));

end

%---------------------------------------
% MAIN BODY
%---------------------------------------

switch length(size(eigvals))
    case 2
        disp('gcohout dims not implemented')
        return
        
    case 3
        for i = 1:size(eigvals,2)
            for j = 1:size(eigvals,3)
                
                Pxy_recon(:,:,i,j) = reconstruct_Pxy(eigvecs_l(:,:,i,j),...
                                                     eigvals(:,i,j),...
                                                     eigvecs_r(:,:,i,j),...
                                                     modes);
                
            end
        end
        
    case 4
        for i = 1:size(eigvals,2)
            for j = 1:size(eigvals,3)
                for k = 1:size(eigvals,4)
                
                    Pxy_recon(:,:,i,j) = reconstruct_Pxy(eigvecs_l(:,:,i,j,k),...
                                                         eigvals(:,i,j,k),...
                                                         eigvecs_r(:,:,i,j,k),...
                                                         modes);
                
                end
            end
        end
        
end

%---------------------------------------
% SUBFUNCTIONS
%---------------------------------------

    function Pxy_recon = reconstruct_Pxy(U,S,V,modes)
        
        S_new = zeros(size(S));
        S_new(modes) = S(modes);
        S_new = eye(length(S)) .* S_new;
        Pxy_recon = U * S_new * V';
        
    end

%---------------------------------------
% END CODE
%---------------------------------------
end
