function cohpower = cohPSD(eigvecs_l, eigvals, eigvecs_r, varargin)
    %COHPSD Compute power spectral density from reconstructed cross-spectrum
    %
    % This function takes the outputs of recon_mtx and computes the power spectral
    % density (PSD) by taking the diagonals of the first two dimensions of the
    % reconstructed cross-spectrum, leaving other dimensions intact.
    %
    % Usage:
    %    PSD = cohPSD(eigvecs_l, eigvals, eigvecs_r)
    %    PSD = cohPSD(eigvecs_l, eigvals, eigvecs_r, modes)
    %
    % Inputs:
    %    eigvecs_l - A matrix containing left eigenvectors.
    %    eigvals - A matrix containing eigenvalues.
    %    eigvecs_r - A matrix containing right eigenvectors.
    %    modes - (optional) A scalar specifying which modes to be used in reconstruction.
    %
    % Outputs:
    %    PSD - The power spectral density matrix.
    %
    % Example: 
    %    output = gcoh_plus(data, epochs, mtparams, gcohparams);
    %    Pxy = recon_mtx(output.eigenvectors_l, output.eigenvalues, output.eigenvectors_r);
    %    PSD = cohPSD(output.eigenvectors_l, output.eigenvalues, output.eigenvectors_r);
    %
    % See also: RECON_MTX, GCOH_PLUS

    if ~exist('modes','var') || isempty(modes)
        modes = 1;
    end
    
    % Compute the cross-spectrum
    Pxy = recon_mtx(eigvecs_l, eigvals, eigvecs_r, modes);

    % Initalize cohpower
    sizePxy = size(Pxy);
    cohpower = zeros([sizePxy(1) sizePxy(3:end)]);
    
    % Get the diagonals of the first two dimensions, leaving other dimensions intact
    if length(sizePxy) == 4
        for i = 1:size(Pxy, 3)
            for j = 1:size(Pxy, 4)
                cohpower(:,i,j) = diag(Pxy(:,:,i,j));
            end
        end

    elseif length(sizePxy) == 5
        for i = 1:size(Pxy, 3)
            for j = 1:size(Pxy, 4)
                for k = 1:size(Pxy,5)
                    cohpower(:,i,j) = diag(Pxy(:,:,i,j));
                end
            end
        end

    end

    cohpower = real(cohpower);

end
