%FDB_MTX Frequency-domain bootstrap for multitaper cross-spectrum
%Optional file header info (to give more details about the function than one line)
%Optional file header info (to give more details about the function than one line)
%
% Usage:  
%    Pxy_out = fdb_mtx(Pxy, Niter, version)
%
% Inputs:
%    Pxy - Cross-spectral estimates, in the form Ntapers x Nchans x Nchans x Nmw
%    Niter - Number of bootstrap iterations to perform
%    version - 'epoch', median across time, or 'window', does not
%
% Outputs:
%    Pxy_out - Averaged cross-spectral replicates. If version is 'epoch',
%       takes the form Nchans x Nchans x Niter. If version is 'window', takes
%       the form Nchanx x Nchans x Nmw*Niter.
%
% Other m-files required: none
% Other requirements: none
%
% See also: bootstrap_mtxsptrm.m (intracranial)
%
% Copyright Apr-2020, David Zhou, dwzhou@mit.edu
% Last revision 06-Apr-2020
%------------------------------------------------
 
function Pxys_out = fdb_mtx(Pxys, Niter, version)

%---------------------------------------
% PROCESS INPUTS
%---------------------------------------

% size parameters
if ~exist('Niter','var') || isempty(Niter)
    Niter = 1000;
end
Nchans = size(Pxys,2);
Ntapers = size(Pxys,1);

% remove nan times
tpos_nonnan = ~isnan(Pxys(1,1,1,:));
Pxys = Pxys(:,:,:,tpos_nonnan);
Nmw = size(Pxys,4);

% init output
switch version
    case 'window'
        Pxys_out = nan(Nchans,Nchans,Nmw*Niter);
    case 'epoch'
        Pxys_out = nan(Nchans,Nchans,Niter);
end

% random seed
rng(1,'twister');

% populate rands
resample_tapers = randi(Ntapers,Ntapers,Niter);

%---------------------------------------
% COMPILE REPLICATES
%---------------------------------------

switch version
    
    case 'window'
        
        c=1;
        for i = 1:Niter
            for mw = 1:Nmw
            
                replPxy = Pxys(resample_tapers(:,i),:,:,mw); % sample Pxys
                replPxy = squeeze(mean(replPxy,1)); % average over tapers

                Pxys_out(:,:,c) = replPxy;
                c=c+1;
            
            end
        end
    
    case 'epoch'
        
        % populate tpos rands
        resample_tpos = randi(Nmw,Nmw,Niter);

        for i = 1:Niter

            replPxy = Pxys(resample_tapers(:,i),:,:,resample_tpos(:,i)); % sample Pxys
            replPxy = squeeze(mean(replPxy,1)); % average over tapers
            replPxy = squeeze(nanmedian(real(replPxy),3))+...
                1i*squeeze(nanmedian(imag(replPxy),3)); % median over time chunks

            Pxys_out(:,:,i) = replPxy;

        end
        
end

%---------------------------------------
% END CODE
%---------------------------------------
end
