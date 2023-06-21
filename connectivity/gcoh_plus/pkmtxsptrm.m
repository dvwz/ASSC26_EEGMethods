% PKMTXSPTRM     Calculate cross spectrum using a multitaper method.
%                DPSS tapers from 0 to floor(2*NW)-1 are averaged to
%                generate coherence spectrum.
%
% Usage:
%   [Pxy,f] = pkmtxsptrm(data,params)
%
% Input:
%   data    in a format of samples x channels/trials -- required
%   params  structure with the following fields -- optional
%    .tapers  precalculated tapers from dpss or in the one of the following
%             forms:
%             (1) A numeric vector [TW K] where TW is the
%                 time-bandwidth product and K is the number of
%                 tapers to be used (less than or equal to
%                 2TW-1).
%             (2) A numeric vector [W T p] where W is the
%                 bandwidth, T is the duration of the data and p
%                 is an integer such that 2TW-p tapers are used. In
%                 this form there is no default i.e. to specify
%                 the bandwidth, you have to specify T and p as
%                 well. Note that the units of W and T have to be
%                 consistent: if W is in Hz, T must be in seconds
%                 and vice versa. Note that these units must also
%                 be consistent with the units of params.Fs: W can
%                 be in Hz if and only if params.Fs is in Hz.
%                 The default is to use form 1 with TW = 3 and K = 5
%
%    .pad     padding factor for the FFT -- optional
%             can take integer values -1,0,1,2...
%             -1 corresponds to no padding,
%              0 corresponds to padding to the next highest power of 2 etc.
%             e.g. For N = 500, if PAD = -1, we do not pad; if PAD = 0, we
%             pad the FFT to 512 points, if pad = 1, we pad to 1024 points
%             etc. Defaults to 0.
%    .Fs      sampling frequency -- optional. Default 1.
%    .fpass   frequency band to be used in the calculation -- optional
%             Default all frequencies between 0 and Fs/2
%    .err     error calculation -- optional. Default 0.
%             [1 p] - Theoretical error bars;
%             [2 p] - Jackknife error bars
%             [0 p] or 0 - no error bars
%  .trialave  average over trials/channels -- optional. Default 0
%             average if 1, don't average if 0
% (Note! units have to be consistent. See chronux.m for more information.)
%
% Output:
%   Pxy     cross spectrum (Nfreq x Nchan x Nchan)
%   f       frequency (Nfreq x 1)
%
%********************************************************************

function [Pxy,f] = pkmtxsptrm(data,params)

if nargin < 1; 
    error('Need data'); 
end
if nargin < 2; 
    params = []; 
end

[taps,pad,Fs,fpass] = getparams(params);

NW = taps(1); 
K = taps(2); 

data = change_row_to_column(data);
[N,Nchan] = size(data);

data = data - repmat(mean(data),N,1);
nfft = max(2^(nextpow2(N)+pad),N);

if ~isfield(params,'df')
    df = Fs/nfft;
else
    df = params.df;
end

[f,findx] = getfgrid(Fs,nfft,fpass); f = f(:); findx = findx(:);
%[f,findx] = getfgrid_df(Fs,nfft,fpass,df); f = f(:); findx = findx(:);

tapers = dpsschk([NW,K],N,Fs); % check tapers

J = mtfftc(data,tapers,nfft,Fs);
J = J(findx,:,:); Nf = length(findx);

Jx = conj(repmat(J,[1,1,1,Nchan]));
Jy = repmat(reshape(J,[Nf,K,1,Nchan]),[1,1,Nchan,1]);
Pxy = Jx.*Jy;

% Pyy = mean(abs(J).^2,2);
% Pxx = repmat(squeeze(Pyy),[1,1,Nchan]);
% Pyy = repmat(Pyy,[1,Nchan,1]);
%
% Cxy = (abs(Pxy).^2)./(Pxx.*Pyy);

% [EOF] pkmtxsptrm.m