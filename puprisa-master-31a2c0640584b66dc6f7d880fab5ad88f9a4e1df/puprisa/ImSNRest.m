function [SNR, SNRdb, noise] = ImSNRest(totalImage)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to estimate the SNR of an image
% Input(s):
%     totalImage - Gray scale image
%
% Output(s):
%     SNRtotEst - Estimate of the image SNR
% 
% Author: FERobles 2012
% Based on JTL Thong et al. "Single-Image Signal to Noise Ration
% Estimation," Scanning Vol. 23 (5), 2001
%
% modified 10/2012 Jesse W. Wilson
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[nx, ny] = size(totalImage);
[ACF]=autocorr(totalImage(:),length(totalImage(:))-1);
reshACF=reshape(ACF,nx,ny);

reshACF=[reshACF(end/2+1:end,1:end/2);reshACF(1:end/2,1:end/2)];

r0 = max(reshACF(:)); 
zerodelay = find(reshACF == r0);
mu = min(reshACF(:));
rhat = (reshACF(zerodelay,2) + reshACF(zerodelay-1,1)/2 + reshACF(zerodelay+1,1)/2)/2;
noise = r0-rhat;
SNR = (rhat - mu^2)/(r0-rhat);
SNRdb = 10*log10(SNR);
end