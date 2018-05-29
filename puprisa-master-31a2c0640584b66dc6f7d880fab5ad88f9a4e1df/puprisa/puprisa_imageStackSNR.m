function SNR = puprisa_imageStackSNR( imageStack )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to estimate the SNR of an delay stack
% Input(s):
%     totalImage - Gray scale image
%
% Output(s):
%     SNRtotEst - Estimate of the image SNR
% 
% Author: Jesse W. Wilson (2012)
% Based on FERobles' implementation of JTL Thong et al. "Single-Image 
% Signal to Noise Ratio Estimation," Scanning Vol. 23 (5), 2001
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    y = squeeze(imageStack(154,846,:));
    figure;
    plot(y);
    
    SNR = pixelSNR(y);
        
end

function SNRtotEst = pixelSNR( y )
    %[ACF]=autocorr(totalImage(:),length(totalImage(:))-1);
    y = y - mean(y);
    ACF = xcorr( y );

    %reshACF=reshape(ACF,nx,ny);

    %reshACF=[reshACF(end/2+1:end,1:end/2);reshACF(1:end/2,1:end/2)];

    
    r0 = max(reshACF(:)); 
    zerodelay = find(reshACF == r0);
    mu = min(reshACF(:));
    rhat = (reshACF(zerodelay,2) + reshACF(zerodelay-1,1)/2 + reshACF(zerodelay+1,1)/2)/2;
    SNRtotEst = 10*log((rhat - mu^2)/(r0-rhat));
end