function X = dftapplyshift( X0, diffphase, row_shift, col_shift )
% function [] = dftapplyshift()
% Takes parameters from Manuel Guizar-Sicairos' dftregistration() and
% applies them to shift an image.
%
% input args derived from output of dftregistration()
%       diffphase = output(2);
%       row_shift = output(3);
%       col_shift = output(4);
%
% Based on:
% Manuel Guizar-Sicairos, Samuel T. Thurman, and James R. Fienup, 
% "Efficient subpixel image registration algorithms," Opt. Lett. 33, 
% 156-158 (2008).

    % copied from dft registration code
    [nr,nc]=size(X0);
    Nr = ifftshift([-fix(nr/2):ceil(nr/2)-1]);
    Nc = ifftshift([-fix(nc/2):ceil(nc/2)-1]);
    [Nc,Nr] = meshgrid(Nc,Nr);
    Greg = fft2(X0).*exp(i*2*pi*(-row_shift*Nr/nr-col_shift*Nc/nc));
    Greg = Greg*exp(i*diffphase);
    X = real(ifft2(Greg));
    
end