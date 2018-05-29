function imfiltd = ifilt_highPass( im )
    f = fftshift(fft2(im));
    [szx, szy] = size(im);
    filt1d = 1-hamming( szx );
    
    filt2d = repmat( filt1d, [1, szy] ).*repmat( filt1d, [1, szy] ).';
    
    f = f .* filt2d;
    
    imfiltd = real( ifft2( ifftshift( f )));
end

