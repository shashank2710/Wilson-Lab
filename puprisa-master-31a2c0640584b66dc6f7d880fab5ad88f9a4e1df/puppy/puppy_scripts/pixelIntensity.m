function out = pixelIntensity(S, header)

    [t,X] = puprisa_getChannelFromSlices(S, 1, header);

    [X, ~] = puprisa_baselineCorrection( t, X, 6 );
    
    I = (sum((X),3));
    Imean = mean(I(:));
    Istd = std(I(:));
    
    mask = uint8(abs(I) >= 3);
    
    Iimage = I / (6*Istd);
    Iimage( Iimage > 1 ) = 1;
    Iimage( Iimage < -1 ) = -1;
    Iimage = (Iimage + 1) / 2;;
    Iimage = 255*Iimage;
    Iimage = uint8(Iimage);
    
    Ithresh = Iimage.*mask + 127*(1-mask);
    
    out.Iimage = Iimage;
    out.I = I;
    out.I_vect = I(:);
    out.Ithresh = Ithresh;  