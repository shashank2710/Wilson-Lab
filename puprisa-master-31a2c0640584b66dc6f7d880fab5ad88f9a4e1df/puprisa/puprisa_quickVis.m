% puprisa_quickVis
% quick visualization of pump-probe data

function IRGB = puprisa_quickVis( imageStack, t )

    % get dimensions of the image
    [nr, nc, nt] = size( imageStack );

    % instrument response function
    tfwhm = 0.3; 
    
    % generate the basis
    Bshort = (t >= -tfwhm) & (t <= tfwhm);
    
    Blong = t > tfwhm;
    
    Boffs = t*0 + 1;
    
    
    % unwrap the image into a list of delay scans
    ds = reshape( imageStack, [nr*nc, nt] );
    
    B = [Bshort; Blong; Boffs];
    size(B)
    
    % find the mixing matrix (least sq fit of the image stack to our
    % basis)
    A = ds * pinv( B );
    size(A)
    
    % generate images from each component
    I_short =  reshape( A(:,1), [nr, nc]);
    I_shortPos = I_short.*( I_short > 0 );
    I_shortNeg = -I_short.*( I_short < 0 );
    
    I_long = reshape( A(:,2), [nr, nc]);
    I_longPos = I_long.*( I_long > 0 );
    I_longNeg = -I_long.*( I_long < 0 );
    
    % color it!
    
    % (old RGB values)
    %c_shortPos = [1, 0.5, 0];
    %c_shortNeg = [0, 0, 1];
    %c_longPos = [1, 0, 0.5];
    %c_longNeg = [0, 1, 0];
    
    % (new HSV colors)
    c_shortPos = hsv2rgb([45/360, 1, 1]);
    c_shortNeg = hsv2rgb([225/360, 1, 1]);
    c_longPos = hsv2rgb([315/360, 1, 1]);
    c_longNeg = hsv2rgb([135/360, 1, 1]);
    
    IRGB_shortPos = colorize( I_shortPos, c_shortPos );
    IRGB_longPos = colorize( I_longPos, c_longPos );
    
    IRGB_shortNeg = colorize( I_shortNeg, c_shortNeg );
    IRGB_longNeg = colorize( I_longNeg, c_longNeg );
    
    IRGB = IRGB_longPos + IRGB_longNeg + IRGB_shortPos + IRGB_shortNeg;
    
%     sc = mean(IRGB(:)) + 5*std(IRGB(:));
%     IRGB2 = IRGB / sc;
%     IRGB2( IRGB2 > 1 ) = 1;
%     
%     figure;
%     imshow(IRGB2);

end

function RGB = colorize( I, c )
    RGB = repmat(I, [1,1,3] );
    RGB(:,:,1) = RGB(:,:,1) * c(1);
    RGB(:,:,2) = RGB(:,:,2) * c(2);
    RGB(:,:,3) = RGB(:,:,3) * c(3);
end