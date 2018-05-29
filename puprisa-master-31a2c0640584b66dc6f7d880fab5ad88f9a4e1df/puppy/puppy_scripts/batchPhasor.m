function out = batchPhasor(S, header)
    [t,X] = puprisa_getChannelFromSlices(S, 1, header);
    
    [X, ~] = puprisa_baselineCorrection( t, X, 6 );
    
    [nr,nc,nt] = size(X);
    
    I = sum(X.^2,3) / nt;
    %[SNR, SNRdb, noise] = ImSNRest(X(:,:,1));
    I_frame1 = X(:,:,1).^2;
    noise_rms = sqrt(std(I_frame1(:)))    % root mean-squared noise
    
    thresh = 2*noise_rms / sqrt(nt);
    
    mask = logical(I >= thresh);
    
%     figure; imagesc(I.*mask); colorbar;
%     cm=colormap(jet);
%     cm(1,:) = [0,0,0];
%     colormap(cm);
%     drawnow;
%     %title(['SNR: ', num2str(SNR)]);
%     %pause(1);
%     close;
    
    H = phasor(X,t,'batch','mask',mask);
    
    Hmean = mean(H(:));
    Hstd = std(H(:));
    
    % make an image
    Himg = H / (Hmean + 4*Hstd);
    Himg( Himg > 1 ) = 1;
    Himg = 255* Himg;
    Himg = uint8(Himg);
    
    out = [];
    out = setfield(out, 'phasorPlot', H);
    out = setfield(out, 'phasorPlotImg', Himg);
    %out = setfield(out, names{2}, F(2));
    
    %out.euOverPheo = F(2) / F(1);
    
    %out
