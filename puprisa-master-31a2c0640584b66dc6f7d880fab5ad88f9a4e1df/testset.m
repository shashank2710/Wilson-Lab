function mnuBiexpTestSet( src, evt )
    % come up with a 32x32 test set
    
    nrow = 16;
    ncol = 16;
    nt = 64;
    t = linspace(0,10,nt);
    
    % generate biexp. decays for the 4 quadrants
    %[A10,T10,A20,T20,off0]
    A_ul = [.5,3,.5,0.2,0];
    A_ll = [0.5,3,-0.5,0.2,0];
    A_ur = [-2.0,5, 1.0, 0.5, 0];
    A_lr = [-0.1, 0.4, 1.0, 6, 0];
    
    fn = @(A,xx) A(1)*exp(-xx./A(2)) + A(3)*exp(-xx./A(4))+A(5);
    
    f_ul = fn(A_ul, t);
    f_ll = fn(A_ll, t);
    f_ur = fn(A_ur, t);
    f_lr = fn(A_lr, t);
    
    % make an image stack
    X(1:8,1:8,:) = repmat(reshape(f_ul,[1,1,nt]),[8,8,1]);
    X(9:16,1:8,:) = repmat(reshape(f_ll,[1,1,nt]),[8,8,1]);
    X(1:8,9:16,:) = repmat(reshape(f_ur,[1,1,nt]),[8,8,1]);
    X(9:16,9:16,:) = repmat(reshape(f_lr,[1,1,nt]),[8,8,1]);
    
    
    % then store appropriately in appdata, etc.
    setappdata(gcbf,'stackType','delay stack');
    setappdata(gcbf,'nChannels',1);
    setappdata(gcbf,'nSlices',nt);
    header.fullHeaderText = 'Biexponential test set';
    header.pixelsperline = 32;
    header.scanrangex = 1.0;
    setappdata(gcbf,'fileHeader',header);
    setappdata(gcbf, 'delays', t);
    setappdata(gcbf,'currentSlice',1);
    imageStackChannels{1} = X;
    setappdata(gcbf,'imageStackChannels',imageStackChannels);
    
    hChannelImages = getappdata(gcbf,'hChannelImages');
    hChannelAxes = getappdata(gcbf,'hChannelAxes');
    
    setappdata(hChannelAxes(1),'imageStack', X);
    
    set(hChannelImages(1),'XData',[1,ncol],'YData',[1,nrow],...
            'CData',squeeze(X(:,:,1)));
    set(hChannelAxes(1),'XLim',[1,ncol],'YLim',[1,nrow]);

    % display all channels for current slice
    updateAll( gcbf );
end