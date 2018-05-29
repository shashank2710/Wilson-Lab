% puprisa_viewChannelZStack
% Jesse Wilson (2011,2012) syrex314@gmail.com
% Duke University
%
% View a single channel as a z stack.
%
% A visualization tool for 3D stacks of images, such as those acquired in
% pump--probe microscopy.
% Shows one z slice at a time, and also allows user to change
% the currently viewed slice
%
% The purpose is to have a tool for visualizing 3D image stacks.
%
% X is an n by m by l - element matrix where
% m is the number of pixels in the y-direction
% n is the number of pixels in the x-direction
% p is the number of time delays (the number of image slices)
%
%
% SOMEDAY: implement callbacks for cbar autoscale context menu
% Then clean up this code before going any further...
%
% instead of using MATLAB slice:
% just make an image plane, change CDATA, and move it up and down
%
function puprisa_viewChannelZStack( X, zPos, fileName, header, displayMode )
    figure('Name','Channel');
    if nargin == 4
        displayMode = '2D';
    end
    setappdata(gcf,'displayMode',displayMode);
    
    % maximize the window vertically
    pos = get(gcf,'OuterPosition');
    screenSize = get(0,'ScreenSize');
    pos(2) = screenSize(2)+100;
    pos(4) = screenSize(4)-100;
    set(gcf,'OuterPosition',pos);
    
    % demo mode
    if nargin == 0
        X = demoDataSet();
        Y = demoDataSet();
        [szx, szy, szz] = size(Y);
        zPos = 0.1*(1:szz);
    elseif isstr( X )
        fileName = X;
        [X,Y,zPos] = loadData( fileName ); 
    end
    
    % get dimensions of the stack
    % assume X and Y channels have same dimensions!
    [m,n,p] = size(X);
       
    clf;
    
    initMenus( gcf );
    
    % set up 2 subplots
    % 1) image for selected time delay
    % 2) depth profile
    
    ik = 1; % ik is slice
    setappdata(gcf,'currentSlice',ik);

    hImageAxes = gca();
    set(hImageAxes,'Position',[0.02,0.33,0.96,0.7])

    displayMode = getappdata(gcf,'displayMode');
    setappdata(gcf,'imageStackX', X);
    
    switch displayMode
        case '3D'
            % use Kevin's 3D viewing code
            drawSurfaces();
        case '2D'
            im = imagesc(X(:,:,ik));
            set(im,'Tag','theImage','ButtonDownFcn',@axesButtonDownFcn);
            axis square;
    end
    
    puprisa_colorMap('pumpProbe');
    colorbar;

    set(gca,'YDir','normal');
    set(gca,'Tag','imageAxes');
    setappdata(gcf,'zPos',zPos);
    
    puprisa_initScaleBar( gca, header );
    
    if nargin >= 3
        title(fileName,'interp','none');
        setappdata(gcf,'fileName',fileName);
    else
        title('viz3dStack');
    end

    % set up scroll wheel fn to be able to flip between slices
    set(gcf,'WindowScrollWheelFcn', @scroll_wheel);
    
    % set up keypress handler for flipping between slices
    set(gcf,'KeyPressFcn', @keyPress);
    
    % attempt auto-clim
    % take the biggest difference between mean and zero
    m = mean(X(:));
    s = 3*std(X(:));
    set(hImageAxes,'clim', [-1,1]*max( [abs(m+s),abs(m-s)] ));
    
    if(strcmp(displayMode,'2D'))
        % display depth
        hTxtDepth = text(0,header.linesperframe,depthString(zPos(ik)),...
            'color',[0,0,0]+0.99,...
            'HorizontalAlignment','left',...
            'VerticalAlignment','top',...
            'FontSize',24,...
            'tag','txtDepth');

        % plot depth trend
        axTrend = axes('Position',[0.05,0.05,0.90,0.3],'tag','axZtrend');
        trend = puprisa_imageStackTrend(X);
        lZTrend = line(zPos,trend,'Tag','lZTrend');
        
        setappdata(gcf,'ZlineAMP',[min(trend),max(trend)]);
        Zlocation =line([zPos(ik),zPos(ik)],[min(trend),max(trend)]...
            ,'Color','g','LineStyle','--','Tag','Zlocation');
        
        xlabel('Z Position (\mum)')
        setappdata(gcf,'trend',trend);
        grid on;
    end
    
    % disable contrast stretching by default unless the user requests it
    setappdata( gcf, 'stretchContrast', 0 );
    
    % this setting used by auto-contrast methods
    setappdata( gcf,'imagePolarity','bipolar' );

    
end

function drawSurfaces()
    X = getappdata(gcf,'imageStackX');
    
    [nx,ny,nz] = size(X);
    
    hSurfs = zeros(1,nz);
    clf;
    axes();
    set(gca,'color','k');
    for iz = 1:nz
        
        hSurfs(iz) = ...
            surface(iz + zeros(nx,ny),X(:,:,iz),...
            'AlphaData',abs(X(:,:,iz)),...
            'FaceColor',[0 1 0],...
            'EdgeColor','none',...
            'CDataMapping','scaled',...
            'FaceAlpha','flat');
    end
    alim([0,0.05]);
    view(-35,15);
    
    setappdata(gcf,'hSurfs',hSurfs);
end

function str = depthString(z)
    str = ['z = ',num2str(z), ' \mum'];
end

%--------------------------------------------------------------------------
% INITIALIZATION ROUTINES
%
% menu initialization
function initMenus( fig )
    % sets up menus in the figure
    f = uimenu('Label','ImageAnalysis');
    uimenu(f,'Label','Open Image Stack...','Callback',@mnuOpenImageStack);
    uimenu(f,'Label','Histogram...','Callback',@mnuShowHistogram);
    uimenu(f,'Label','Multiexponential fit...','Callback',@mnuMultiexpFit);
    uimenu(f,'Label','Principal component composite image...','Callback',@mnuPCA);
    uimenu(f,'Label','Principal component cluster analysis...','Callback',@mnuPCAClusters);
    uimenu(f,'Label','K-clustering analysis...','Callback',@mnuPCAClusters2);
    uimenu(f,'Label','Print stack on page...','Callback',@mnuPrintStack);
    uimenu(f,'Label','Subract background value...','Callback',@mnuBgSubtract);
    uimenu(f,'Label','Auto background','Callback',@mnuBgSubtractAuto);
    uimenu(f,'Label','Export TIFF stack...','Callback',@mnuExportStack);
    uimenu(f,'Label','Export animated GIF','callback',@mnuExportGIF);
    uimenu(f,'Label','Montage','Callback',@mnuExportMontage);
    uimenu(f,'Label','Depth offset...','callback',@mnuDepthOffset);
    uimenu(f,'Label','Subtract average frame','callback',@mnuSubtractAverage);
    uimenu(f,'Label','Average stack','Callback',@mnuAverageStack);
    uimenu(f,'Label','High pass filter','Callback',@mnuHighPassFilter);
    uimenu(f,'Label','Normalize stack','Callback',@mnuNormalizeStack);
    uimenu(f,'Label','Invert stack','Callback',@mnuInvertStack);
    uimenu(f,'Label','Apply exponential corection','Callback',@mnuExponentialCorrection);
    uimenu(f,'Label','Max Projections','Callback',@mnuMaxProjections);
    
    f = uimenu('Label','Display');
    % this menu contains checkboxes on what to display
    uimenu(f,'Label','Depth','Callback',@mnuDisplayDepth,'checked','on');
    uimenu(f,'Label','3D','Callback',@mnuDisplay3D);
      
    % set up mode menu for selecting colormaps
    puprisa_menu_mode( );
end


%==========================================================================
% CALLBACK FUNCTIONS
%

function scroll_wheel(src,evnt )
    ik = getappdata(gcbf,'currentSlice') + evnt.VerticalScrollCount;
    
    updateDisplay(ik)
end


function keyPress( src, evnt )
    % arrow keys?
    incr = 0;
    
    if( evnt.Character == char(28) )
        incr = -1;
    elseif( evnt.Character == char(29) )
        incr = +1;
    end
    
    if incr ~= 0
        ik = getappdata(gcbf,'currentSlice') + incr;
        updateDisplay(ik);
    end
end
   
function updateDisplay( ik )
    X = getappdata(gcbf,'imageStackX');
    [m,n,p] = size(X);
    
    Y = getappdata(gcbf,'imageStackY');
    
    if nargin == 0
        ik = getappdata(gcbf,'currentSlice');
    else
        if ik < 1
            ik = 1;
        elseif ik > p
            ik = p;
        end

        setappdata(gcbf,'currentSlice', ik);
    end
        
    h = findobj(gcf,'Tag','imageAxes');
    
    if gca ~= h
        axes(h);
    end
    
    fn = getappdata(gcbf,'fileName');
    title([fn, num2str(ik)],'interp','none');
    zPos = getappdata(gcbf,'zPos');
    %[x,y,z] = meshgrid([1:256],[1:256],zPos);
    %h=slice(x,y,z,X, 256, 256, zPos(ik));
    
    cdata = X(:,:,ik);
    
    switch getappdata(gcf,'displayMode')
        case '3D'

            hsliceobj = getappdata(gcf,'hsliceobj');
            cdata = X(:,:,ik);
            zdata = cdata*0 + zPos(ik);

            set(hsliceobj,'CData',cdata, 'ZData', zdata);
        case '2D'
            im =findobj(gcbf,'Tag','theImage');
            set(im,'cdata',cdata);
            
            % optionally stretch contrast
            if( getappdata( gcf, 'stretchContrast' ) )
                % find caxis limits for high-contrast display of the image.
                mn = mean( cdata(:) );
                st = std( cdata(:) );
                N = 4;
                
                switch( getappdata( gcf, 'imagePolarity' ) )
                    case 'positive'
                        caxis( [0,1*N*st + mn] );
                    case 'negative'
                        caxis( [-1*N*st + mn, max(cdata(:))] );
                    case 'bipolar'
                        caxis( [-1,1]*N*st + mn );
                    otherwise
                        error('imagePolarity must be positive, negative, or bipolar');
                end
                
                
                
                % an alternate methd is to use stretchlim() 
                %
                % Because stretchlim expects a positive image, we pass to 
                % it an image with an offset
                %minval = min(cdata(:));
                %clim = stretchlim( cdata - minval, 0.0001 );
                %clim = clim + minval;
                %caxis( clim );
                
            end
    end
    
    hTxtDepth = findobj(gcbf,'tag','txtDepth');
    set(hTxtDepth,'string',depthString(zPos(ik)));
    
    Zlineamp=getappdata(gcf,'ZlineAMP');
    zlinelocation = findobj(gcbf,'Tag','Zlocation');
    set(zlinelocation,'xdata',[zPos(ik),zPos(ik)],'ydata',Zlineamp)
    
    drawnow;
end

function button_up(src,evnt)
    set(gcf,'WindowButtonMotionFcn','');
end

%==========================================================================
% MENU CALLBACKS
%

% Main (top) menu callbaks
function mnuAutoScaleStack(src, evt)
    get(get(get(src,'Parent'),'Parent'))
end

function mnuDepthOffset(src, evt)
    % subtract arbitrary offset from depth scans
    
    % get z position vector and current slice
    zPos = getappdata(gcbf,'zPos');
    currentSlice = getappdata(gcbf,'currentSlice');
    currentZ = zPos(currentSlice);
    
    % ask for offset
    prompt = {'Z Offset:'};
    dlg_title = 'Z Offset (microns)';
    num_lines = 1;
    def = {num2str(currentZ)}; % default is current depth
    answer = inputdlg(prompt,dlg_title,num_lines,def);
    
    if ~isempty(answer)

        z0 = str2num(answer{1});
        if isempty(z0)
            error('invalid entry');
        end

        % apply the offst
        zPos = zPos - z0;
        setappdata(gcbf, 'zPos', zPos);

        % update depth gauge readout
        ik = getappdata(gcbf,'currentSlice');
        hTxtDepth = findobj(gcbf,'Tag','txtDepth');
        set(hTxtDepth,'string',depthString(zPos(ik)));

        % update x-axis of trend plot
        lZTrend = findobj(gcbf,'Tag','lZTrend');
        set(lZTrend,'XData',zPos);

        drawnow;
    end
end

function mnuSubtractAverage( src, evt )
    % average stack, then subtract the average from each frame
    X = getappdata( gcbf, 'imageStackX' );
    
    [~,~,nSlices] = size(X);
    
    avg = sum(X,3) / nSlices;
    
    X = X - repmat(avg,[1,1,nSlices]);
    
    setappdata( gcbf,'imageStackX',X);
    
    updateDisplay();

end

function mnuExponentialCorrection( obj, evt )
    % apply correction factor for exponential signal decay
    
    % ask for correction
    prompt = {'Exponential factor:'};
    dlg_title = 'Apply exponential correction:';
    num_lines = 1;
    def = {'0.02'}; % default is current depth
    answer = inputdlg(prompt,dlg_title,num_lines,def);
    
    if ~isempty(answer)

        a = str2num(answer{1})
        if isempty(a)
            error('invalid entry');
        end
    
    
        X = getappdata( gcbf, 'imageStackX' );
        X_original = getappdata(gcbf,'imageStackX_original');
        if isempty( X_original )
            X_original = X;
            setappdata(gcbf,'imageStackX_original',X_original);

        end
        
        X = X_original;
        [nx,ny,nz] = size(X);

        z = getappdata(gcbf,'zPos');
    
        exp1D = exp(-z * a);
        exp3D = repmat( shiftdim(exp1D,-1), [nx, ny, 1] );
        X = X.* exp3D;
        

        setappdata( gcbf,'imageStackX',X);

        updateDisplay();
    end

end


function mnuAverageStack(obj, evt )
    hAxes = findobj(gcbf,'tag','imageAxes');
    cm = colormap();
    cl = get(hAxes,'clim');
    X = getappdata(gcbf,'imageStackX');
    [~,~,nSlices] = size(X);
    Xs = sum(X,3) / nSlices;
    figure('Name','averaged');
    imagesc(Xs);
    axis image;
    set(gca,'ydir','normal');
    colormap(cm);
    caxis(cl);
end

function mnuHighPassFilter( src, evt )
    % transform to Fourier domain, and high-pass filter the entire image
    % stack
    % in order to really be useful, this should have an adjustable cut-off
    
    % iterate through image stack
    X = getappdata(gcbf,'imageStackX');
    [~,~,nSlices] = size(X);
    
    k = [-1,-1,-1;-1,8,-1;-1,-1,-1]; %high-pass convolution kernel
    
    for iSlice = 1:nSlices
        im = X(:,:,iSlice);
        im = conv2(im, k, 'same');
        X(:,:,iSlice) = im;
    end

    setappdata(gcbf,'imageStackX',X);
    updateDisplay();
end

function mnuNormalizeStack( src, evt )
    % normalize each slice to the same intensity level
    
    X = getappdata(gcbf,'imageStackX');
    [nr,nc,nz] = size(X);
    
    I = abs(squeeze(sum(sum(X,1),2)) / (nr*nc));
    
    for iz = 1:nz
        X(:,:,iz) = X(:,:,iz)/I(iz);
    end
    
    setappdata(gcbf,'imageStackX',X);
    updateDisplay();
end

function mnuInvertStack( src, evt )
    % invert stack (for PMT)
    
    X = getappdata(gcbf,'imageStackX');
    setappdata(gcbf,'imageStackX',-1.0*X);
    updateDisplay();
end

    

function mnuPCAClusters2( src, evt )
    % get the image stack
    X = getappdata( gcbf, 'imageStackX' );    
    
    % reshape to a list of delay scans
    [nrows,ncols,ndelays] = size(X);
    delayScans = zeros(nrows*ncols, ndelays);
    for irow = 1:nrows
        minIndex = (irow-1)*ncols+1;
        maxIndex = (irow)*ncols;
        delayScans( minIndex:maxIndex, : ) = squeeze(X( irow, 1:ncols, : ));
    end
    
    % pick the strongest n signal points
    nStrongestSignals = 5000;
    strengths = sum(abs(delayScans),2);
    [b,ii] = sort(strengths,'Descend');

    delayScans = delayScans( ii(1:nStrongestSignals), : );
    delayScanSums = repmat(sum(abs(delayScans),2),[1,ndelays]);
    %delayScans = delayScans ./ delayScanSums;
    

    
    [pc, zscores, pcvars] = princomp(delayScans);
    figure('name','princomp');
    scatter(zscores(:,1), zscores(:,2));
    
    %pcclusters = clusterdata(delayScans,2);
    pcclusters = kmeans(delayScans,2);
    gscatter(zscores(:,1), zscores(:,2), pcclusters);
    
    figure
    foo = zeros(1,nrows*ncols);
    foo(ii(pcclusters == 3)) = 1;
    img = reshape(foo,[nrows,ncols]);
    imagesc(img);
end 

function mnuPCAClusters(src, evt)
    
    % perform principle component analysis,
    
    % get the image stack
    X = getappdata( gcbf, 'imageStackX' );    
    
    % reshape to a list of delay scans
    [nrows,ncols,ndelays] = size(X);
    delayScans = zeros(nrows*ncols, ndelays);
    for irow = 1:nrows
        minIndex = (irow-1)*ncols+1;
        maxIndex = (irow)*ncols;
        delayScans( minIndex:maxIndex, : ) = squeeze(X( irow, 1:ncols, : ));
    end
    
    puprisa_clusterAnalysis( delayScans );
    
    
end

function mnuPCA( src, evt )
    % open dialog to prompt user where to get PCs from,
    % and then show image, with layers of selectable colors for PCs
    % autoscale should have separate scales for each PC
    % perhaps give multiple views, allowing selection of PCs in each view
    
    X = getappdata(gcbf,'imageStackX');
    puprisa_componentImage( X );
    
end

function mnuOpenImageStack( src, evt )
    % prompt user for an image stack to load
    
    % recall working directory
    wd = getappdata(gcbf,'workingDirectory');
    
    % or use MATLAB's working directory if none has been set yet
    if isempty(wd)
        wd = pwd;
    end
    
    [fileName, pathName] = uigetfile([wd,'\*.dat'],...
        'Select a ScanImage Stack' );
    
    if ~isequal(fileName,0)
        % if not canceled, load the file
        [X,Y, zPos] = loadData([pathName,'\',fileName]);
        
        % remember the working directory for later
        setappdata(gcbf,'workingDirectory',pathName);
        
        if ~isempty(X)
            % restart this program with the new data
            viz3dStack( X, Y, zPos, fileName ) ;
        end
    end
end

function mnuShowHistogram( src, evt)
    % show histogram of image stack
    % display distribution of pixel intensity values
    X = getappdata( gcbf, 'imageStackX' );
    figure('Name','histogram');
    hist(X(:),100);
    xlabel('Value');
    ylabel('Pixel count');
end

function mnuMultiexpFit(src, evt)
    % perform a multiexponential fit of the currently-displayed time series
    
    % open dialog to prompt for:
    % - min/max time
    % - no. exponentials
    % - initial guess (implement this later?)
    prompt = {'Min probe delay (ps):','Max probe delay (ps):', 'No. exponentials'};
    dlg_title = 'Parameters for multi-exponential fit';
    num_lines = 1;
    def = {'0.3', 'Inf', '2'};
    answer = inputdlg(prompt,dlg_title,num_lines,def);
    
    % parse user input
    minDelay = str2num( answer(1) );
    maxDelay = str2num( answer(2) );
    nExp     = str2num( answer(3) );
    
    % check for valid inputs
    if( isempty(minDelay) || isempty(maxDelay) || isempty(nExp) )
        error('Invalid inputs to multiexponential fit');
    end
    
    % !!! BOOKMARK
    % call the exponential fitting routine
    
    % then plot the results
    
end

function mnuPrintStack(src,evt)
    
    % get the data
    X = getappdata(gcbf,'imageStackX');
    axes(findobj(gcbf,'Tag','imageAxes'));
    climits = get(gca,'clim');
    cm = colormap();
    
    
    [~,~,nSlices] = size(X);

    % first figure out row/column break-down
    nPanels = nSlices;
    ff = factor(nPanels);
    if length(ff) == 1
        nPanels = nPanels + 1;
        ff = factor(nPanels);
    end
    
    nCols = ff(end);
    nRows = nPanels / ff(end);
    
    figure('Name','Stack layout');
    for iSlice = 1:nSlices
        subplot(nRows, nCols, iSlice);
        imagesc(X(:,:,iSlice));
        colormap(cm);
        set(gca,'clim',climits);
        grid on;
    end
    
    
    % eventually replace subplot() with my own pixel-based positioning
    
end

function mnuBgSubtract( src, evt )
    % ask for constant background level
    prompt = {'Background level:'};
    dlg_title = 'Constant background offset correction';
    num_lines = 1;
    def = {'0.00'};
    answer = inputdlg(prompt,dlg_title,num_lines,def);
    bg = str2num(answer{1});
    if isempty(bg)
        error('invalid entry');
    end
    
    % get the image stack and subtract background
    X = getappdata( gcbf, 'imageStackX' );
    X = X - bg;
    setappdata( gcbf, 'imageStackX', X );
    
    
    % refresh the image
    im = findobj(gcbf,'Tag','theImage');
    ik = getappdata(gcf,'currentSlice');
    set(im,'CData', X(:,:,ik));
    drawnow;
    
    disp('subtracted')
    
end

function mnuBgSubtractAuto( src, evt )
    % get the image stack
    X = getappdata( gcbf, 'imageStackX' );
    
    % background from average of first and last slices
    Xbg = X(:,:,[1,end]);
    bg = sum(Xbg(:)) / length(Xbg(:));
    
    X = X - bg;
    setappdata( gcbf, 'imageStackX', X );
end

% display menu callbacks
function mnuDisplayDepth(src, evt)
    checked = get(src,'checked');
    
    % toggle checkbox
    switch checked
        
        case 'on'                       % if it was previously checked
            set(src,'checked','off');   %    uncheck the menu item
                                        %    and hide depth gauge
            set( findobj(gcbf,'tag','txtDepth'),...
                 'visible', 'off' );
             
        case 'off'                       % if it was previously unchecked
             set(src,'checked','on');   %    check the menu item
                                         %    and show depth gauge
            set( findobj(gcbf,'tag','txtDepth'),...
                 'visible', 'on' );
    end
            
end

function mnuDisplay3D(src, evt)
    % change display to 3D stack view
    X = getappdata(gcbf,'imageStackX');
    zPos = getappdata(gcbf,'zPos');
    fileName = getappdata(gcbf,'fileName');
    header = getappdata(gcbf, 'header');
    displayMode = '3D';
    
    puprisa_viewChannelZStack( X, zPos, fileName, header, displayMode );
end

% ROI context menu callbacks
function mnuNewROI( src, evt )
    newROI(findobj(gcbf,'Tag','imageAxes'),...
                    findobj(gcbf,'Tag','theImage'),...
                    findobj(gcbf,'Tag','axZtrend'));
end

function mnuRemoveROI( src, evt )
    iROI = getappdata(gco,'iROI');
    
    removeROI( iROI );
    
end

function mnuSetROIColor( src, evt )
    % pick a new color
    c = uisetcolor(gco);
    
    iROI = getappdata(gco,'iROI');
    ROIs = getappdata(gcbf,'ROIs');
    
    ROI = ROIs(iROI);
    hColoredObjs = [ ...
        ROI.lCrossX, ...
        ROI.lCrossY, ...
        ROI.lCrossCenter, ...
        ROI.lBox, ...
        ROI.lBoxCorner, ...
        ROI.lZTrend ];
    
    for iObj = 1:length(hColoredObjs)
        set(hColoredObjs,'Color',c);
    end
    
    drawnow;
    
    
end

 

%==========================================================================
% DATA DISPLAY/REFRESH/DRAWING
%
% These routines update the various displays of the data

function X = demoDataSet()
    % demo dataset: random noise multiplied by exp. decay in the time
    % domain
    m = 256; n = 256; p = 32;
    
    X = rand( m, n, p );
    
    decay = exp(-(1:p)/20);
    
    decay = repmat(reshape(decay,1,1,p),[m,n,1]);
    
    X = X.*decay;
end


function [X,Y, zPos] = loadData( fileName )
    hAxImage = findobj(gcf,'Tag','imageAxes');

    slices = readImageStack( fileName, hAxImage );
    
    
    % allocate image array
    [mx,nx] = size(slices(1).imageData{1});
    p = length(slices);
    X = zeros(mx,nx,p);
    
    if length(slices(1).imageData) > 1
        [my,ny] = size(slices(1).imageData{2});
        Y = zeros(my,ny,p);
    else
        Y = [];
    end

    
    % just pull out x-channel data for now
    for iSlice = 1:length(slices)
        X(:,:, iSlice) = slices(iSlice).imageData{1};
        if ~isempty(Y)
            Y(:,:, iSlice) = slices(iSlice).imageData{2};
        end
    end
    
    zPos = [slices.posZ];
end

function mnuExportStack(evt,obj)
    % exports raw stack data as a TIFF

    if ispref( 'puprisa', 'tiffExportDir' )
        wd = getpref( 'puprisa', 'tiffExportDir' );
    else
        wd = pwd;
        addpref( 'puprisa', 'tiffExportDir', wd );
    end
    
    % prompt for filename
    [fileName, pathName] = uiputfile('*.tif','Save TIFF stack as:',...
        [wd,filesep,'zstack.tif']);
    
    if fileName == 0; return; end;  % stop if user clicked 'cancel'
    
    setpref( 'puprisa', 'tiffExportDir', pathName );
    
    stack = getappdata(gcbf,'imageStackX');
    [~,~,nslices] = size(stack);
    
    % go to first slice
    updateDisplay( 1 );
    
    % write first slice
    imwrite(getimagergb(findobj(gcbf,'tag','imageAxes')),...
        [pathName,filesep,fileName],'tif',...
        'WriteMode','overwrite');
    
    % append the rest of the slices
    for iSlice = 2 : nslices
        
        updateDisplay( iSlice);
        
        imwrite(getimagergb(findobj(gcbf,'tag','imageAxes')),...
            [pathName,filesep,fileName],'tif',...
            'WriteMode','append');
    end
    
end

function mnuExportGIF(evt,obj)
    % get image stack data
    stack = getappdata(gcbf,'imageStackX');
    [nrows,ncols,nslices] = size(stack);
    
    % export stack to animated GIF
    %
    % prompt user for parameters:
    % frame from, frame to, loop, and delay
    exportParm = puprisa_exportGIFDialog( nslices );
    
    nSlicesToAnimate = exportParm.endSlice - exportParm.startSlice + 1;
    
    % if animation will run forwards then backwards,
    % double the no. of frames.
    if ~exportParm.backForth
        nFrames = nSlicesToAnimate;
    else
        nFrames = 2*nSlicesToAnimate;
    end
    
    gifStack = zeros(nrows,ncols,1,nFrames);
    
    % set up colormap
    % we need to add white to the image colormap under use
    cmInUse = colormap();
    cmInUse(1,:) = [1,1,1];
    
    
    
    for iFrame = 1 : nSlicesToAnimate
        
        updateDisplay( iFrame + exportParm.startSlice - 1);
        
        frame = getframe(gca);
        
        % set colormap on first frame
        if iFrame == 1
            [gifStack,cm] = ...
                rgb2ind(frame.cdata,cmInUse);
        else
            gifStack(:,:,1,iFrame) = ...
                rgb2ind(frame.cdata,cmInUse);
        end
    end
    
    % go backwards
    if exportParm.backForth
        for iFrame = 1:nSlicesToAnimate
            
            updateDisplay( exportParm.endSlice - iFrame + 1);

            frame = getframe(gca);
            
            gifStack(:,:,1,iFrame+nSlicesToAnimate) = ...
                rgb2ind(frame.cdata,cmInUse);
        end
    end
    
    
    if exportParm.loop == 1
        loopCount = Inf;
    else
        loopCount = 0;
    end
    
    imwrite(gifStack,cmInUse,exportParm.fileName,'gif',...
        'DelayTime',exportParm.delay,'LoopCount',Inf);
        
    %web(exportParm.fileName);
    %close(ff);
end

function mnuExportMontage(evt,obj)
    % get image stack data
    stack = getappdata(gcbf,'imageStackX');
    [~,~,nFrames] = size(stack);
       
    cmInUse = colormap();
    
    updateDisplay(2);
    frame = getframe(gca);
    [nrows,ncols,~] = size(frame.cdata);
    gifStack = zeros(nrows,ncols,1,nFrames);
        
    for iFrame = 1 : nFrames
        updateDisplay(iFrame);
        frame = getframe(gca);
        gifStack(:,:,:,iFrame) = rgb2ind(frame.cdata,cmInUse);
    end
    prompt = {'rows','columns'};
    dlg_title = 'Size of montage:';
    num_lines = 1;
    def = {'5', '5'};
    answer = inputdlg(prompt,dlg_title,num_lines,def);
    
    m = str2num( answer{1} );
    n = str2num( answer{2} );
    
    figure
    h = montage(gifStack,cmInUse, 'Size', [m n]);
    imageOut  = get(h,'cdata');
    imwrite(imageOut,cmInUse,'montage.png','png');
   
end


function mnuMaxProjections(evt,obj,xypt)

if nargin<3, xypt = []; end

            X = getappdata( gcbf, 'imageStackX' );
            zPos = getappdata(gcbf,'zPos');
            m = mean(X(:));
            s = 3*std(X(:));
            
             ZP = squeeze(max(X(:,:,:),[],3))';
%             ZPp = squeeze(max(X(:,:,:),[],3))';
%             ZPn = squeeze(max(-X(:,:,:),[],3))';
%             ZP = (ZPp>ZPn).* ZPp + (ZPn>ZPp).*(-1*ZPn);
             
            YP = squeeze(max(X(:,:,:),[],1))';
%             YPp = squeeze(max(X(:,:,:),[],1))';
%             YPn = squeeze(max(-X(:,:,:),[],1))';
%             YP = (YPp>YPn).* YPp + (YPn>YPp).*(-1*YPn);
            
            XP = squeeze(max(X(:,:,:),[],2))';
%             XPp = squeeze(max(X(:,:,:),[],2))';
%             XPn = squeeze(max(-X(:,:,:),[],2))';
%             XP = (XPp>XPn).* XPp + (XPn>XPp).*(-1*XPn);
            
            
            
            figure,
%             set(gcf,'Position',[416 546 1118 281])
                subplot(4,4,[1 2 3 5 6 7 9 10 11])
                imagesc([],[],ZP',[-1,1]*max([abs(m+s),abs(m-s)]) );
                title('Z Projection, X-Y image ')
                axis xy
                set(gca,'XTick',[])
%                 xlabel('X (pixels)'),
                ylabel('Y (pixels)')
                
                subplot(4,4,[13 14 15])
                imagesc([],zPos,YP, [-1,1]*max( [abs(m+s),abs(m-s)]) );
                title('Y Projection, X-Z image ')
                xlabel('X (pixels)'),ylabel('Z (µm)')
                
                subplot(4,4,[4 8 12])
                imagesc(zPos,[],XP', [-1,1]*max( [abs(m+s),abs(m-s)]) );
                title('X Projection, Y-Z image ')
                set(gca,'YTick',[])
%                 ylabel('Y (pixels)'),
                xlabel('Z (µm)')
                axis xy
                colorbar; %axis square;
                puprisa_colorMap('pumpProbe');
                
                
                
                 if isempty(xypt)~=1;
                    subplot(4,4,[1 2 3 5 6 7 9 10 11])
                        line([xypt(1) xypt(1)], [1 size(ZP,1)], ...
                            'Color', [1 1 1],'LineWidth', 1, 'LineStyle', '--')
                        line([1 size(ZP,2)],[xypt(2) xypt(2)], ...
                            'Color', [1 1 1],'LineWidth', 1, 'LineStyle', '--')
                    
                    subplot(4,4,[13 14 15])
                        line([xypt(1) xypt(1)], [zPos(1) zPos(end)], ...
                            'Color', [1 1 1],'LineWidth', 1, 'LineStyle', '--')
                    
                    subplot(4,4,[4 8 12])
                        line([zPos(1) zPos(end)],[xypt(2) xypt(2)], ...
                                'Color', [1 1 1],'LineWidth', 1, 'LineStyle', '--')
                end    
                
end


function axesButtonDownFcn( src, evt )
   sel_typ = get(gcbf,'SelectionType');
   
   switch sel_typ
        case 'open'
            
            [x,y] = ginput(1);
             X = getappdata( gcbf, 'imageStackX' );
             zPos = getappdata(gcbf,'zPos');
             m = mean(X(:));
             s = 3*std(X(:));
             
             mnuMaxProjections([],[],round([x,y]))
             
             figure,
             set(gcf,'Position',[12 38 1028 201])
                subplot(121)
                imagesc([],zPos,squeeze(X(round(y),:,:))',[-1,1]*max( [abs(m+s),abs(m-s)] ));
                title(['X-Z image, at y = ', num2str(round(y)),' pixel'])
                xlabel('X (pixels)'),ylabel('Z (µm)')
                colorbar; %axis square;
                subplot(122)
                imagesc([],zPos,squeeze(X(:,round(x),:))',[-1,1]*max( [abs(m+s),abs(m-s)] ));
                title(['Y-Z image, at x = ', num2str(round(x)),' pixel'])
                xlabel('Y (pixels)'),ylabel('Z (µm)')
                colorbar; %axis square;
                puprisa_colorMap('pumpProbe');
                
                
                           
   end
   
   
end
