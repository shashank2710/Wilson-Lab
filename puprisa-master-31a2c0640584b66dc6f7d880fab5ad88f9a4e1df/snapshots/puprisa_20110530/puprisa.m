% puprisa
% pump-probe image analysis
% Jesse Wilson (2011) syrex314@gmail.com
% Duke University
%
% A visualization tool for pump--probe delay stacks of images.
% Shows one slice at a time, and allows user to view the time delay trace
% associated with a selected pixel in the slice. Also allows user to change
% the currently viewed slice
%
% The purpose is to have a tool for visualizing time delay scans at each
% pixel.
%
% assume at most 4 channels
%
% X is an n by m by l - element matrix where
% m is the number of pixels in the y-direction
% n is the number of pixels in the x-direction
% p is the number of time delays (the number of image slices)
%
%
% CHANGELOG
% 04/28/2011: Modified double-click behavior to open up 
%             channel in either z-stack or delay-stack viewer, depending
%             on how the stack was recorded

function puprisa( X, Y, zPos, fileName )
    %close all;
    
    % set up the figure with 4 panels and a statusbar
    f = figure('Name','puprisa','Units','pixels');
    initFigure( f );
    
    initMenus( f );
    
    updateStatus('puprisa ready.');
    
end

%--------------------------------------------------------------------------
% INITIALIZATION ROUTINES
%
function initFigure( f )
    set(f,'Units','pixels');
    fpos = get(f,'Position');
    fwidth = fpos(3);
    
    % other figure options
    set(f,'Toolbar','figure');

    % set up channel views
    nChanRow = 2; 
    setappdata(gcf,'nChanRow',nChanRow);
    nChanCol = 2;
    setappdata(gcf,'nChanCol',nChanCol);
    setappdata(gcf,'maxNChannels',nChanRow*nChanCol);
    
    hChannelAxes = zeros(1,nChanRow*nChanCol);
    hChannelImages = zeros(1,nChanRow*nChanCol);

    for iChannel = 1:(nChanRow*nChanCol);
        hChannelAxes(iChannel) = axes('Units','pixels');
        hChannelImages(iChannel) = image(0,'CDataMapping','scaled',...
                                    'ButtonDownFcn',@axesButtonDownFcn);
        set(hChannelAxes(iChannel), 'PlotBoxAspectRatio',[1 1 1],...
                                    'DataAspectRatioMode','auto',...
                                    'CLimMode','auto',...
                                    'YDir','normal',...
                                    'ButtonDownFcn',@axesButtonDownFcn);
                                
        title(['Channel ',num2str(iChannel)]);
    end
    
    setappdata(gcf,'hChannelAxes',hChannelAxes);
    setappdata(gcf,'hChannelImages',hChannelImages);
    
    % set up status bar
    hStatusText = uicontrol('style','text','String','ImageAnalysisJ',...
        'Units','pixels','tag','hStatusText','position',[0,0,fwidth,20]);

    % resize everything
    setChannelAxesPos(gcf);
    % set up callbacks
    set(gcf,'ResizeFcn',@windowResizeFcn)
       
    % set up scroll wheel fn to be able to flip between slices
    set(gcf,'WindowScrollWheelFcn', @scroll_wheel);
end

% set up colorbar for axis and its associated image stack data
function  initColorBar( ax, appDataTag_imageStack )
    % ensure requested axis is current
    if gca ~= ax
        axes( ax )
    end
    
    cb = colorbar('Location','WestOutside');
    drawnow;
    cbpos = get(cb,'position');
    cbpos = cbpos - [0.1,0,0,0];
    set(cb,'position',cbpos);
    
    % set up callbacks for colorbar
    %set(cb,'ButtonDownFcn',@(src,evt) buttonDown_colorbar( src,evt,ax,  appDataTag_imageStack) );
    h = uicontextmenu();
    uimenu(h,'Label','Autoscale whole stack','Callback',@mnuAutoScaleStack);
    uimenu(h,'Label','Autoscale individual slices');
    uimenu(h,'Label','Manual scale...');
    uimenu(h,'Label','Zero-balanced symmetric autoscale','Separator','on',...
        'Checked','on');
    uimenu(h,'Label','Full autoscale','Separator','on');
    uimenu(h,'Label','Histogram std. dev. autoscale...','Checked','on');
    
    set(cb,'UIContextMenu',h);
end

% Menu set-up
function initMenus( fig )
    % sets up menus in the figure
    f = uimenu('Label','ImageAnalysis');
    uimenu(f,'Label','Open Image Stack...','Callback',@mnuOpenImageStack);
    uimenu(f,'Label','Frame shift correction...','Callback',@mnuFrameShiftCorrection);
    uimenu(f,'Label','Export previews...','Callback',@mnuExportPreviews);
    uimenu(f,'Label','About puprisa','Callback',@mnuAbout);
end

%--------------------------------------------------------------------------
% CALLBACK FUNCTIONS
%
function windowResizeFcn( hObj, evt )
    % resize function
    
    % arrange into 4 panels with info text bar at bottom
    
    % get info on new window size
    set(gcbf,'Units','pixels');
    fpos = get(gcbf,'Position');
    fwidth = fpos(3);
    
    % set statusbar width
    hStatusText = findobj('tag','hStatusText');
    statPos = get(hStatusText, 'Position');
    set(hStatusText, 'Position', [0, 0, fwidth, statPos(4)]);
    
    
    setChannelAxesPos(gcbf);
    
end

function scroll_wheel(src,evnt )
    ik = getappdata(gcbf,'currentSlice') + evnt.VerticalScrollCount;
    
    nSlices = getappdata(gcbf,'nSlices');
    
    if ik < 1
        ik = 1;
    elseif ik > nSlices
        ik = nSlices;
    end


    setappdata(gcbf,'currentSlice', ik);
    updateAll( gcbf );
end

function axesButtonDownFcn( src, evt )
   sel_typ = get(gcbf,'SelectionType');
    switch sel_typ
        case 'open'
            % double-click -- open analysis window on selected channel
            h = src;
            
            % get parent axes if the user clicked on the image instead
            if(strcmp(get(h,'Type'),'image'))
                h = get(h,'Parent');
            end
            
            ax = h;
            
            fileName = getappdata(gcbf,'fileName');

            imageStack = getappdata(ax,'imageStack');
            [~,~,p] = size(imageStack);

            % determine whether a z stack or delay stack
            stackType = getappdata(gcbf,'stackType');
            
            switch stackType
                case 'delay stack'
                    zpos = getappdata(gcbf,'delays');
                    puprisa_viewChannel( imageStack, zpos, fileName );
                case 'z stack'
                    zpos = getappdata(gcbf,'zPos');
                    puprisa_viewChannelZStack(imageStack, zpos, fileName);
            end
    end
end

% menu callbacks

function mnuAbout(src, evt)
    msgbox('PUPRISA: PUmp PRobe Image Stack Analysis. (2011 Jesse Wilson)',...
        'About PUPRISA', 'help');
end

function mnuAutoScaleStack(src, evt)
    get(get(get(src,'Parent'),'Parent'))
end

function buttonDown_colorbar( src,evt,ax,  appDataTag_imageStack )
    disp('clicked on color bar!');
end

function mnuOpenImageStack( src, evt )
    % prompt user for an image stack to load
    
    % recall working directory
    wd = getappdata(gcbf,'workingDirectory');
    
    % or use MATLAB's working directory if none has been set yet
    if isempty(wd)
        wd = 'C:\Users\jw295\Data';
    end
    
    [fileName, pathName] = uigetfile([wd,'\*.dat'],...
        'Select a ScanImage Stack' );
    
    if ~isequal(fileName,0)
        % remember the working directory for later
        setappdata(gcbf,'workingDirectory',pathName);
        
        % if not canceled, load the file
        loadData(gcbf, [pathName,'\',fileName]);
        
    end
end

function mnuFrameShiftCorrection( src, evt )
    % prompt for alignment channel
    alignChanStr = inputdlg('Frame shift registration.',...
        'Channel number for calculating alignment:', 1, {'1'});
    alignChan = str2num(alignChanStr{1});
    
    if isempty( alignChan )
        error('invalid channel number');
    end
    
    % align slices based on cross-correlation of y-channel
    hChannelAxes = getappdata(gcbf,'hChannelAxes');

    X = getappdata(hChannelAxes(1),'imageStack');
    Y = getappdata(hChannelAxes(alignChan),'imageStack');
    %Y = Y(1:256,1:256,:);      % reduce size
    
    figure();
    [nrows,ncols,nslices] = size(Y);
    
    %for islice = 1:nslices
    %    imagesc(Y(:,:,islice));
    %    drawnow;
    %end
    
    for islice = 2:nslices
        X2 = X(:,:,islice);
        Y2 = Y(:,:,islice);
        %Y2 = circshift(Y1,[20,13]);
        
        

        [output, shifted] = dftregistration(...
            fft2(Y(:,:,1)), fft2(Y2),100);
        
        Y2 = abs(ifft2(shifted));
        Y(:,:,islice) = Y2;
        
        % shift X accordingly
        %X2 = circshift(X2,[shifty,shiftx]);
        %
        % copied from dft registration code
        diffphase = output(2);
        row_shift = output(3);
        col_shift = output(4);
        [nr,nc]=size(X2);
        Nr = ifftshift([-fix(nr/2):ceil(nr/2)-1]);
        Nc = ifftshift([-fix(nc/2):ceil(nc/2)-1]);
        [Nc,Nr] = meshgrid(Nc,Nr);
        Greg = fft2(X2).*exp(i*2*pi*(-row_shift*Nr/nr-col_shift*Nc/nc));
        Greg = Greg*exp(i*diffphase);
        X2 = real(ifft2(Greg));
        X(:,:,islice) = X2;
        
        
    end
    
    %for islice = 1:nslices
    %    imagesc(Y(:,:,islice));
    %    drawnow;
    %end
    
    % save the results
    setappdata(hChannelAxes(alignChan),'imageStack',Y);
    setappdata(hChannelAxes(1),'imageStack',X);

end

function mnuExportPreviews(src, evt)
    previewDir = 'C:\Users\jw295\projects\imageDatabase\previews';
    hdr = getappdata(gcbf,'fileHeader');
    date = hdr.date;
    
    previewDir = [previewDir,'\',date];
    
    if exist(previewDir) == 0
        mkdir(previewDir);
    end
    
    if exist(previewDir) ~= 7
        error('Could not make directory');
    end
    
    filename = hdr.filename;
    [s,t]=regexp(filename,'\\([^\\]+)\.dat','tokens','match');
    filename = [previewDir,'\',s{:}{1},'.png']
    
    print('-dpng','-r72', filename);
    
end

%--------------------------------------------------------------------------
% CALLBACK SUPPORT FUNCTIONS
%
function setChannelAxesPos( f )
    hChannelAxes = getappdata(f,'hChannelAxes');
    hStatusText = findobj('tag','hStatusText');
    statusTextPos = get(hStatusText,'position');
    statusTextHeight = statusTextPos(4);
    
    nChanRow = getappdata(f,'nChanRow');
    nChanCol = getappdata(f,'nChanRow');

    % assume figure units are pixels
    pos = get(f,'Position');
    figWidth = pos(3);
    figHeight = pos(4) - statusTextHeight;
    bottom = statusTextHeight;
    
    chanWidth = round(figWidth)/nChanCol;
    chanHeight = round(figHeight)/nChanRow;
    
    nChan = length( hChannelAxes );
    
    for iChan = 1:nChan
        hAx = hChannelAxes(iChan);
        
        % start at top left, go right, then down
        row = floor((iChan+nChanCol-1) / nChanCol);
        col = mod(iChan-1, nChanCol) + 1;
        
        tightInset = get(hAx, 'tightInset');
        vMargin = tightInset(2)+tightInset(4);
        hMargin = tightInset(1)+tightInset(3);
        
        axHeight = chanHeight - vMargin;
        axWidth = chanWidth - hMargin;
        
        posX = (col-1)*chanWidth + tightInset(1);
        posY = (nChanRow - row)*chanHeight + tightInset(2) + bottom;
        
        
        set(hAx, 'position', [posX, posY, axWidth, axHeight]);
    end
    
    
end

function loadData( f, fileName )
% load data, then parse into figure data structures, and call update
% routines
%
% ACTION: Modify this to read header too; look at scanaxis to tell whether
%         it was a z stack or a delay stack. 
%         DelayStack: scanaxis == 0
%         ZStack:     scanaxis == 1
    setappdata(f,'fileName',fileName);

    set(f,'pointer','watch');
    drawnow;
    figLoading = figure('WindowStyle','modal','Name','Loading...');
    hAxImage = axes();
    [slices, header] = readImageStack( fileName, hAxImage );
    nSlices = length(slices);
    close(figLoading);
    
    if isempty( slices )
        set(f,'pointer','arrow');
        return;
    end
    
    
    % allocate image stack arrays
    nChannels = length(slices(1).imageData);
    setappdata(f,'nChannels',nChannels);
    setappdata(f,'nSlices',nSlices);
    setappdata(f,'fileHeader',header)
    
    [nRows,nCols] = size(slices(1).imageData{1});
    
    imageStackThisChannel = zeros(nRows, nCols);
    for iChannel = 1:nChannels
        % pull out each channel
        for iSlice = 1:nSlices
            imageStackThisChannel(:,:, iSlice) ...
                = slices(iSlice).imageData{iChannel};
        end
        imageStackChannels{iChannel} = imageStackThisChannel;
    end
        
    % save all image stacks
    setappdata(f,'imageStackChannels',imageStackChannels);
    
    % save delay axes
    setappdata(f, 'delays', cell2mat({slices.delays}));
    setappdata(f, 'zPos', cell2mat({slices.posZ}));
    
    % determine whether this was a delay stack or z stack
    if isfield( header, 'scanaxis' )
        % newer header format
        isZStack = header.scanaxis;
    elseif isfield( header, 'variableaxisZt' )
        % old header format
        isZStack = 1 - header.variableaxisZt;
    else
        error('Scan axis not specified in header! Header should contain scanaxis or variable axis (Z/t).')
    end
    
    if isZStack
        stackType = 'z stack';
    else 
        stackType = 'delay stack';
    end
    setappdata(f,'stackType',stackType);
        

    % default -- view first slice
    setappdata(f,'currentSlice',1);
    
    % change sizes of arrays appropriately
    hChannelImages = getappdata(f,'hChannelImages');
    hChannelAxes = getappdata(f,'hChannelAxes');

    nChannels = length(hChannelImages);
    for iChannel = 1:nChannels
        if iChannel <= length(imageStackChannels)
            imageStackThisChannel = imageStackChannels{iChannel};
        else
            imageStackThisChannel = zeros(nRows,nCols,nSlices);
        end
        
        % store image stack in app data for this axes object
        setappdata(hChannelAxes(iChannel), ...
            'imageStack', imageStackThisChannel);
        
        set(hChannelImages(iChannel),'XData',[1,nCols],'YData',[1,nRows],...
            'CData',squeeze(imageStackThisChannel(:,:,1)));
        set(hChannelAxes(iChannel),'XLim',[1,nCols],'YLim',[1,nRows]);
        
    end

    % display all channels for current slice
    updateAll( gcf );
    set(f,'pointer','arrow');
end


%--------------------------------------------------------------------------
% UPDATE/REDRAW FUNCTIONS
%
function updateStatus( str )
    set(findobj('tag','hStatusText'), 'String', str);
    drawnow;
end

function updateAll( f )
    % find out the current slice and change cdata on all axes to match
    % default -- view first slice
    currentSlice = getappdata(f,'currentSlice');
    nSlices = getappdata(f,'nSlices');

    
    % change sizes of arrays appropriately
    hChannelImages = getappdata(f,'hChannelImages');
    hChannelAxes = getappdata(f,'hChannelAxes');

    nChannels = length(hChannelImages);
    for iChannel = 1:nChannels
        imageStack = getappdata(hChannelAxes(iChannel),'imageStack');
        imageSlice = squeeze(imageStack(:,:,currentSlice));
        
        set(hChannelImages(iChannel),'CData',imageSlice);
    end
    
    fileName = getappdata(f,'fileName');
    updateStatus([fileName,': slice ', num2str(currentSlice), ' / ', ...
        num2str(nSlices)]);
    
    drawnow;
end