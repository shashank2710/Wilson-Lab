% puprisa_viewChannel
% Jesse Wilson (2010) syrex314@gmail.com
% Duke University
%
% View a single channel
%
% A visualization tool for 3D stacks of images, such as those acquired in
% pump--probe microscopy.
% Shows one slice at a time, and allows user to view the time delay trace
% associated with a selected pixel in the slice. Also allows user to change
% the currently viewed slice
%
% The purpose is to have a tool for visualizing time delay scans at each
% pixel.
%
% X is an n by m by l - element matrix where
% m is the number of pixels in the y-direction
% n is the number of pixels in the x-direction
% p is the number of time delays (the number of image slices)
%
%
% SOMEDAY: make callback stubs for ROI context menu items
% SOMEDAY: implement callbacks for cbar autoscale context menu
% Then clean up this code before going any further...
% TODO : check bounds on ROI movement
% TODO: make ROI tool with uitoolbar()
% consider making ROI an object?
% start work on multiexponential analysis 
% PCA analysis and basis fitting

function puprisa_viewChannel( X, zPos, fileName )
    figure('Name','Channel');
    
    % maximize the window horizontally
    pos = get(gcf,'OuterPosition');
    screenSize = get(0,'ScreenSize');
    pos(1) = screenSize(1);
    pos(3) = screenSize(3);
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
    
    % set up 4 subplots
    % 1) image for selected time delay
    % 2) delay trace for selected pixel
    % 3) alternate channel image
    % 4) unused
    
    ik = 1; % ik is slice
    setappdata(gcf,'currentSlice',ik);

    
    axes('Position',[0.1,0.1,0.4,0.8]);
    hImageAxes = gca();
    axXChanImage = gca();
    im = imagesc(X(:,:,ik));
    set(im,'Tag','theImage');
    axis square;
    %colormap(cmap(256,[0,0.8,0],[0,0,0],[1,0.0,0.2]));
    %colormap(cmap(256,[0.5,0.5,1],[0,0,1],[0,0,0],[1,0,0],[1,0.5,0.5]));
    colormap(cmap(256,[0,1,1],[0,0,1],[0,0,0],[1,0,0],[1,1,0]));
    initColorBar( axXChanImage, 'imageStackX' );

    set(gca,'YDir','normal');
    set(gca,'Tag','imageAxes');
    setappdata(gcf,'imageStackX', X);
    setappdata(gcf,'zPos',zPos);
    
    if nargin == 3
        title(fileName,'interp','none');
        setappdata(gcf,'fileName',fileName);
    else
        title('viz3dStack');
    end
    
    axes('Position',[0.55,0.1,0.4,0.8],'tag','axTransient');
    setappdata(gcf,'zPos',zPos);
    hTransientAxes = gca();
    
    
    
    % set up controls for region of interest
    newROI(hImageAxes,im, hTransientAxes);
    newROI(hImageAxes,im, hTransientAxes);

     

   
    % set up scroll wheel fn to be able to flip between slices
    set(gcf,'WindowScrollWheelFcn', @scroll_wheel);
    
    % attempt auto-clim
    % take the biggest difference between mean and zero
    m = mean(X(:))
    s = 4*std(X(:))
    %s = max(abs(X(:)));
    set(hImageAxes,'clim', [-1,1]*max( [abs(m+s),abs(m-s)] ));
    
end

%--------------------------------------------------------------------------
% INITIALIZATION ROUTINES
%
function newROI( hImageAxes, hImage, hLineoutAxes )
    % set up controls
    % assumes ROI goes in the current figure
    
    % get relevant app data
    t  = getappdata(gcf,'zPos');
    X  = getappdata(gcf,'imageStackX');
    [m,n,p] = size(X);

    ik = getappdata(gcf,'currentSlice');
    
        
    % get roi array and find out this ROI's index
    ROIs = getappdata(gcf,'ROIs');
    ROI.index = 1+length(ROIs)

    
    % a crosshair with point in the middle
    axes(hImageAxes);
    
    % initial pixel and delay
    % center pixel,
    % first delay slice
    ii = round(n/2);
    ij = round(m/2);
    setappdata(gcf,'roiCenter',[ii,ij]);

    lCrossXHL = line([0,0],[0,0],'Color','white','LineWidth',1.5);
    lCrossX = line(xlim(), [0,0]+ij, 'Color', [0,0,1],'Tag','lCrossX','LineWidth',0.5);
    hl = linkprop([lCrossX,lCrossXHL],{'XData','YData'});
    set(lCrossX,'UserData',hl);
    
    lCrossYHL = line([0,0],[0,0],'Color','white','LineWidth',1.5);
    lCrossY = line([0,0]+ii, ylim(), 'Color', [0,0,1],'Tag','lCrossY','LineWidth',.5);
    hl = linkprop([lCrossY,lCrossYHL],{'XData','YData'});
    set(lCrossY,'UserData',hl);
    
    lCrossCenter = line( ii, ij, 'Color',[0,0,1], 'Marker','o',...
        'ButtonDownFcn',@buttonDown_moveRoi,...
        'Tag','lCrossCenter','UIContextMenu',getappdata(gcf,'cmnuROI'));
    
    % label this UI item with the ROI index
    setappdata(lCrossCenter,'iROI',ROI.index);
    
    setappdata(gcf,'ROImode','box'); % may be either box or pixel
    
    % highlight box under ROI box 
    lRoiBoxHL = line([0,0,0,0,0],[0,0,0,0,0],'Color','white','LineWidth',1.5);
    
    % a box around the ROI
    roiWidth = 20;
    roiHeight = 8;
    setappdata(gcf,'roiBox', [roiWidth, roiHeight]);
    lRoiBox = line( [ii-roiWidth/2, ii-roiWidth/2, ...
                     ii+roiWidth/2, ii+roiWidth/2, ii-roiWidth/2],...
                    [ij-roiHeight/2, ij+roiHeight/2, ...
                     ij+roiHeight/2, ij-roiHeight/2, ij-roiHeight/2],...
                 'Color', [1,1,1], 'Tag','lRoiBox','LineWidth',0.5 );
             
    % linke ROI box with its underlying hightlight box
    hl = linkprop( [lRoiBox,lRoiBoxHL], {'XData','YData'});
    set(lRoiBox,'UserData',hl); % store link in userdata
                     
    set(hImage,'ButtonDownFcn',{@(s,e)buttonDown_moveRoi(s,e,ROI.index)});
    
    lRoiBoxCorner = line( ii+roiWidth/2, ij-roiHeight/2, 'Color',[1,1,1],'Marker','o',...
        'Tag','lRoiBoxCorner','ButtonDownFcn', @buttonDown_resizeRoi );
    setappdata(lRoiBoxCorner,'iROI',ROI.index);

    
    % set up line to display the ROI
    y = squeeze(X(ii,ij,:));
    axes( hLineoutAxes );
    lTransient = line(t,y,'Tag','lTransient','Marker','.');
    line(t,0*t,'Color','k','LineStyle','--');
    xlim([min(t),max(t)]);
    xlabel('Pump--probe delay (ps)');
    ylabel('Probe transient absorption (\mu V)');
    lTimeSliceMark = line( [1,1]*t(ik), ...
        get(hLineoutAxes,'ylim'),...
        'color','green', 'tag','lTimeSliceMark');
    
    % gather all relevant data into a roi structure
    ROI.center = [ii,ij];
    ROI.box = [roiWidth, roiHeight];
    ROI.mode = 'box';
    ROI.lCrossX = lCrossX;
    ROI.lCrossY = lCrossY;
    ROI.lCrossCenter = lCrossCenter;
    ROI.lBox = lRoiBox;
    ROI.lBoxCorner = lRoiBoxCorner;
    ROI.lTransient = lTransient;

    % add this ROI
    
    if isempty(ROIs)
        clear ROIs;
    end
    
    
    ROIs(ROI.index) = ROI;
    
    % save the ROIs
    setappdata(gcf,'ROIs',ROIs);
    
end

function removeROI( ROIIndex )
    
    % get the ROI structure array
    ROIs = getappdata(gcf,'ROIs');
    
    ROI = ROIs(ROIIndex);

    % delete lines belonging to the element in question
    delete(ROI.lCrossX);
    delete(ROI.lCrossY);
    delete(ROI.lCrossCenter);
    delete(ROI.lBox);
    delete(ROI.lBoxCorner);
    delete(ROI.lTransient);
    
    
    switch ROIIndex 
        case length( ROIs )
            % special case: last element
            ROIs = ROIs(1:(end-1));
        case 1
            % special case: first element
            ROIs = ROIs(2:end);
        otherwise
            % remove a middle element
            ROIs = [ROI(1:(ROIIndex-1)), ROIs((ROIIndex+1):end)];
    end
    
    % re-index the ROIs
    % and update functoin handles to work with the correct index
    for iROI = 1:length(ROIs)
        ROIs(iROI).index = iROI;
        
        % set iROI of UI elements so we can tell which ROI was clicked on
        setappdata(ROIs(iROI).lCrossCenter,'iROI', iROI);
        setappdata(ROIs(iROI).lBoxCorner,'iROI', iROI);
    end
    
    % save the new ROI array
    setappdata(gcf,'ROIs',ROIs);
    
    % re-draw
    updateCursor();
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
    uimenu(f,'Label','SVD de-noise...','Callback',@mnuSVDdeNoise);
    uimenu(f,'Label','Print stack on page...','Callback',@mnuPrintStack);
    uimenu(f,'Label','Subract background value...','Callback',@mnuBgSubtract);
    uimenu(f,'Label','Auto background','Callback',@mnuBgSubtractAuto);
    uimenu(f,'Label','Export stack','Callback',@mnuExportStack);
    uimenu(f,'Label','Export animated GIF','callback',@mnuExportGIF);
    uimenu(f,'Label','Export image stack to workspace...',...
        'callback',@mnuExportToWorkspace);
    
    % context menu for ROI
    cmnuROI = uicontextmenu();
    uimenu(cmnuROI,'Label','New ROI', 'Callback', @mnuNewROI);
    uimenu(cmnuROI,'Label','Delete ROI', 'Callback', @mnuRemoveROI);
    uimenu(cmnuROI,'Label','Color...','Callback',@mnuSetROIColor);
    uimenu(cmnuROI,'Label','Label...');
    setappdata(fig,'cmnuROI',cmnuROI);

end


%==========================================================================
% CALLBACK FUNCTIONS
%
function buttonDown_colorbar( src,evt,ax,  appDataTag_imageStack )
    disp('clicked on color bar!');
end

function scroll_wheel(src,evnt )
    ik = getappdata(gcbf,'currentSlice') + evnt.VerticalScrollCount;
    
    X = getappdata(gcbf,'imageStackX');
    [m,n,p] = size(X);
    
    Y = getappdata(gcbf,'imageStackY');
    
    if ik < 1
        ik = 1;
    elseif ik > p
        ik = p;
    end
    
    setappdata(gcbf,'currentSlice', ik);
    
    h = findobj(gcf,'Tag','imageAxes');
    
    if gca ~= h
        axes(h);
    end
    
    fn = getappdata(gcf,'fileName');
    title([fn, num2str(ik)],'interp','none');
    
    
    im = findobj(gcbf,'Tag','theImage');
    set(im,'CData', X(:,:,ik));
    
    if ~isempty(Y)
        im = findobj(gcbf,'Tag','imY');
        set(im,'CData', Y(:,:,ik));
    end
    
    lTimeSliceMark = findobj(gcf,'Tag','lTimeSliceMark');
    t = getappdata(gcbf,'zPos');
    axTransient = findobj(gcf,'Tag','axTransient');
    set(lTimeSliceMark,'xdata',[0,0]+t(ik),'ydata',get(axTransient,'ylim'));
    
    
    drawnow;
end

function buttonDown_moveRoi(src,evnt )
% src - the object that is the source of the event
% evnt - empty for this property
    sel_typ = get(gcbf,'SelectionType');
    
    % get ROI struct array
    ROIs = getappdata(gcbf,'ROIs');
    
    % which ROI?
    ROIindex = getappdata(src,'iROI');
    
    % get the ROI that was selected
    ROI = ROIs(ROIindex);

    switch sel_typ
        case 'normal'
            % single click with left button-- track mouse motions
            set(gcf,'WindowButtonMotionFcn',{@(s,e)buttonMotion_moveRoi(s,e, ROI)});
            set(gcf,'WindowButtonUpFcn',@button_up);
        case 'open'
            % double-click -- toggle ROI mode
            % needs to be fixed
            switch ROI.mode
                case 'box'
                    ROI.mode = 'point';
                case 'point'
                    ROI.mode = 'box';
                otherwise % invalid mode-- change back to 'point'
                    disp(['invalid ROI mode: ', ROI.mode, ' -- changing back to point']);
                    ROI.mode = 'point';
            end
            
            % update the ROI structure in appdata
            ROIs(ROIindex) = ROI;
            setappdata(gcbf,'ROIs', ROIs);
            
            % refresh the display after changing ROI mode
            updateCursor();
        case 'extend' % shift click -- make a new ROI
           newROI(findobj(gcbf,'Tag','imageAxes'),...
                findobj(gcbf,'Tag','theImage'),...
                findobj(gcbf,'Tag','axTransient'));

            updateCursor();
    end
end

function buttonMotion_moveRoi(src,evnt, ROI)
    % find out where mouse is
    cp = get(gca,'currentPoint');
    ii = round(cp(1,1));
    ij = round(cp(1,2));
    
    
    ROIs = getappdata(gcbf,'ROIs');
    ROI.center = [ii,ij];
    ROIs(ROI.index) = ROI;
    setappdata(gcbf,'ROIs', ROIs);
    
    updateCursor();
   
end

% buttonDown_resizeRoi()
% Called when user clicks on ROI corner handle to resize.
function buttonDown_resizeRoi(src,evnt )
% src - the object that is the source of the event
% evnt - empty for this property
   sel_typ = get(gcbf,'SelectionType');
   
    % which ROI?
    ROIindex = getappdata(src,'iROI');
    
   
    % get ROI struct array
    ROIs = getappdata(gcbf,'ROIs');
    
    % get the ROI that was selected
    ROI = ROIs(ROIindex);

   switch sel_typ 
      case 'normal'
         set(gcf,'WindowButtonMotionFcn',{@(s,e)buttonMotion_resizeRoi(s,e, ROI)});
         set(gcf,'WindowButtonUpFcn',@button_up);
   end
end

function buttonMotion_resizeRoi(src,evnt, ROI)
    % calculate new width and height
    roiCenter = ROI.center;
    iCenter = roiCenter(1);
    jCenter = roiCenter(2);
    cp = get(gca,'currentPoint');
    ii = round(cp(1,1));
    ij = round(cp(1,2));
    
    roiWidth = 2*(ii-iCenter);
    roiHeight = -2*(ij-jCenter);
    
    % check bounds
    if roiWidth < 1
        roiWidth = 1;
    end
    if roiHeight < 1
        roiHeight = 1;
    end
    
    roiBox = [roiWidth, roiHeight];
    ROI.box = roiBox;
    
    % save the new ROI
    ROIs = getappdata(gcbf,'ROIs');
    ROIs(ROI.index) = ROI;
    setappdata(gcbf,'ROIs',ROIs);

    updateCursor();
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

function mnuSVDdeNoise( src, evt )
    % perform principle component analysis,
    % keeping only the first N principal components!
    
    % ask how many components to retain
    prompt = {'Number of principal components to retain:'};
    dlg_title = 'SVD de-noise parameters';
    num_lines = 1;
    def = {'3'};
    answer = inputdlg(prompt,dlg_title,num_lines,def);
    nPCs = str2num(answer{1});
    if isempty(nPCs)
        error('invalid entry');
    end
    
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
    
    
    [U,S,V] = svd( delayScans, 0 );
    size(U)
    size(S)
    size(V)
    
    s = diag(S);
    
    sumPCStack = zeros(nrows, ncols, ndelays);
    
    %nPCs = 3;
    
    for iPC = 1:nPCs
        % make a matrix of this principal component
        PC_mtrx = permute(repmat( V(:,iPC), [1,nrows,ncols]),[2,3,1]);
    
        % project the image stack onto this principal component
        PCImage = sum(X.*PC_mtrx ,3);
        
        % add it to the sum image
        sumPCStack = sumPCStack + repmat(PCImage,[1,1,ndelays]).*PC_mtrx;
        
    end
    
    % save back in main image
    setappdata (gcbf,'imageStackX', sumPCStack);
    
    src= [];
    evnt.VerticalScrollCount = -1;
    scroll_wheel(src,evnt);
    
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
            puprisa_viewChannel( X, Y, zPos, fileName ) ;
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

% ROI context menu callbacks
function mnuNewROI( src, evt )
    newROI(findobj(gcbf,'Tag','imageAxes'),...
                    findobj(gcbf,'Tag','theImage'),...
                    findobj(gcbf,'Tag','axTransient'));
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
        ROI.lTransient ];
    
    for iObj = 1:length(hColoredObjs)
        set(hColoredObjs,'Color',c);
    end
    
    drawnow;
    
    
end

%==========================================================================
% DATA DISPLAY/REFRESH/DRAWING
%
% These routines update the various displays of the data
function updateCursor()
    % get ROI structure
    ROIs = getappdata(gcf,'ROIs');
    nROIs = length(ROIs);
    
    minVal = inf;
    maxVal = -inf;
    
    for iROI = 1:nROIs
        ROI = ROIs(iROI);
        
    % get relevant handles
    lCrossX = ROI.lCrossX;
    lCrossY = ROI.lCrossY;
    lCrossCenter = ROI.lCrossCenter;
    lRoiBox = ROI.lBox;
    lRoiBoxCorner = ROI.lBoxCorner;
    lTransient = ROI.lTransient;
    ROImode = ROI.mode;
    roiCenter = ROI.center;
    roiBox = ROI.box;


    ii = roiCenter(1);
    ij = roiCenter(2);
    
    % find out whether we operate a boxed ROI or a single pixel
    
    
    switch ROImode
        case 'box'
            % ensure box is all visible
            set(lRoiBox,'visible','on');
            set(lRoiBoxCorner,'visible','on');
        case 'point'
            % hide box, leave only center point
            set(lRoiBox,'visible','off');
            set(lRoiBoxCorner,'visible','off');
        otherwise
            error(['Invalide ROI mode: ', ROImode]);
    end
    
    % update cursor
    set(lCrossX,'ydata',[0,0]+ij);
    set(lCrossY,'xdata',[0,0]+ii);
    set(lCrossCenter,'xdata',ii,'ydata',ij);
    
    % update ROI box
    
    roiWidth = roiBox(1);
    roiHeight = roiBox(2);
    set(lRoiBox,'xdata',[ii-roiWidth/2, ii-roiWidth/2, ...
                     ii+roiWidth/2, ii+roiWidth/2, ii-roiWidth/2]);
    set(lRoiBox,'ydata',[ij-roiHeight/2, ij+roiHeight/2, ...
                     ij+roiHeight/2, ij-roiHeight/2, ij-roiHeight/2]);
                 
    set(lRoiBoxCorner,'xdata',ii+roiWidth/2, 'ydata', ij-roiHeight/2);
                    
    % redraw decay curve on right
    X = getappdata(gcbf, 'imageStackX');
    
    % find lineout of selected region (or pixel)
    switch ROImode
        case 'box'
            % find X within ROI box
            iLower = round(ii - roiWidth/2);
            iUpper = round(ii + roiWidth/2);
            jLower = round(ij - roiHeight/2);
            jUpper = round(ij + roiHeight/2);

            % bound checking
            % do this later-- be careful for now
            Xroi = X(  jLower:jUpper,iLower:iUpper, : );
            XroiSum = squeeze(sum(sum(Xroi,1),2));
            lineout = XroiSum /((jUpper-jLower)*(iUpper-iLower)); % normalize
        case 'point'
            % find transient lineout under point ROI
            Xroi = squeeze(X(ij,ii,:));
            lineout = Xroi;
    end
    
    % record min/max values
    if( max(lineout) > maxVal )
        maxVal = max(lineout);
    end
    if( min(lineout) < minVal )
        minVal = min(lineout);
    end
    
    %set(lTransient,'YData', squeeze(X(ij,ii,:)));
    set(lTransient,'YData', lineout);
    
    lTimeSliceMark = findobj(gcf,'Tag','lTimeSliceMark');
    axTransient = findobj(gcf,'Tag','axTransient');
    %set(axTransient,'ylim',[min(lineout),max(lineout)]);
    set(lTimeSliceMark,'ydata',[minVal,maxVal]);
    
    end
    drawnow;
    
end

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
    stack = getappdata(gcbf,'imageStackX');
    clim = get(findobj(gcbf,'Tag','imageAxes'),'clim');
    
    % subtract offset
    stack = stack - clim(1);
    
    % scale
    stack = stack / (clim(2)-clim(1));
    
    % truncate at 0-255
    stack( stack < 0 ) = 0;
    stack( stack > 1 ) = 1;
    
    [nrows,ncols,nslices] = size(stack);
    
    for iSlice = 1:nslices
        imwrite(stack(:,:,iSlice),['foo_',num2str(iSlice),'.tif'],'tif');
    end
        
    
end

function mnuExportGIF(evt,obj)
    stack = getappdata(gcbf,'imageStackX');
    clim = get(findobj(gcbf,'Tag','imageAxes'),'clim');
    cm = colormap();
    % subtract offset
    stack = stack - clim(1);
    
    % scale
    stack = 255*stack / (clim(2)-clim(1));
    
    % truncate at 0-255
    stack( stack < 0 ) = 0;
    stack( stack > 255 ) = 255;
    
    [nrows,ncols,nslices] = size(stack);
    
    gifStack = zeros(nrows,ncols,1,2*nslices);
    ff = figure(50);
    colormap(cm);
    for iSlice = 1:nslices
        gifStack(:,:,1,iSlice) = stack(:,:,iSlice);
        imagesc(stack(:,:,iSlice));
        drawnow;
    end
    
    for iSlice = 1:nslices
        gifStack(:,:,1,(nslices+iSlice)) = stack(:,:,(nslices-iSlice+1));
        imagesc(stack(:,:,iSlice));
        drawnow;
    end
    
    imwrite(gifStack,cm,['exported','.gif'],'gif',...
        'DelayTime',0,'LoopCount',Inf);
        
    close(ff);
end

function mnuExportToWorkspace( evt, obj )
    % export image stack to workspace variable specified by user
    
    % prompt user for variable:
    varName = inputdlg('Variable name for export:',...
        'Export image stack', 1, {'X'});
    
    stack = getappdata(gcbf,'imageStackX');
    
    assignin('base',varName{1},stack);
    
    
end