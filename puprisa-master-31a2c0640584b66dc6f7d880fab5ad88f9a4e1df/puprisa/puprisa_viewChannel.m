% puprisa_viewChannel
% Jesse Wilson (2010--2013) syrex314@gmail.com
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
% INPUT ARGUMENTS:
%
%   X is the image stack, an n by m by l - element matrix where
%     m is the number of pixels in the y-direction
%     n is the number of pixels in the x-direction
%     p is the number of time delays (the number of image slices)
%
%   zPos is the z position (or time delay)
% 
%   fileName is the .dat file from which the stack was loaded
%
%   header is the text of the .dat file header
%
%   chan is the channel number (e.g. 1 for lock-in X, 2, for lock-in Y,
%       etc)
%
% SOMEDAY: make callback stubs for ROI context menu items
% SOMEDAY: implement callbacks for cbar autoscale context menu
% Then clean up this code before going any further...
% TODO : check bounds on ROI movement
% TODO: make ROI tool with uitoolbar()
% consider making ROI an object?
%
% CHANGELOG:
% 06-07-2011: Sped up user interaction and redrawing by setting figure
%             BusyAction property to 'cancel' while handling motion events
% 06-22-2011: Draw scale bar.

function puprisa_viewChannel( X, zPos, fileName, header, chan )
    
    if nargin > 2
        % set title of figure window if user supplied enough information to
        % do so
        switch filesep
            case '\'
                s = regexp(fileName,'\\([^\\]+$)','match');
            case '/'
                s = regexp(fileName,'/([^/]+$)','match');
        end
        if isempty(s)
            s = fileName;
        else
            s = s{1};
            s = [s(2:end), ' (Ch ', num2str(chan),')'];
        end

        figure('Name',s);
    else
        header = [];
    end
    
    % save header info
    setappdata(gcf,'header',header);
    
    % set print size to fit a lab notebook
    printSize = [6+3/8, 3.0];
    lmargin = (8.5 - printSize(1))/2;
    tmargin = (11 - printSize(2))/2;
    set(gcf,'PaperUnits','inches',...
            'PaperSize',[8.5 11],...
            'PaperPosition',[lmargin, tmargin, printSize]);
    
    
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

    %-----------------------------------------------
    % Make an axis in which to display the image
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
    
    if ~isempty(header)
        puprisa_initScaleBar( gca, header );
    end

    set(gca,'YDir','normal');
    set(gca,'Tag','imageAxes');
    setappdata(gcf,'imageStackX', X);
    setappdata(gcf,'zPos',zPos);
    
    if nargin >= 3
        title(fileName,'interp','none');
        setappdata(gcf,'fileName',fileName);
    else
        title('viz3dStack');
    end
    
    % -----------------------------------------------
    % Make an axes in which to display ROI lineouts
    axes('Position',[0.55,0.1,0.4,0.8],'tag','axTransient');
    setappdata(gcf,'zPos',zPos);    % in picoseconds
    hTransientAxes = gca();
    
    
    % set up controls for region of interest
    newROI(gcf, hImageAxes,im, hTransientAxes);
    % newROI(gcf, hImageAxes,im, hTransientAxes); 

   
    % set up scroll wheel fn to be able to flip between slices
    set(gcf,'WindowScrollWheelFcn', @scroll_wheel);
    
    
    % set up keypress handler for flipping between slices
    set(gcf,'KeyPressFcn', @keyPress);
    
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
function newROI( fig, hImageAxes, hImage, hLineoutAxes )
    % set up controls
    % assumes ROI goes in the current figure
    
    % get relevant app data
    t  = getappdata(fig,'zPos');
    X  = getappdata(fig,'imageStackX');
    [m,n,p] = size(X);

    ik = getappdata(fig,'currentSlice');
    
        
    % get roi array and find out this ROI's index
    ROIs = getappdata(fig,'ROIs');
    ROI.index = 1+length(ROIs)

    
    % a crosshair with point in the middle
    axes( hImageAxes );
    
    % initial pixel and delay
    % center pixel,
    % first delay slice
    ii = round(n/2);
    ij = round(m/2);
    setappdata(fig,'roiCenter',[ii,ij]);

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
    
    setappdata(fig,'ROImode','box'); % may be either box or pixel
    
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
    y = squeeze(X(ij,ii,:));
    axes( hLineoutAxes );
    lTransient = line(t,y,'Tag','lTransient','Marker','.');
    line(t,0*t,'Color','k','LineStyle','--');
    xlim([min(t),max(t)]);
    type = getappdata(fig,'stackType');
   
    if strcmp(type,'stackType');
        xlabel('z-position (\mu m)');
    else
        xlabel('Pump-probe delay (ps)');
    end
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
    ROI.lCrossXHL = lCrossXHL;
    ROI.lCrossYHL = lCrossYHL;
    ROI.lCrossCenter = lCrossCenter;
    ROI.lBox = lRoiBox;
    ROI.lBoxCorner = lRoiBoxCorner;
    ROI.lTransient = lTransient;
    ROI.lRoiBoxHL = lRoiBoxHL;
    
    ROI.QBest = []; % blank field to store biexp. fit results

    % add this ROI
    
    if isempty(ROIs)
        clear ROIs;
    end
    
    
    ROIs(ROI.index) = ROI;
    
    % save the ROIs
    setappdata(fig,'ROIs',ROIs);
    
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
    delete(ROI.lCrossXHL);
    delete(ROI.lCrossYHL);
    delete(ROI.lRoiBoxHL);
    
    
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
    uimenu(f,'Label','Histogram...','Callback',@mnuShowHistogram);
    uimenu(f,'Label','Baseline','Callback',@mnuShowBackground);
    uimenu(f,'Label','Sum-squared signal intensity image', 'Callback',@mnuSumSqIntensity);
    uimenu(f,'Label','Signal Variance Histogram...','Callback',@mnuShowSigVarHistogram);
    uimenu(f,'Label','Multiexponential fit...','Callback',@mnuMultiexpFit);
    uimenu(f,'Label','Global biexponential fit...','Callback',@mnuGlobBiexpFit);
    uimenu(f,'Label','Time const image...','Callback',@mnuTConstImage);
    uimenu(f,'Label','Time const Histogram...','Callback',@mnuTConstHist);
    uimenu(f,'Label','Principal component composite image...','Callback',@mnuPCA);
    uimenu(f,'Label','Linear decomposition...','Callback',@mnuLinearDecomp);
    uimenu(f,'Label','Principal component cluster analysis...','Callback',@mnuPCAClusters);
    uimenu(f,'Label','K-clustering analysis...','Callback',@mnuPCAClusters2);
    uimenu(f,'Label','SVD de-noise...','Callback',@mnuSVDdeNoise);
    uimenu(f,'Label','Print stack on page...','Callback',@mnuPrintStack);

    uimenu(f,'Label','Export stack','Callback',@mnuExportStack);
    uimenu(f,'Label','Export animated GIF','callback',@mnuExportGIF);
    uimenu(f,'Label','Export image stack to workspace...',...
        'callback',@mnuExportToWorkspace);
    uimenu(f,'Label','Project onto ROI','Callback',@mnuProjectOntoROI);
    uimenu(f,'Label','Project onto standard spectra...', 'Callback', @mnuProjectOntoSpectra);
    uimenu(f,'Label','Average stack','Callback',@mnuAverageStack);
    uimenu(f,'Label','PCA HetMap', 'callback', @mnuOverlapSim);
    uimenu(f,'Label','Correlation Map', 'Callback', @mnuCorrMap);
    uimenu(f,'Label','MultiExpFit (grid)','Callback',@mnuMultiExpFitGrid);
    uimenu(f,'Label','View Phasor Plot', 'Callback', @mnuViewPhasors);
    uimenu(f,'Label','Quick Visualization','Callback',@mnuQuickVis);
    
    f = uimenu('Label','EditImage');
    % IMAGE EDITING FUNCTIONS
    % anything that modifies the image data goes here
    uimenu(f,'Label','Subract background value...','Callback',@mnuBgSubtract);
    uimenu(f,'Label','Auto background','Callback',@mnuBgSubtractAuto);
    uimenu(f,'Label','Auto background each pixel','Callback',@mnuBgSubtractAutoPx);
    uimenu(f,'Label','Auto background (Reduced rank estimate)','Callback',@mnuBgSubtractReducedRank);
    uimenu(f,'Label','Normalize by ROI lineout','Callback',@normalizeByROI);
    uimenu(f,'Label','Bin Image...','Callback',@binImage);
    uimenu(f,'Label','Delay offset','Callback',@mnuDelayOffset);
    uimenu(f,'Label','Delete region','Callback',@deleteRegion);
    uimenu(f,'Label','Threshold image...','Callback',@thresholdImage);
    uimenu(f,'Label','Crop to ROI','Callback',@cropToROI);
    uimenu(f,'Label','Gaussian Bin','Callback',@mnuGaussBin);
    uimenu(f,'Label','Phasor Mask Pu/Pr: 720/810','Callback',@mnuPhasorMask720);
    uimenu(f,'Label','Phasor Mask Custom','Callback',@mnuPhsaoeMaskCustom);
    uimenu( f, 'Label', 'Pump-probe model filter', ...
        'Callback', @mnu_puprModelFilt );
    uimenu( f, 'Label','Append 3x3 neighborhood', ...
                                      'Callback', @mnuAppendNeighborhood );
    uimenu( f, 'Label','NLmeans de-noise',...
            'Callback',@mnu_nlmeans );

    
    % context menu for ROI
    cmnuROI = uicontextmenu();
    uimenu(cmnuROI,'Label','New ROI', 'Callback', @mnuNewROI);
    uimenu(cmnuROI,'Label','Delete ROI', 'Callback', @mnuRemoveROI);
    uimenu(cmnuROI,'Label','Color...','Callback',@mnuSetROIColor);
    uimenu(cmnuROI,'Label','Label...');
    setappdata(fig,'cmnuROI',cmnuROI);
    
    % set up mode menu for selecting colormaps
    puprisa_menu_mode();

end


%==========================================================================
% CALLBACK FUNCTIONS
%
function buttonDown_colorbar( src,evt,ax,  appDataTag_imageStack )
    disp('clicked on color bar!');
end


function scroll_wheel( src, evnt )
    incrementSlice(evnt.VerticalScrollCount);
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
        incrementSlice( incr );
    end
end


function incrementSlice(incr )
    ik = getappdata(gcbf,'currentSlice') + incr;
    
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
    %set(gcf,'interruptible','off')
    % Set event handling to cancel pending callback executions with new
    % events. Since the user can drag the mouse around faster than we can
    % redraw, this prevents a bunch of redraws for queueing up, and only
    % focses on the most recent event.
    set(gcf,'busyAction','cancel')
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
    set(gcf,'busyAction','cancel')

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
    
    % Restore busyAction to 'queue' so we don't miss any events
    set(gcf,'busyAction','queue')
end

%==========================================================================
% MENU CALLBACKS
%

% Main (top) menu callbaks
function mnuAutoScaleStack(src, evt)
    get(get(get(src,'Parent'),'Parent'));
    
    im = findobj(gcbf,'Tag','theImage');
    imAx = get(im, 'Parent');
    X = getappdata(gcbf, 'imageStackX');
    
    m = mean(X(:));
    s = 4*std(X(:));
    set(imAx,'clim', [-1,1]*max( [abs(m+s),abs(m-s)] ));

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
    
    X = getappdata(gcbf,'imageStackX'); % the image stack
    header = getappdata(gcbf,'header'); % header information
    t = getappdata(gcbf,'zPos');        % delays

    puprisa_componentImage( X, t, header );
    
end

function mnuMultiExpFitGrid( src, evt )
    % open the delay stack with a GUI for running multiexp fit 
    X = getappdata(gcbf,'imageStackX'); % the image stack
    t = getappdata(gcbf,'zPos');        % delays

    puprisa_exponentialFit( X, t );
end

function mnuLinearDecomp( src, evt )
    % open the delay stack with a GUI for running linear spectral
    % decomposition
    X = getappdata(gcbf,'imageStackX'); % the image stack
    t = getappdata(gcbf,'zPos');        % delays

    puprisa_linearSpectralDecomp( X, t );
end

function mnuShowHistogram( src, evt )
    % show histogram of image stack
    % display distribution of pixel intensity values
    X = getappdata( gcbf, 'imageStackX' );
    figure('Name','histogram');
    Y = sum(X.^2,3);
    hist(Y(:),1000);
    xlabel('Value');
    ylabel('Sum-squared signal');
end

function mnuShowBackground( src, evt )
    % get the image stack
    X = getappdata( gcbf, 'imageStackX' );
    
    % the z axis is the time delay
    t = getappdata(gcbf,'zPos');
    
    % find the first slice
    [~,ii] = min(t);
    
    % use that slice for background values
    bg = X(:,:,ii);
    
    % show a black-and-white image of that first slice
    figure('Name','Baseline map');
    imshow(bg);
    title('Baseline map');
    caxis([-.1,.1]);
end

function mnuSumSqIntensity( src, evt )
    % calculate sum-square intensity of signal for each pixel, display as
    % an image
    X = getappdata( gcbf, 'imageStackX' );
    figure('Name','Sum-square signal intensity');
    Y = sum(X.^2,3);
    imagesc(Y);
    set(gca,'ydir','normal');
    axis square;
    colorbar;
    colormap(jet);
end

function mnuShowSigVarHistogram( src, evt )
    % calculate histogram of delay scan signal variance
    X = getappdata( gcbf, 'imageStackX' );
    figure('Name','Signal variance histogram');
    
end

function mnuMultiexpFit(src, evt)
    % perform a multiexponential fit of the currently-displayed time series
        
    %get paramters from dialog
    fitParm = puprisa_MultiExpFitDialog();
    
    % Parameters for multi-exponential fit
    nExp  = fitParm.NoOfExp; %No. of exponentials
        
    if (nExp == 1)
        LB = [fitParm.InstRespL,fitParm.A1L,fitParm.T1L];
        UB = [fitParm.InstRespU,fitParm.A1U,fitParm.T1U];
        InCon = [fitParm.InstResponse,fitParm.A1,fitParm.T1];
    else
        %bounds for multi-exponential fit
        LB = [fitParm.InstRespL,fitParm.A1L,fitParm.T1L,fitParm.A2L,fitParm.T2L];
        UB = [fitParm.InstRespU,fitParm.A1U,fitParm.T1U,fitParm.A2U,fitParm.T2U];
        %Initial conditions for multi-exponential fit
        InCon = [fitParm.InstResponse,fitParm.A1,fitParm.T1,fitParm.A2,fitParm.T2];
    end
    
    %fit options
    options = optimset('MaxIter',12000,'TolFun',1.0e-9,'MaxFunEvals',8000);
    
       
    l = findobj(gcbf,'selected','on'); % get selected object for the line
    x = get(l,'XData');
    y = get(l,'YData');
    
    xfwhm = 0.25;%cross-correlation width in ps
    x0 = 0;%delay offset
    dx_min = min(diff(x));
    xNew = linspace(min(x) - (10*dx_min),max(x) + (10*dx_min),...
        abs(round(max(x) - min(x)/(0.1*dx_min))));%create new equally spaced time axis
    xLen = length(xNew);
    [~,ix0] = min(abs(xNew - x0));
    delFunc = 1*(xNew==xNew(ix0));
    stepFunc = 1*(xNew>=x0);
    
    I = zeros(xLen,xLen);
    for ik = 1:xLen
        % calculte impulse response at this lag 
        I(ik,:) = exp(-4*log(2)*((xNew-x0-xNew(ik))/xfwhm).^2);
        I(ik,:) = I(ik,:)/sqrt(sum(I(ik,:).^2)); %unit vector
    end
    
    if (nExp == 1)
        fn = @(A,x) SingleExpFitFunc(A,x,xNew,I,delFunc,stepFunc);
        [AA,RSS,Residual,~,~,~,J] = lsqcurvefit(fn,InCon,x,y,LB,UB,options);
        ci = nlparci(AA,Residual,'jacobian',J);%95% confidence intervals
        AAerror = (ci(:,2)-ci(:,1))/(2*1.96); %standard error
        tau = sprintf('tau = %6.3f %c %6.4f',AA(3),177,AAerror(3));
        amp = sprintf('Amp = %6.3f %c %6.4f',AA(2),177,AAerror(2));
        sprintf('%s \n%s \nInst Response = %6.4f',tau,amp,AA(1))
        ee = fn(AA,x);
    else
        fn = @(A,x) multiExpFitFunc(A,x,xNew,I,delFunc,stepFunc);
        [AA,RSS,Residual,~,~,~,J] = lsqcurvefit(fn,InCon,x,y,LB,UB,options);
        ci = nlparci(AA,Residual,'jacobian',J);%95% confidence intervals
        AAerror = (ci(:,2)-ci(:,1))/(2*1.96); %standard error
        
        tau = sprintf('tau = %6.3f %c %6.4f\t %6.3f %c %6.3f',AA(3),177,AAerror(3),AA(5),177,AAerror(5));
        amp = sprintf('Amp = %6.3f %c %6.4f\t %6.3f %c %6.3f',AA(2),177,AAerror(2),AA(4),177,AAerror(4));
        sprintf('%s \n%s \nInst Response = %6.4f\n',tau,amp,AA(1))
        ee = fn(AA, x);
    end
    line(x,ee);
   
end

function y = SingleExpFitFunc(A,t,tNew,I,delFunc,stepFunc)

yy = A(1)*delFunc...
    + (A(2)*exp(-(tNew.*stepFunc)/A(3))).*stepFunc;

yyConv = (yy*I);

y = interp1(tNew,yyConv,t,'spline');

end

function y = multiExpFitFunc(A,t,tNew,I,delFunc,stepFunc)

yy = A(1)*delFunc...
    + (A(2)*exp(-(tNew.*stepFunc)/A(3))...
    + A(4)*exp(-(tNew.*stepFunc)/A(5))).*stepFunc;

yyConv = (yy*I);

y = interp1(tNew,yyConv,t,'spline');

end


function mnuGlobBiexpFit( src, evt )
    % fits biexponential to selected ROIs,
    % using restart algorithm to find global optimum fits
    
    % prompt user for a few constratins on start points:
    % min, max time constants
    % min, max amplitudes
    % time delay limit (so as not to include time-overlap)
    % how many restarts should we allow?
    % perhaps some other convergence criteria?
    
    tMin = 0.5;     % min delay to examine
    tauMax = 100;   % max time const for initial guess, in ps
    tauMin = 0.1;   % min time const for initial guess
    AMax = 2;       % max amplitude (microvolts)  
    AMin = -2;      % min amplitude (microvolts)
    offMin = -0.1;  % min offset
    offMax = 0.1;   % max offset
    maxRestarts = 1000;
    
    % start with ROI 1 for now
    % (later make this a loop through all ROIs)
    ROIs = getappdata(gcbf,'ROIs');
    
    for iROI = 1:length(ROIs)

        % extract time axis and the decay data
        y = get(ROIs(iROI).lTransient,'YData');
        t = get(ROIs(iROI).lTransient,'XData');

        mask = t >= tMin;
        t = t(mask);
        y = y(mask);

        % variables to store results
        chi2Best = inf;    % goal is to minimize chi2
        QBest = [];         % store vector in Q
        chi2BestHistory = [];

        % function to which we fit
        % takes as search parameter Q:
        % Q = [amplitude1, amplitude2, timeconst1, timeconst2, offset]
        fn = @(Q,t) Q(1)*exp(-t./Q(3)) + Q(2)*exp(-t./Q(4)) + Q(5);

        figure;
        % loop over restart iterations
        iter = 1;
        done = 0;
        while( iter <= maxRestarts && ~done )
            % generate a new starting point
            Q0 = rand(1,5) ...
                 .* [(AMax-AMin),(AMax-AMin),...
                     (tauMax - tauMin),(tauMax - tauMin), (offMax-offMin)] ...
                 + [AMin, AMin, tauMin, tauMin, offMin];

            % find parameters Q that minimize fit error
            [Q, chi2] = lsqcurvefit(fn, Q0, t, y);

            % if this is the best fit so far, then save the result
            if( chi2 < chi2Best )
                % stop when only 1/10 percent improvement
                if( (chi2Best - chi2)/chi2Best < 0.001 )
                    done = 1;
                end

                QBest = Q;
                chi2Best = chi2;
                
                if isempty(chi2BestHistory)
                    chi2BestHistory = chi2Best;
                else
                    chi2BestHistory(end+1) = chi2Best;
                end

            end

            plot(chi2BestHistory);
            xlabel('i. best fit');
            ylabel('\chi^2');
            title('Global fit convergence');
            drawnow;

            iter = iter+1;

        end

        plot(chi2BestHistory);
        xlabel('i. best fit');
        ylabel('\chi^2');
        title('Global fit convergence');
        drawnow;

        % print the result
        figure();
        plot(t,y,t,fn(QBest,t));
        xlabel('probe delay (ps)');
        ylabel('Signal');
        legend('Data','Biexp. Fit');

        % store the results
        ROIs(iROI).QBest = QBest;
        QBest
    end
    
    % store results in appdata
    setappdata(gcbf,'ROIs',ROIs);

end

function mnuTConstImage( src, evt )
    % display with false colors based on fitted time constant
    X = getappdata(gcbf,'imageStackX');
    t = getappdata(gcbf,'zPos');
    puprisa_globalModel(t,X);
end


function mnuTConstHist( src, evt )
    
    clear all;
    prompt = {'Min probe delay (ps):','Max probe delay (ps):', 'No. exponentials'};
    dlg_title = 'Parameters for multi-exponential fit';
    num_lines = 1;
    def = {'0.098', 'Inf', '2'};
    answer = inputdlg(prompt,dlg_title,num_lines,def);

    % parse user input
    minDelay = str2num( answer{1} );
    maxDelay = str2num( answer{2} );
    nExp     = str2num( answer{3} );

    % check for valid inputs
    if( isempty(minDelay) || isempty(maxDelay) || isempty(nExp) )
        error('Invalid inputs to multiexponential fit');
    end

    % open dialog to prompt for initial guess:

    if (nExp == 1)
        prompt = {'A1:','T1 (ps):'};
        dlg_title = 'Initial conditions for exponential fit';
        num_lines = 1;
        def = {'-1', '1'};
        answer = inputdlg(prompt,dlg_title,num_lines,def);

        A1 = str2num( answer{1} );
        T1 = str2num( answer{2} );
    end

    if (nExp > 1)
        prompt = {'A1:','T1 (ps):','A2:','T2 (ps):'};
        dlg_title = 'Initial conditions for multi-exponential fit';
        num_lines = 1;
        def = {'-1', '1','-0.5', '5'};
        answer = inputdlg(prompt,dlg_title,num_lines,def);

        A1 = str2num( answer{1} );
        T1 = str2num( answer{2} );
        A2 = str2num( answer{3} );
        T2 = str2num( answer{4} );
    end


    Y = getappdata(gcf,'imageStackX');
    x = getappdata(gcf,'zPos');

    mask = ((x >= minDelay) & (x <= maxDelay));
    x = x(mask);

    [m n ~] = size(Y);
    [dummy sizeOFx] = size(x);
    index = 1;

    %fit options
    options = optimset('MaxIter',1000,'TolFun',1.0e-8,'MaxFunEvals',1000);

    for i2 = 1:m
        for j = 1:n
            y = (squeeze(Y(i2,j,:)))';
            y = y(mask);
            if (nExp == 1)
                sum = [0 0];
                fn = @(A,xx) A(1)*exp(-xx./A(2))+ A(3);
                AA = lsqcurvefit(fn, [A1,T1,0.01], x, y,[],[],options);
                if (AA(1) < 0) %select only fits with negative signal
                    sum(1) = sum(1)+AA(2);
                    sum(2) = sum(2)+1;
                end
                if (sum(2) > 0 )
                    tauShort(index) = sum(1)/sum(2);
                    index= index+1;
                end
            else
                fn = @(A,xx) A(1)*exp(-xx./A(2)) + A(3)*exp(-xx./A(4))+A(5);
                [minI ~] = min(y);
                if (minI < -0.1) %select only specific pixels with good Signal
                    [AA,Residual] = lsqcurvefit(fn, [A1,T1,A2,T2,0.01], x, y,[],[],options);
                    if (AA(1) < 0 && AA(3) < 0) %select only fits with negative signal
                        if(AA(2) < AA(4))%assign to short or long component group
                            tauLong(index) = AA(4);
                            tauShort(index) = AA(2);
                            ampLong(index) = AA(3);
                            ampShort(index) = AA(1);
                        else
                            tauLong(index) = AA(2);
                            tauShort(index) = AA(4);
                            ampLong(index) = AA(1);
                            ampShort(index) = AA(3);
                        end
                        fitError(index) = sqrt(Residual/sizeOFx);
                        index= index+1;
                    end
                end
            end
        end
    end

    name = getappdata(gcf,'fileName');
    FileName_S = [name(1:(end-4)),'_tauShort.txt'];
    FileName_T = [name(1:(end-4)),'_tauLong.txt'];
    FileName_AS = [name(1:(end-4)),'_ampShort.txt'];
    FileName_AT = [name(1:(end-4)),'_ampLong.txt'];
    FileName_FER = [name(1:(end-4)),'_fitError.txt'];
    %write to file
    fid = fopen(FileName_S, 'w');
    fprintf(fid, '%12.6f\n', tauShort);
    fclose(fid);
    fid = fopen(FileName_T, 'w');
    fprintf(fid, '%12.6f\n', tauLong);
    fclose(fid);
    fid = fopen(FileName_AS, 'w');
    fprintf(fid, '%12.6f\n', ampShort);
    fclose(fid);
    fid = fopen(FileName_AT, 'w');
    fprintf(fid, '%12.6f\n', ampLong);
    fclose(fid);
    fid = fopen(FileName_FER, 'w');
    fprintf(fid, '%12.6f\n', fitError);
    fclose(fid);

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
    
    % refresh ROI displays
    updateCursor();
    
end

function thresholdImage( src, evt )
    % ask for threshold level
    prompt = {'Threshold level (sum-squared delay trace value):'};
    dlg_title = 'Zero out pixels below threshold';
    num_lines = 1;
    def = {'0.00'};
    answer = inputdlg(prompt,dlg_title,num_lines,def);
    threshold = str2num(answer{1});
    if isempty(threshold)
        error('invalid entry');
    end
    
    % get the image stack and apply thresholding mask
    X = getappdata( gcbf, 'imageStackX' );
    [nx,ny,nt] = size(X);
    M = repmat(sum(X.^2,3) >= threshold, [1,1,nt]);
    X = X.*M;
    setappdata( gcbf, 'imageStackX', X );
    
    % refresh the image
    im = findobj(gcbf,'Tag','theImage');
    ik = getappdata(gcf,'currentSlice');
    set(im,'CData', X(:,:,ik));
    drawnow;
    
    % refresh ROI displays
    updateCursor();
end

function cropToROI( src, evt )
    zPos = getappdata(gcbf,'zPos');
    fileName = getappdata(gcbf,'fileName');
    header = getappdata(gcbf,'header');
    chan = getappdata(gcbf,'chan');
    X = getappdata(gcbf,'imageStackX');
    
    % locate ROI
    ROIs = getappdata(gcbf,'ROIs');
    center = ROIs(1).center;
    box = ROIs(1).box;
    
    imin = center(1)-round(box(1)/2);
    imax = center(1)+round(box(1)/2);
    jmin = center(2)-round(box(2)/2);
    jmax = center(2)+round(box(2)/2);
    
    
    XX = X(  jmin:jmax,imin:imax, : );

    puprisa_viewChannel( XX, zPos, fileName, header, chan );
end

function mnu_puprModelFilt( src, evt )
    % filter image stack by fitting it to a pump-probe model
    
    % recall preferences
    tfwhm = getPrefWithDefault( 'modelFit_tfwhm', 0.25 );
    tau_min = getPrefWithDefault( 'modelFit_tau_min', 0.1 );
    tau_max = getPrefWithDefault( 'modelFit_tau_max', 30 );
    
    % prompt user for options
    prompt = { 'Instrument response FWHM (ps)', ...
                    'Min. decay time (ps)', ...
                    'Max decay time (ps)' };
    a = inputdlg( prompt, 'Model parameters', ...
                    [1,1,1], ...
                    { num2str(tfwhm), ...
                      num2str(tau_min), ...
                      num2str(tau_max) } );
    
    % parse and validate response
    tfwhm = validate_doublePositiveNonzero( a{1}, prompt{1} );
    tau_min = validate_doublePositiveNonzero( a{2}, prompt{2} );
    tau_max = validate_doublePositiveNonzero( a{3}, prompt{3} );
    
    % generate model and fit
    imageStack = getappdata( gcbf, 'imageStackX' );
    t = getappdata( gcbf, 'zPos' );
    
    [nr, nc, nt] = size( imageStack );
    delayScans = reshape( imageStack, [nr*nc, nt] );
    [Yfit, err] = puprModelFit( delayScans, t, tfwhm, ...
        tau_min, tau_max );
    imageStack_fit = reshape( Yfit, [nr, nc, nt] );
    setappdata( gcbf, 'imageStackX', imageStack_fit );
    
    
    % refresh the image
    im = findobj(gcbf,'Tag','theImage');
    ik = getappdata(gcbf,'currentSlice');
    set(im,'CData', imageStack_fit(:,:,ik));
    drawnow;
    
    % refresh ROI displays
    updateCursor();
    
    % save preferences
    setPref( 'modelFit_tfwhm', tfwhm );
    setPref( 'modelFit_tau_min', tau_min );
    setPref( 'modelFit_tau_max', tau_max );
end

function val = validate_doublePositiveNonzero( str, label )
    % converts string to a double, and throws an error if it is 
    % not > 0. the input 'label' is used to make a more informative 
    % error message
    
    val = str2num( str );
    
    if isempty( val )
        error( [label, ' is not a valid number.'] );
    end
    
    if val < 0
        error( [label, ' must be a nonzero, positive number.'] );
    end
end
        

function val = getPrefWithDefault( prefname, default )
    % get preference, optionally specifying a default value
    if ispref( 'puprisa', prefname )
        val = getpref( 'puprisa', prefname );
    else
        val = default;
    end
end

function setPref( prefname, val )
    % set preference, making a new entry in preferences structure if
    % necessary
    if ispref( 'puprisa', prefname )
        setpref( 'puprisa', prefname, val );
    else
        addpref( 'puprisa', prefname, val );
    end
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

function mnuBgSubtractAutoPx( src, evt )
    % get the image stack
    X = getappdata( gcbf, 'imageStackX' );
    
    % use first pixels to subtract background
    t = getappdata(gcbf,'zPos');
    
    
    % find slice with most negative delay
    [~,ii] = min(t);
    
    % use that slice for background values
    bg = X(:,:,ii);
   
    % repeat the background values to match the size of the entire image
    % stack
    BG = repmat(bg, [1,1,length(t)]);
    
    % subtract background values
    X = X - BG;
    
    % store the result
    setappdata( gcbf, 'imageStackX', X );
    
    % refresh the image
    im = findobj(gcbf,'Tag','theImage');
    ik = getappdata(gcf,'currentSlice');
    set(im,'CData', X(:,:,ik));
    drawnow;
    
    % refresh ROI displays
    updateCursor();
   
    
end

function mnuBgSubtractReducedRank( src, evt )
% use reduced-rank approximation to estimate background level for each
% pixel

    % get the image stack
    X = getappdata( gcbf, 'imageStackX' );
    
    % use first pixels to subtract background
    t = getappdata(gcbf,'zPos');
    
    [X,~] = puprisa_baselineCorrection( t, X, 6 );

    % store the result
    setappdata( gcbf, 'imageStackX', X );
    
    % refresh the image
    im = findobj(gcbf,'Tag','theImage');
    ik = getappdata(gcf,'currentSlice');
    set(im,'CData', X(:,:,ik));
    drawnow;
    
    % refresh ROI displays
    updateCursor();
end

% ROI context menu callbacks
function mnuNewROI(src, evt )
    newROI(gcbf, findobj(gcbf,'Tag','imageAxes'),...
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

function mnuDelayOffset(src, evt)
    % subtract arbitrary offset from delay scans
    
    % get delay vector and current slice
    delayPos = getappdata(gcbf,'zPos');
    currentSlice = getappdata(gcbf,'currentSlice');
    currentDelay = delayPos(currentSlice);
    
    % ask for offset
    prompt = {'Delay Offset:'};
    dlg_title = 'Delay Offset (ps)';
    num_lines = 1;
    def = {num2str(currentDelay)}; % default is current delay
    answer = inputdlg(prompt,dlg_title,num_lines,def);
    
    if ~isempty(answer)

        delay0 = str2num(answer{1});
        if isempty(delay0)
            error('invalid entry');
        end

        % apply the offst
        delayPos = delayPos - delay0;
        setappdata(gcbf, 'zPos', delayPos);

        % update x-axis of trend plot
        incrementSlice(0);
        
         % update x-axis of trend plot
        lTransient = findobj(gcbf,'Tag','lTransient');
        set(lTransient,'XData',delayPos);
        
        axTransient = findobj(gcbf,'Tag','axTransient');
        set(axTransient,'XLim',[min(delayPos) max(delayPos)]);

    end
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
    
    % draw a scalebar
    header = getappdata(gcbf,'header');
    
    if ~isempty(header)
        puprisa_initScaleBar( gca, header );
    end
end

function mnuProjectOntoROI( obj, evt )
    % project image onto delay trace given by the ROI
    
    % retrieve info on ROI
    hROIs = getappdata(gcbf, 'ROIs');
    
    if( length(hROIs) > 1 )
        error('Only one ROI allowed for projection. Delete unneeded ROIs first.');
    end
    
    % get curve for the ROI
    h = hROIs(1).lTransient;
    
    y = get(h,'YData');
    
    % now project image stack onto y
    X = getappdata(gcf,'imageStackX');  % get the image stack
    [nx,ny,~] = size(X);
    
    % duplicate y
    Y = shiftdim(repmat(y,[ny,1,nx]),2);
    
    % multiply each delay scan against y
    XY = X.*Y;
    
    % sum the overlap with y
    Xprj = sum(XY,3);
    
    figure('name','ROI projection');
    imagesc(Xprj);
    axis image;
    set(gca,'ydir','normal');
    colormap(gray);
    
end

function mnuProjectOntoSpectra( src, evt )
    % Warning: this is still highly sensitive to linear offset. Will need
    % to change from linear projection to linear least-squares fitting with
    % an offset.

    % project image onto delay trace loaded from a file
    puprisaDir = which('puprisa.m');
    standardsDir = [puprisaDir(1:end-9),'..',filesep,'standards'];
    %standardsDir = 'C:\Users\jw295\projects\principalComponent_pumpProbeImages\pumpProbeStandards';
    % prompt user for standard
    [FileName,PathName,FilterIndex] = ...
        uigetfile('*.mat','Select pump-probe standards',[standardsDir,filesep,'*.mat']);
    
    dat = load([PathName,filesep,FileName]);
    
    t = getappdata(gcbf,'zPos');
    
    % Determine whether we need to interpolate the standard reference
    % spectra by comparing the reference spectrum's time delay points
    % to the image stack.
    doInterp = 1;
    if( length(t) == length(dat.t) )    % first check lengths
        if( ~any( t - dat.t ) )         % then check individual values
            
            % don't interpolate if everything matches up perfectly
            doInterp = 0;               
        end
    end
    
    if doInterp
        % Perform linear interpolation. Leave extrapolated values
        % at zero; this has the effect of disregarding time delays from 
        % the image stack that were not sampled in the reference spectrum.
        
        y_eu = interp1( dat.t, dat.eumelanin, t, 'linear', 0 );
        y_pheo = interp1( dat.t, dat.pheomelanin, t, 'linear', 0 );
    else
        % Otherwise, don't interpolate; use raw data from reference
        % spectra. Convert row vectors to column vectors by transponsing
        % to get everything to line up.
        
        y_eu = dat.eumelanin.';
        y_pheo = dat.pheomelanin.';
    end
    
    % now project the image stack onto y
    X = getappdata(gcf,'imageStackX');  % get the image stack
    [nx,ny,~] = size(X);
    
    % duplicate y
    Y_eu = shiftdim(repmat(y_eu,[ny,1,nx]),2);
    Y_pheo = shiftdim(repmat(y_pheo,[ny,1,nx]),2);

    
    % multiply each delay scan against y
    XY_eu = X.*Y_eu;
    XY_pheo = X.*Y_pheo;
    
    % sum the overlap with y
    Xprj_eu = sum(XY_eu,3);
    Xprj_eu( Xprj_eu < 0 ) = 0;
    
    Xprj_pheo = sum(XY_pheo,3);
    Xprj_pheo( Xprj_pheo < 0 ) = 0;
    
    mag = max(max( [Xprj_eu(:), Xprj_pheo(:)] ));
    
    figure('name','Spectroscopic projection');
    
    XprjRGB = zeros( nx, ny, 3);
    XprjRGB(:,:,1) = Xprj_eu / mag;
    XprjRGB(:,:,2) = Xprj_pheo / mag;
    
    imshow(XprjRGB*2);
    set(gca,'ydir','normal');
    
    euFrac = sum(Xprj_eu(:)) / sum(Xprj_pheo(:)+Xprj_eu(:));
    
    disp(['Eumelanin fraction = ', num2str(euFrac*100), '%']);
end

function binImage(obj, evt)
    % bin image, lowering resolution to improve SNR
    binFactor = 2;
    
    X = getappdata(gcbf,'imageStackX');
    [nr, nc, ns] = size(X);

    Xnew = zeros(nr/binFactor, nc/binFactor, ns);
    
    h = waitbar(0,'Binning');
    % for each slice
    for is = 1:ns
        for ir = 1:(nr/binFactor)
            for ic = 1:(nc/binFactor)
                selected = X( ((1:2)+((ir-1)*binFactor)),...
                    ((1:2)+((ic-1)*binFactor)),...
                    is );
                Xnew(ir, ic, is) = ...
                    mean(mean(selected));
            end
        end
        waitbar( is/ns, h );
        drawnow;
    end
    
    close(h);
    
    
    zPos = getappdata(gcbf,'zPos');
    header = getappdata(gcbf,'header');
    fileName = getappdata(gcbf,'fileName');
    puprisa_viewChannel( Xnew, zPos, fileName, header )

    
end

function deleteRegion( obj, evt )
    % prompt user for a region to delete
    % useful before prin comp. analysis to get rid of surgical ink, etc.
    
    % set current axes to the displayed image
    h = findobj(gcf,'tag','imageAxes');
    hImg = findobj(h,'tag','theImage');
    axes(h);                    
    
    % get a region of interest
    hROI = impoly;              % then make a new ROI
    
    m = createMask(hROI,hImg);
    
    % get the image stack
    X = getappdata(gcbf,'imageStackX');
    [nx,ny,nz] = size(X);
    
    M = repmat(m, [1,1,nz]);
    
    X = X.*(1-M);
    
    setappdata(gcbf,'imageStackX', X);
    
    delete(hROI);
    % redraw image
    redrawImage();
    
    updateCursor();
end

function normalizeByROI(obj, evt)
    r=getappdata(gcbf,'ROIs');
    y = get( r(1).lTransient, 'YData' );
    
    X = getappdata(gcbf,'imageStackX');
    %X(X < -0.1) = 0;
    %Y(1,1,:) = y;
    
    
    [nRows,nCols,nSlices] = size(X);
    
    y = (squeeze(sum(sum(abs(X),1),2))) / (nRows*nCols);
    Y(1,1,:) = y;
    
    Y = repmat(Y,[nRows, nCols,1]);
    
    X = X ./ Y;
    
    setappdata(gcbf,'imageStackX', X);
   
    % redraw image
    redrawImage();
end

function redrawImage()
    ik = getappdata(gcbf,'currentSlice');
    
    h = findobj(gcbf,'Tag','imageAxes');
    
    X = getappdata(gcbf,'imageStackX');

    im = findobj(gcbf,'Tag','theImage');
    set(im,'CData', X(:,:,ik));
    
    drawnow;
end

function mnuOverlapSim( src, evt )
    X = getappdata(gcbf,'imageStackX');
    
    puprisa_PCA_hetMap( X );
    
   
end

function mnuCorrMap( src, evt )
% Time signal correlation map.
% Highlights pixels which are highly correlated in time.
% could be used to remove noise or separate signals that vary randomly
% in time from signals that are correlated.
%
% See Jonathan, Enfield, and Leahy. J. Biophotonics. v4 no 9, p 583-587
% (2011).
    X = getappdata(gcbf,'imageStackX');
    
    [nx, ny, nt] = size(X);
    
    % no. pixels in neighborhood
    k = 7;   % must be an odd number
    
    C = zeros(nx, ny);
    
    % calculate correlation image
    wb = waitbar(0,'Calculating correlation image...')
    for ix = 1:nx
        for iy = 1:ny
            sig = squeeze(X(ix,iy,:));
            C(ix,iy) = sum(xcorr(sig, sig,'coeff'));
        end
        waitbar(ix/nx,wb);
    end
    close(wb);
    
    figure
    imagesc(C);
    
    
end

function mnuGaussBin( src, evt )

 % get the image stack
    X = getappdata( gcbf, 'imageStackX' );
    
    % use first pixels to subtract background
%     t = getappdata(gcbf,'zPos');
  
%     h_filt = fspecial('gaussian',[3 3],1);%default are [3 3] square matrix with stdev 0.5
%     for ii = 1:length(t)
%         X(:,:,ii) = imfilter(X(:,:,ii),h_filt,'replicate');
%     end 
      
     X = smooth3(X,'gaussian',[3 3 1],1);

    % store the result
    setappdata( gcbf, 'imageStackX', X );
    
    % refresh the image
    im = findobj(gcbf,'Tag','theImage');
    ik = getappdata(gcbf,'currentSlice');
    set(im,'CData', X(:,:,ik));
    drawnow;
    
    % refresh ROI displays
    updateCursor();
    
    display('binned spatial pixels using Gaussian filter')
end
    

function ComputePhasors()
 
 %      f = input('input frequency (in THz):  ');
      
        f = 0.25;
      
    display(['Phaosrs calculated at freq: ', num2str(f)]);

    % get the image stack
    X = getappdata( gcbf, 'imageStackX' );
    
    % use first pixels to subtract background
     t = getappdata(gcbf,'zPos');
     
    [nr,nc,tcomp] = size(X);

    omega = 2*pi*f;

    X = reshape(X,[nr*nc,tcomp]);

    Xint = (sum(abs(X),2));

    g = sum(X.*repmat(cos(omega*t),[nr*nc,1]),2)./Xint;
    s = sum(X.*repmat(sin(omega*t),[nr*nc,1]),2)./Xint;

    g(Xint == 0) = 0;
    s(Xint == 0) = 0;

            NN = 1;
            X2sum = sum(abs(X),2);
            X2mean = mean(X2sum(:));
            X2stdev = std(X2sum(:));
            IntImMask = (X2sum > (X2mean + NN*X2stdev));

            x=linspace(-1,1,256);
            ctrs{1} = x;
            ctrs{2} = x;
            [N,~] = hist3([g(IntImMask>0),s(IntImMask>0)],ctrs);

     
    setappdata(gcbf,'g',g);
    setappdata(gcbf,'s',s);
    setappdata(gcbf,'NPhasorCounts',N);
    
 end % ComputePhasors()
 
function mnuViewPhasors( src, evt)

    ComputePhasors()
    
    N = getappdata(gcbf,'NPhasorCounts');
    x = linspace(-1,1,256);

    figure('Name', 'Phasor Plot');
    title(['Histogram using method: threshold above 1 StDev'])
    imagesc('xdata',x,'ydata',x,'cdata',N.'); 
    load 'contourColorMap'
    colormap(contourColorMap);
    xlim([-1,1]);ylim([-1,1]);
    xlabel('g(\omega)'); ylabel('s(\omega)');
    xx = linspace(0,1,512);
    semicircle = sqrt(xx.*(1-xx));
    line(xx,semicircle,'LineStyle','--','Color','k');
    line(-xx,-semicircle,'LineStyle','--','Color','k');

    for drawlines = -.8:.2:.8;
        line(drawlines,-1:.05:1,'Color','k','Linewidth',1);
        line(-1:.05:1,drawlines,'Color','k','Linewidth',1);
    end
       
    
end

function mnuQuickVis( src, evt )
    % provide quick, model-free visualization
    t = getappdata( gcbf, 'zPos' );
    imageStack = getappdata( gcbf, 'imageStackX' );
    
    IRGB = puprisa_quickVis( imageStack, t );
    
    sc = mean(IRGB(:)) + 4*std(IRGB(:));
    IRGB2 = IRGB / sc;
    
    % clip saturated pixels
    %
    % The naive way to do this is:
    % 
    IRGB2( IRGB2 > 1 ) = 1;
    % 
    % But that changes the color of saturated pixels, and they tend to
    % whiten.
    % 
    % The procedure below preserves color
    
%     [nr,nc,~] = size(IRGB2);
%     lRGB = reshape(IRGB2,[nr*nc,3]);
%     
%     for satChan = 1:3
%     
%         % pick out pixel that are saturated in red
%         rsat = lRGB(:,satChan) .* (lRGB(:,satChan) > 1) + 1*(lRGB(:,satChan) <= 1);
%         % and rescale those
%         lRGB(:,satChan) = lRGB(:,satChan) ./ rsat;
%         lRGB(:,satChan) = lRGB(:,satChan) ./ rsat;
%         lRGB(:,satChan) = lRGB(:,satChan) ./ rsat;
%     end

    
    figure;
    imshow( IRGB2 );
    set(gca,'ydir','normal');
    
end

function GenPhasorMask()

    % get the image stack
    X = getappdata( gcbf, 'imageStackX' );    
    % use first pixels to subtract background
    % t = getappdata(gcbf,'zPos');
    
    nt = size(X,3);
    
    g = getappdata(gcbf,'g');
    s = getappdata(gcbf,'s');

    PhasorMask = getappdata(gcbf,'PhasorMask');
    
    [py,px] = size(PhasorMask);
    npx = sqrt(size(g,1));% g is nr*nc,  and nr=nc
    G = reshape(g,[npx, npx]);
    S = reshape(s,[npx, npx]);

    ImMask = zeros(npx,npx);
    for pxx = 1:npx;
        for pyy = 1:npx;
            gp = round((G(pyy,pxx) * px/2) + px/2) ; %maps from phasor val to pixel val of phasor plot
            sp = round((S(pyy,pxx) * py/2) + py/2) ;
            if gp==0; gp=1;elseif sp==0; sp=1;end
            ImMask(pyy,pxx) = PhasorMask( sp, gp);

        end
    end


    X = X.*repmat(1.0.*ImMask,[1,1,nt]);    
    
    % store the result
    setappdata( gcbf, 'imageStackX', X );
    
    % refresh the image
    im = findobj(gcbf,'Tag','theImage');
    ik = getappdata(gcbf,'currentSlice');
    set(im,'CData', squeeze(X(:,:,ik)));
    drawnow;
    
    % refresh ROI displays
    updateCursor();
    
    display('Masked images using Phasor mask at Pu/Pr: 720/810')
    
end

function mnuPhasorMask720( src, evt )
    load PhasorMasks\SkinEuPheoPhasorROI.mat
    setappdata(gcbf,'PhasorMask',PhasorMask);

    mnuViewPhasors( src, evt)
    line(xi,yi,'LineWidth',3,'color','r');
    
    GenPhasorMask       
    
    display('Masked images using Phasor mask at Pu/Pr: 720/810')    
end

function mnuPhsaoeMaskCustom( src, evt )

    mnuViewPhasors( src, evt)
    [PhasorMask,xi,yi] = roipoly; 
    line(xi,yi,'LineWidth',3,'color','r');
    
    setappdata(gcbf,'PhasorMask',PhasorMask);

    GenPhasorMask       
    
    display('Masked images using custom Phasor Mask')    
end


function mnuAppendNeighborhood( src, evt )
% append 3x3 neighborhood to the image stack
    
    % busy cursor
    set(gcbf,'Pointer','watch');
    drawnow;

    % get the image stack
    % convert to single-precision to minimize memory usage
    X = single(getappdata( gcbf, 'imageStackX' ));
    
    % append neighborhood
    Xn = puprisa_appendNeighborhood(X);
    
    % modify delay vector to accomodate larger stack
    t = getappdata(gcbf,'zPos');
    t = repmat(t,[1,9]);
    
    % open the modified stack and delay position vector
    header = getappdata(gcbf,'header');
    fileName = getappdata(gcbf,'fileName');
    puprisa_viewChannel( Xn, t, fileName, header, 1 );
    
    % change cursor back
    set(gcbf,'Pointer','arrow');
    
end


% nonlocal means de-noise
function mnu_nlmeans( src, evt )
    disp('Nonlocal means de-noising...');
    disp('MATLAB implementation by Dirk-Jan Kroon');
    disp('See https://www.mathworks.com/matlabcentral/fileexchange/27395-fast-non-local-means-1d--2d-color-and-3d');
    disp('Be sure to cite the following paper:')
    disp('A. Buades, B. Coll, and J.-M. Morel, “A non-local algorithm for');
    disp('    image denoising,” in Computer Vision and Pattern Recognition,');
    disp('    2005, vol. 2, pp. 60–65.');
    disp('');
    disp('(this may take a while...)')
    
    

    % check for nonlocal means install
    if isempty( which('NLMF') )
        error('NLMF not installed. Download from MATLAB central, and ensure NLMF.m is in the path.');
    end
    
    % get the image stack and convert to single-precision ranging from 0--1
    X = getappdata( gcbf, 'imageStackX' );
    Xmin = min(X(:));
    Xmax = max(X(:));
    Xs = single( (X - Xmin)/(Xmax - Xmin) );
    
    % de-noise
    Xdn = NLMF(Xs);
    
    % re-scale
    Xdn = single((Xdn * (Xmax - Xmin)) + Xmin);
    
    % store the result
    setappdata( gcbf, 'imageStackX', Xdn );
    
    % refresh the image
    im = findobj(gcbf,'Tag','theImage');
    ik = getappdata(gcf,'currentSlice');
    set(im,'CData', X(:,:,ik));
    drawnow;
    
    % refresh ROI displays
    updateCursor();
    
end