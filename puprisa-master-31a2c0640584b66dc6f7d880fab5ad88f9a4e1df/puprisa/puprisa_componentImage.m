% MAYBE: extend to more than 3PCS, give option for primary color segmenting
%       (classify pixel by strongest fit)
% maybe give option of scaling the PCs individually
% hint: use gcbf to interact with calling figure

function puprisa_componentImage( X, t, header )
    % generate a princpal component image
    % overlay components with selected colors for each component
    % cast to single-precision float for speed
    
    f = initFigure(single(X), t, header);
    initMenus(f);
    
end   

%==========================================================================
% INITIALIZATION ROUTINES
%
function initOptions(fig)
    pcaOptions.scaleByEigenvalues = 0;      % scale eigenvectors by 1/eigval
    pcaOptions.subtractChannelMeans = 0;    % subtract mean of each channel (delay) to get a 'proper' PC decomposition
    setappdata(fig,'pcaOptions',pcaOptions);
    
    setappdata(fig,'useBipolarColoring',false);
    
    setappdata(fig,'gamma',1.0);
    
end

 function fig = initFigure( X, t, header )
    fig = figure('Name','Component Image','Units','pixels',...
        'ResizeFcn',@resizeCallback);
    
    initOptions(fig);

    
    layoutData = initLayout( fig );
    
    % lay out axes locations
    layout(fig);
    
    % store data
    setappdata(fig,'header',header);    % file header
    setappdata(fig,'ImageStack',X);     % the image stack
    setappdata(fig, 't', t );           % probe delays
      
    colors = [1,0,0; 0,1,0; 0,0,1]; % RGB
    setappdata(fig,'PCColors',colors);
    whichPCs = [1,2,3]; % which principal components to use
    setappdata(fig,'whichPCs', whichPCs);
    signs = [1,1,1];
    setappdata(fig,'signs',signs);
    enabledPCs = [1,1,1];   % which PC channels are enabled
    setappdata(fig,'enabledPCs',enabledPCs);
    
    [V,S] = generatePCsFromImageStack( X, fig );
    setappdata(fig,'eigenvectors',V);
     
    % generate RGB composite image from singular values
    PCImageRGB = generatePCImageRGB( X, V, whichPCs, colors, signs,  0*whichPCs + 1, fig );
    
    setappdata(fig, 'PCImageRGB', PCImageRGB);
    
    renormalize(fig);
    

 end
 
 function renormalize( fig )
    PCImageRGB = getappdata(fig,'PCImageRGB');
    layoutData  = getappdata( fig, 'layoutData');
    header = getappdata(fig,'header');
    
    % normalize
    PCDisplayImageRGB = PCImageRGB/max(PCImageRGB(:));
    PCDisplayImageRGB(PCDisplayImageRGB>1) = 1;
    %PCImageRGB(:,:,1) = PCImageRGB(:,:,1) / max(max(PCImageRGB(:,:,1)));
    %PCImageRGB(:,:,2) = PCImageRGB(:,:,2) / max(max(PCImageRGB(:,:,2)));
    %PCImageRGB(:,:,3) = PCImageRGB(:,:,3) / max(max(PCImageRGB(:,:,3)));
    
    % image plot
    axes(layoutData.imageAxes);
    im = image(0*PCDisplayImageRGB);
    set(gca,'Tag','imageAxes');
    set(im,'tag','componentImage');
    set(gca,'YDir','normal');
    axis square;
    
    % put in a scale bar
    if ~isempty(header)
        puprisa_initScaleBar( gca, header );
    end

    
    % histograms
    axes(layoutData.histogramAxes);
    nBins = 512;
    
    [histR,histG,histB, bins] = makeHistRGB(PCImageRGB, nBins);
    
    % throw away the zero
    histMax = max( [histR(2:end), histG(2:end),histB(2:end)] );
    
    histAlpha = 0.5;
    patch( [bins(1),bins,bins(end)], ....
        [0,histR,0], 'red', 'EdgeColor','none','FaceAlpha',histAlpha,...
        'Tag','histRPatch');
    patch( [bins(1),bins,bins(end)], ....
        [0,histG,0], 'green', 'EdgeColor','none','FaceAlpha',histAlpha,...
        'Tag','histGPatch' );
    patch( [bins(1),bins,bins(end)], ....
        [0,histB,0], 'blue', 'EdgeColor','none','FaceAlpha',histAlpha,...
        'Tag','histBPatch' );
    
    xlim([bins(2),bins(end)]);
    
    setAutoColorScale( fig );
    colorScale = getappdata(fig, 'colorScale');
    line([0,0]+colorScale,[0,histMax],'Color','magenta','LineWidth',3,...
        'Tag','histColorScaleLine','ButtonDownFcn',@buttonDown_histColorScaleLine);
    
    setappdata(gcf, 'colorScale',colorScale);
    
    updatePCImage( fig );
    update( fig );
 end

 function [histR, histG, histB, bins] = makeHistRGB( imageRGB, nBins )
    
    bins = linspace(min(imageRGB(:)), max(imageRGB(:)), nBins );
    
    R = imageRGB(:,:,1);
    G = imageRGB(:,:,2);
    B = imageRGB(:,:,3);
    
    histR = hist( R(:), bins )/length(R(:))*100;
    histG = hist( G(:), bins )/length(G(:))*100;
    histB = hist( B(:), bins )/length(B(:))*100;
    
 end
 
 % menus
 function initMenus(f)
    h = uimenu(f,'Label','PrincipalComponents');
    uimenu(h, 'Label', 'Compute from image SVD', 'Callback',@mnuComputePCFromSVD );
    uimenu(h, 'Label', 'Load Principal Components...', 'Callback',@mnuLoadPCs );
    uimenu(h, 'Label', 'Compute from ROIs', 'Callback',@mnuPCsFromROI );
    uimenu(h, 'Label', 'Switch PC Colors 1 and 2', 'Callback',@mnuSwitchPCColor1and2);
    uimenu(h, 'Label', 'Switch PC Colors 2 and 3', 'Callback',@mnuSwitchPCColor2and3);
    
    h = uimenu(f,'Label','Options');
    uimenu(h,'Label','Scale Components by Eigenvalues',...
        'Callback',@mnuOptScaleByEigval );
    uimenu(h,'Label','Subtract Channel Means',...
        'Callback',@mnuOptSubtractChanMeans )
    uimenu(h,'Label','Change Scale Bar Size','Callback',@mnuChangeScaleBar);
    uimenu(h,'Label','Recolor to Eumelanin Percentage','Callback',@mnuEumelaninPercentage);
    uimenu(h,'Label','Show Fourier Transform','Callback',@mnuFourierTransform);
    uimenu(h,'Label','Show Bipolar Components','Callback',@mnuBipolarComp);
end
 
 % GUI layout initialization
 % creation of all UI elements, callbacks
 function layoutData = initLayout( fig )
    % set up axes and areas
    layoutData.imageAxes = axes('Units','pixels','Tag','imageAxes');
    layoutData.histogramAxes = axes('Units','pixels','Tag','histogramAxes');
    
    layoutData.nPCsDisplayed = 3;      % no. PCs displayed on right
    for ii = 1 : layoutData.nPCsDisplayed
        layoutData.PCDisplay(ii).PCAxes = axes('Units','pixels','Tag','PCAxes');
        layoutData.PCDisplay(ii).enabledCheckbox ...
            = uicontrol('Style','checkbox','String','Enabled',...
                        'Callback',@enabledCheckboxCallback,...
                        'UserData',ii, 'Value', 1);
        layoutData.PCDisplay(ii).invertCheckbox ...
            = uicontrol('Style','checkbox','String','Invert',...
                        'Callback',@invertCheckboxCallback,...
                        'UserData',ii, 'Value', 0);
    end
    
    layoutData.rightColumnWidth = 256;
    layoutData.bottomRowHeight = 128;
    
    layoutData.PCAxesHeight = 128;
    
    layoutData.edGamma = uicontrol('Style','edit','String','1.0','Callback',@gammaChanged);

    setappdata( fig, 'layoutData', layoutData );
    
end
    
% GUI layout function
% places UI elements where they belong
% this is called by the ResizeFcn callback whenever the window is resized
function layout( fig )
    % retrieve structure containing layout information
    layoutData = getappdata( fig, 'layoutData' );
    
    imageAxes = findobj(fig,'Tag','imageAxes');
    histogramAxes = findobj(fig,'Tag','histogramAxes');
    
    windowPos = get(fig,'Position');
    windowWidth = windowPos(3);
    windowHeight = windowPos(4);
    
    % Should we enforce some minimum size?
    
    % set position of image axes
    set(imageAxes,'Position',...
        [0,...
         layoutData.bottomRowHeight,...
         windowWidth-layoutData.rightColumnWidth,...
         windowHeight-layoutData.bottomRowHeight]);
     
    % set position of histogram axes
    set(histogramAxes,'Position',...
        [0,...
         0,...
         windowWidth-layoutData.rightColumnWidth,...
         layoutData.bottomRowHeight]);
     
     % set position of principal component axes
     for ii = 1 : layoutData.nPCsDisplayed
         % fill in from the top down
         set( layoutData.PCDisplay( ii ).PCAxes, 'Position', ...
             [ windowWidth - layoutData.rightColumnWidth, ...
               windowHeight - ii*layoutData.PCAxesHeight, ...
               layoutData.rightColumnWidth,...
               layoutData.PCAxesHeight ] );
         
         % enabled checkbox
         pos = get( layoutData.PCDisplay(ii).enabledCheckbox, 'Position');
           
         set( layoutData.PCDisplay(ii).enabledCheckbox, 'Position',...
            [ windowWidth - pos(3), ...
              windowHeight - ii*layoutData.PCAxesHeight, ...
              pos(3), pos(4)] );
          
         % invert checkbox
         pos = get( layoutData.PCDisplay(ii).invertCheckbox, 'Position');
           
         set( layoutData.PCDisplay(ii).invertCheckbox, 'Position',...
            [ windowWidth - pos(3), ...
              windowHeight - ii*layoutData.PCAxesHeight+pos(4)*1.2, ...
              pos(3), pos(4)] );
             
     end
     
     % set position of gamma edit control
     pos = get( layoutData.edGamma, 'Position' );
     set( layoutData.edGamma, 'Position', ...
         [ windowWidth - pos(3), ...
           10+ pos(4), ...
           pos(3), pos(4) ] );
    
    drawnow;
end


function [V,S] = generatePCsFromImageStack( X, fig )
    % reshape to a list of delay scans
    [nrows,ncols,ndelays] = size(X);
    
    delayScans = reshape(X,nrows*ncols, ndelays);
    
    
    % get options
    pcaOptions = getappdata(fig,'pcaOptions');
    
    
    if pcaOptions.subtractChannelMeans
        chanMeans = mean( delayScans );
        delayScans = delayScans - repmat( chanMeans, nrows*ncols, 1 );
        
        % store mean-subtracted version back in X
        
    end
    X = reshape(delayScans,nrows,ncols,ndelays);
    setappdata(fig,'imageStackX',X);
    
    [U,S,V] = svd( delayScans, 0 );
   
    if pcaOptions.scaleByEigenvalues
        V = (S*(V'))';    % multiply by eigenvalues.
        %V = S*V;  % incorrect version that still yields interesting
        %results (see p. 106, 107 of Wilson lab notebook V, 9/20/2011)
    end
    
    
    setappdata(gcf, 'eigenvectors', V);
    setappdata(gcf, 'eigenvalues',diag(S));
    setappdata(gcf,'delayScans',delayScans);
    
    updatePCDisplay( fig );
    
end

function [V,S] = generatePCsFromROIs( X, fig )
    % reshape to a list of delay scans
    % prompt user for an image stack to load
    
    % recall working directory
    wd = getappdata(gcbf,'workingDirectory');
    
    % or use MATLAB's working directory if none has been set yet
    if isempty(wd)
        wd = 'C:/Users/jw295/Data';
    end
    
    [fileName, pathName] = uigetfile([wd,'/*.mat'],...
        'Select a ROI File' );
    
    if ~isequal(fileName,0)
        % remember the working directory for later
        setappdata(gcbf,'workingDirectory',pathName);
        
        % if not canceled, load the file
        load( [pathName,'/',fileName]);
        
    end
   delayScans = data;
   
   [nROIs,ndelays] = size(delayScans);
   
   % get options
    pcaOptions = getappdata(fig,'pcaOptions');
    
    
    if pcaOptions.subtractChannelMeans
        chanMeans = mean( delayScans );
        delayScans = delayScans - repmat( chanMeans, nROIs, 1 );
        
        % store mean-subtracted version back in X
        
    end
   
   
    [U,S,V] = svd( delayScans, 0 );
    
    setappdata(gcf, 'eigenvectors', V);
    setappdata(gcf, 'eigenvalues',diag(S));
    setappdata(gcf,'delayScans',delayScans);
    
    updatePCDisplay( fig );
    
end
 
function PCImageRGB = generatePCImageRGB( X, V, whichPCs, colors, signs,...
                                            enabledPCs, fig )
    [nrows,ncols,~] = size(X);
    
    
    useBipolarColoring = getappdata(fig,'useBipolarColoring');
    
    PCImageRGB = zeros( nrows, ncols, 3 );
    
    % fill in arguments for legacy calls that don't give all args
    if nargin == 4
        signs = 0*whichPCs + 1;
    elseif nargin == 5
        enabledPCs = 0*whichPCs + 1;
    end
        
    for iPC = 1:length( whichPCs )
    %iPC = 3;
        if( enabledPCs(iPC) == 1 )
            % make a matrix of this principal component
            PC_mtrx = permute(repmat( V(:,iPC), [1,nrows,ncols]),[2,3,1]);

            % project the image stack onto this principal component
            PCImage = signs(iPC)*sum(X.*PC_mtrx ,3);

            % eliminate negative values for the moment
            %PCImage( PCImage < 0 ) = 0;

            PCcolor = colors(iPC,:);
            PCImageR = PCcolor(1)*PCImage.*(PCImage > 0);
            PCImageG = PCcolor(2)*PCImage.*(PCImage > 0);
            PCImageB = PCcolor(3)*PCImage.*(PCImage > 0);

            PCImageRGB(:,:,1) = PCImageRGB(:,:,1) + PCImageR;
            PCImageRGB(:,:,2) = PCImageRGB(:,:,2) + PCImageG;
            PCImageRGB(:,:,3) = PCImageRGB(:,:,3) + PCImageB;
            
            if useBipolarColoring
                % mix in the negative projection of this component
                c = [1 1 1] - PCcolor;
                
                PCImageR = c(1)*(-1*PCImage.*(PCImage < 0));
                PCImageG = c(2)*(-1*PCImage.*(PCImage < 0));
                PCImageB = c(3)*(-1*PCImage.*(PCImage < 0));
                
                PCImageRGB(:,:,1) = PCImageRGB(:,:,1) + PCImageR;
                PCImageRGB(:,:,2) = PCImageRGB(:,:,2) + PCImageG;
                PCImageRGB(:,:,3) = PCImageRGB(:,:,3) + PCImageB;
            end
            
            
        end
    end
    

    
    
    setappdata(gcf, 'eigenvectors', V);
end

function setAutoColorScale( fig )
    % set color scale, but do not update image
    % should be called whenever the basis vectors change entirely
    PCImageRGB = getappdata(fig, 'PCImageRGB');
    
    colorScale = max(PCImageRGB(:))*0.9;
    setappdata(fig,'colorScale',colorScale);
    
end

function updatePCImage( fig )
    % re-generate PC image after eigenvectors, color assignments, etc have
    % changed
    
    
    X = getappdata(fig,'ImageStack');
    V = getappdata(fig,'eigenvectors');
     
    colors = getappdata(fig,'PCColors');
    whichPCs = getappdata(fig,'whichPCs');
    signs = getappdata(fig,'signs');
    enabledPCs = getappdata(fig,'enabledPCs');
    
    % generate RGB composite image from singular values
    PCImageRGB = generatePCImageRGB( X, V, whichPCs, colors, signs, ...
        enabledPCs, fig );
    
    setappdata(gcbf, 'PCImageRGB', PCImageRGB);
    
    % update with new histogram
    [histR,histG,histB, bins] = makeHistRGB(PCImageRGB, 512);
    
    histRPatch = findobj( fig, 'Tag', 'histRPatch' );
    set(histRPatch, 'YData', [0,histR,0]);
    
    histGPatch = findobj( fig, 'Tag', 'histGPatch' );
    set(histGPatch, 'YData', [0,histG,0]);
    
    histBPatch = findobj( fig, 'Tag', 'histBPatch' );
    set(histBPatch, 'YData', [0,histB,0]);
    
    
    update( gcf );
end
    

function update( fig )
    hImage = findobj( fig,'Tag', 'componentImage' );
    colorScale = getappdata( fig, 'colorScale' );
    PCImageRGB = getappdata(gcf, 'PCImageRGB');
    
    % scale the image
    PCImageDisplayRGB = PCImageRGB / colorScale;
    
    % apply gamma correction
    gamma = getappdata( fig, 'gamma' );
    PCImageDisplayRGB = PCImageDisplayRGB .^ gamma;
    
    % make sure there are no out-of-bound colors
    PCImageDisplayRGB(PCImageDisplayRGB < 0 ) = 0;
    PCImageDisplayRGB(PCImageDisplayRGB > 1 ) = 1;
    
    set( hImage,'CData', PCImageDisplayRGB );
    
    
    set(findobj(fig,'Tag','histColorScaleLine'),'XData',[0,0]+colorScale);

    
    updatePCDisplay( fig );
    
    drawnow;
end

function updatePCDisplay( fig )
    % redraw displayed principal components on the right-hand side
    layoutData = getappdata( fig, 'layoutData' );
    V = getappdata(fig, 'eigenvectors' );
    colors = getappdata(fig,'PCColors');
    t = getappdata( fig, 't' );             % probe delays
    
    if ~isempty(V)
    
        for ii = 1:layoutData.nPCsDisplayed
            ax = layoutData.PCDisplay( ii ).PCAxes;
            h = findobj(ax,'Tag','l_eig');

            eigenvector = V(:,ii);

            if isempty( h )
                % if these lines haven't yet been displayed, make new ones
                axes(ax);
                hline = line(t, eigenvector,'color',colors(ii,:),'tag','l_eig');
                xlim( [min(t), max(t)] );
                
                % make a line to mark y=0
                line([min(t),max(t)],[0,0],'color','k','linewidth',0.5);
            else
                % update existing line with new data
                set(h(1),'XData', t, 'YData', eigenvector,'color',colors(ii,:));
            end
        end
    end
end


%==========================================================================
% FIGURE CALLBACKS
function resizeCallback(src, evt)
    layout(gcbf);
end

%==========================================================================
% MOUSE DOWN INTERACTION CALLBACKS
%
function buttonDown_histColorScaleLine(src, evt)
% src - the object that is the source of the event
% evnt - empty for this property
    sel_typ = get(gcbf,'SelectionType');
    
    switch sel_typ
        case 'normal'
            % single click with left button-- track mouse motions
            set(gcf,'WindowButtonMotionFcn',{@(s,e)buttonMotion_histColorScaleLine(s,e,src)});
            set(gcf,'WindowButtonUpFcn',@button_up);
    end
end

function buttonMotion_histColorScaleLine(src,evnt,histColorScaleLine)
    % find out where mouse is
    cp = get(gca,'currentPoint');
    x = cp(1,1);
    y = cp(1,2);
    
    setappdata( gcbf, 'colorScale',  x);

    update(gcbf);
end

function button_up(src,evnt)
    set(gcf,'WindowButtonMotionFcn','');
end

%==========================================================================
% MENU CALLBACKS
%
function mnuComputePCFromSVD(s, e)
    set(gcbf,'Pointer','watch');
    drawnow;

    X = getappdata(gcf,'ImageStack');
    [V,S] = generatePCsFromImageStack( X, gcbf );
    setappdata(gcbf,'eigenvectors',V);
     
    colors = [1,0,0; 0,1,0; 0,0,1]; % RGB
    setappdata(gcbf,'PCColors',colors);
    whichPCs = [1,2,3]; % which principal components to use
    setappdata(gcbf,'whichPCs', whichPCs);
    signs = [1,1,1];
    setappdata(gcbf,'signs',signs);
    
    updatePCImage(gcbf);    % re-generate the principal component image
    update(gcbf);
    
    renormalize(gcbf);

    set(gcbf,'Pointer','arrow');
end

function mnuLoadPCs( s, e )
    % recall working directory
    wd = getappdata(gcbf,'workingDirectory');
    
    if ispref('puprisa','PCWorkingDirectory')
        wd = getpref('puprisa','PCWorkingDirectory');
    else
        % or use MATLAB's working directory if none has been set yet
        wd = pwd();
    end
    
    
    % or use MATLAB's working directory if none has been set yet
    if isempty(wd)
        wd = 'C:\Users\jw295\Data';
    end
    
    [fileName, pathName] = uigetfile([wd,'\*.mat'],...
        'Select a Principal Components File' );
    
    if ~isequal(fileName,0)
        % remember the working directory for later
        setpref('puprisa','PCWorkingDirectory',pathName);
  
        % if not canceled, load the file
        dat = load( [pathName, '\',fileName], 'coeff' );
        
    end
    
    if isempty(dat)
        error('invalid file-- must contain the variable coeff');
    else
        coeff = dat.coeff;
    end
    V = coeff;
    
    set(gcbf,'Pointer','watch');
    drawnow;

   setappdata(gcbf,'eigenvectors',V);
     
    colors = [1,0,0; 0,1,0; 0,0,1]; % RGB
    setappdata(gcbf,'PCColors',colors);
    whichPCs = [1,2,3]; % which principal components to use
    setappdata(gcbf,'whichPCs', whichPCs);
    signs = [1,1,1];
    setappdata(gcbf,'signs',signs);
    
    updatePCImage(gcbf);    % re-generate the principal component image
    
    set(gcbf,'Pointer','arrow');
    update( gcbf );
    
  
end

function mnuPCsFromROI( s, e )
    set(gcbf,'Pointer','watch');
    drawnow;

    X = getappdata(gcf,'ImageStack');
    [V,S] = generatePCsFromROIs( X, gcbf );
    setappdata(gcbf,'eigenvectors',V);
     
    colors = [1,0,0; 0,1,0; 0,0,1]; % RGB
    setappdata(gcbf,'PCColors',colors);
    whichPCs = [1,2,3]; % which principal components to use
    setappdata(gcbf,'whichPCs', whichPCs);
    signs = [1,1,1];
    setappdata(gcbf,'signs',signs);
    
    updatePCImage(gcbf);
    update(gcbf);
    
    renormalize(gcbf);

    
    %updatePCImage(gcbf);    % re-generate the principal component image
    set(gcbf,'Pointer','arrow');
    %update( gcbf );
end

function mnuSwitchPCColor1and2 ( src, evt )

pcc = getappdata(gcf,'PCColors')

pccNew(1,:) = pcc(2,:);
pccNew(2,:) = pcc(1,:);
pccNew(3,:) = pcc(3,:);

setappdata(gcf,'PCColors',pccNew)

updatePCImage(gcbf);
    update(gcbf);
    
end
function mnuSwitchPCColor2and3 ( src, evt )

pcc = getappdata(gcf,'PCColors')

pccNew(1,:) = pcc(1,:);
pccNew(2,:) = pcc(3,:);
pccNew(3,:) = pcc(2,:);

setappdata(gcf,'PCColors',pccNew)

updatePCImage(gcbf);
    update(gcbf);
    
end
function mnuOptScaleByEigval( src, evt )
    % get options
    pcaOptions = getappdata(gcbf,'pcaOptions');
    
    % toggle setting to scale by eigenvalues
    pcaOptions.scaleByEigenvalues = 1 - pcaOptions.scaleByEigenvalues;
    setappdata(gcbf,'pcaOptions',pcaOptions);
    
    % change checkbox on this menu item to reflect
    if pcaOptions.scaleByEigenvalues
        set(src,'Checked','on');
    else
        set(src,'Checked','off');
    end
end

function mnuOptSubtractChanMeans(src, evt )
    % get options
    pcaOptions = getappdata(gcbf,'pcaOptions');
    
    % toggle setting to subtract channel means
    pcaOptions.subtractChannelMeans = 1 - pcaOptions.subtractChannelMeans;
    setappdata(gcbf,'pcaOptions',pcaOptions);
    
    % change checkbox on this menu item to reflect
    if pcaOptions.subtractChannelMeans
        set(src,'Checked','on');
    else
        set(src,'Checked','off');
    end
end

%==========================================================================
% CALLBACKS FOR CONTROLS/SETTINGS FOR PCS
function enabledCheckboxCallback( s, e )
   % enable/disable a particular principal component channel
   iPCChannel = get(s,'UserData');  % which channel is stored in userdata
   enableChannel = get(s,'Value');        % 0 unchecked, 1 checked
   
   enabledPCs = getappdata(gcbf,'enabledPCs');
   
   enabledPCs( iPCChannel ) = enableChannel;
   
   setappdata( gcbf,'enabledPCs',enabledPCs );
   
   updatePCImage( gcbf );
end

function invertCheckboxCallback( s, e )
   % toggle invert of a particular principal component channel
   iPCChannel = get(s,'UserData');  % which channel is stored in userdata
   invertChannel = get(s,'Value');        % 0 unchecked, 1 checked
   
   signs = getappdata(gcbf,'signs');
   
    if( invertChannel == 0 )
        signs( iPCChannel ) = 1;
    else
        signs( iPCChannel ) = -1;
    end
   
   setappdata( gcbf,'signs',signs);
   
   updatePCImage( gcbf );
end

function mnuChangeScaleBar( src, evt )
    % allows user to change the size of the scale bar to a specified size
    % prompt the size of the scale bar 
    prompt = {'Scale Bar Size:'};
    dlg_title = 'Scale Bar Size';
    num_lines = 1;
    def = {'100'};
    answer = inputdlg(prompt,dlg_title,num_lines,def);
    sb = str2num(answer{1});
    if isempty(sb)
        error('invalid entry');
    end
    
    % retrieve scale bar information
    
    header = getappdata( gcbf,'header');
    if ispref('puprisa','scaleX_volts_per_micron')
        scaleX_volts_per_micron = ...
            getpref('puprisa','scaleX_volts_per_micron');
%         scaleY_volts_per_micron = ...
%             getpref('puprisa','scaleY_volts_per_micron');
        
        nx = header.pixelsperline;
        scanRangeX = header.scanrangex;
        
        scaleX_px_per_micron = scaleX_volts_per_micron * nx / scanRangeX;
        
        setappdata(gcbf,'scaleX_px_per_micron',scaleX_px_per_micron);
        % delete old scale bar
         h = findobj(gcf,'tag','scaleBarGroup');
         delete(h(1));
         % switch axes to imageAxes so the scalebar is drawn on the image only
    
         h = findobj(gcf,'tag','imageAxes');
         axes(h(1));
        % draw scale bar of the specified size
        puprisa_scaleBar(sb,scaleX_px_per_micron);
    end
end


function mnuEumelaninPercentage ( src, evt )

    % make a mask to exclude ink from calculation
    h2 = findobj(gcf,'tag','componentImage');
    axes(get(h2,'Parent'));     % set current axes to the displayed image
    
    % get a region of interest to exclude the ink
    hROI = getappdata(gcbf,'ROI');  % get stored handle to region of interest (ROI)
    if isempty( hROI )              % if the ROI has not yet been created...
        hROI = impoly;              % then make a new ROI
        setappdata(gcbf,'ROI',hROI);% and remember its handle for next time
    end
    
    m = createMask(hROI,h2);

    % calculate eumelanin and pheomelanin quantities
    c = getappdata(gcbf,'PCImageRGB');
    PC1image = c(:,:,1);
    PC2image = c(:,:,2);
    
    %method = 'ours';
    %method = 'tomOld';
    method = 'Paco';
    
    switch( method )
    
        case 'ours'
            % our method
            euImage = (0.1243.*(PC1image))+(0.0345.*PC2image);
            pheoImage = (-0.0954.*(PC1image))+(0.1775.*(PC2image));
            
        case 'tomOld'
            % from Tom's code (old data set)
            euImage = (PC1image + (0.537 .* PC2image)) ./ 0.143;
            pheoImage = (PC1image - (3.603 .* PC2image)) ./ -0.735;
            % correction from re-projection of principal components onto eu/phoe
            % cuvette spectra (See p. 295-299 of Wilson lab Notebook V, 1/11/2012)
              pheoImage = pheoImage * 0.735 / 0.2397;
    
        case 'tomNew'
            % from Tom's code (new data set)
            euImage = (PC2image - (9.27 .* PC1image)) ./ -37.88;
            pheoImage = (PC2image - (-0.175 .* PC1image)) ./ 8.058;
            % correction from re-projection of principal components onto eu/phoe
            % cuvette spectra (See p. 295-299 of Wilson lab Notebook V, 1/11/2012)
              pheoImage = pheoImage * 0.735 / 0.2397;
        case 'Paco'   
            
            euImage = 20.2750.*PC1image +  7.9010.* PC2image ;
            pheoImage = -5.0367*PC1image + 13.6216 .* PC2image ;
            
    end
    
    % consider pixels with negative values to be zero
    validMask = ((euImage >= 0) & (pheoImage >= 0));
    euImage( validMask == 0 ) = 0;
    pheoImage( validMask == 0 ) = 0;
    
    % eliminate negative eu and pheo components
    %euImage( euImage < 0 ) = 0;
    %pheoImage( pheoImage < 0 ) = 0;
    

    % make an image with eumelanin percent as the color scale
    % fracImage = (euImage)./((pheoImage)+(euImage));
    % figure;
    % imagesc(fracImage);

    % subtract masked area
    euImageMasked = euImage.*(1-m);
    pheoImageMasked = pheoImage.*(1-m);
    fracImageMasked = (euImageMasked)./((pheoImageMasked)+(euImageMasked));
    totalImageMasked = pheoImageMasked+euImageMasked;
    
    % figure('name','eumelanin');
    % imagesc(euImageMasked);
    % set(gca,'ydir','normal');
    % colormap(cmap(256,[0 0 0],[1 0 0]));
    
    % figure('name','pheomelanin');
    % imagesc(pheoImageMasked);
    % set(gca,'ydir','normal');
    % colormap(cmap(256,[0 0 0], [0 1 0]));
    
     % A = size(pheoImageMasked)
     % B = A(1)
    %  M = zeros(B,B,3);
     % M(:,:,1)=euImageMasked;
     % M(:,:,2)=pheoImageMasked;
     % figure('name','Eumelanin and Pheomelanin Concentration Image');
     % imshow(M)
    
    figure('Name','Fractional Eumelanin','inverthardcopy','off');
    im = imagesc(fracImageMasked*100);
    set(gca,'ydir','normal','xtick',[],'ytick',[]);
    totalImage = ((pheoImageMasked)+(euImageMasked));
    set(gca,'color','k');
    %colormap(cmap(256,[1 0 1], [0 1 1], [1 1 0]));
    %caxis([0,100]);
    colormap(cmap(256,[0 1 0], [1 1 0], [1 0 0]));
    caxis([0,100]);
    axis image;
    h=colorbar;
    ylabel(h, 'Eumelanin percentage');
    
    % set up transparency to indicate melanin concentration
    set(im,'alphadata',totalImage);
    set(im,'alphadatamapping','scaled');
    alim([min(totalImageMasked(:)),max(totalImageMasked(:))]);
    % set up user interface controls for adjusting transparency
    makeTransparencyControls(); 
    
    % calculate bulk eumelanin fraction
    melFrac = sum(euImageMasked(:))/(sum(euImageMasked(:))+sum(pheoImageMasked(:)));
    disp(['Bulk eumelanin fraction: ', num2str(100*melFrac), '%']);
    
    % calculate histogram
    % figure('name','histogram');
    % hist(fracImageMasked(:), linspace(0,1));
    % title('Histogram, fractional eu')
    % ylabel('N pixels');
    % xlabel('Eumelanin Fraction');
    
    % calculate weighted histogram
    figure('name','histogram (weighted)');
    
    [n,bin] = histc(fracImageMasked(:), linspace(0,1));
    
    % now accumulate weights
    w = 0*n;
    for ii = 1:length(n)-1;
        w(ii) = sum(totalImageMasked( bin == ii) );
    end
    
    h = bar(linspace(0,1),w);
    set(h,'tag','histogramBarGraph');
    xlabel('Eumelanin fraction');
    ylabel('Sum pixel intensity');
    title('Intensity-weighted Histogram')
    axis tight
    %hpb = uicontrol('style','pushbutton',...
    %    'string','Delete 0% and 100% bins','position',[10,10,150,20],...
    %    'callback',@delete0And100HistogramBins);
    
    
    % calculate stats

    xd = linspace(0,1).';
    yd = w;
    
    %xd = xd(2:(end-1));
    %yd = yd(2:(end-1));
    
    %set(h,'xdata',xd);
    %set(h,'ydata',yd);
    
    f = yd / sum(yd);   % normalized probability distribution
    
    %delete(src);
    
    %title('Intensity-weighted Histogram (0% and 100% bins deleted)');
    %drawnow;
    
    theMean = sum(xd.*f);
    
    %theStdDev = sqrt( (1/length(xd)) * sum( (xd - theMean).^2) );
    theStdDev = sqrt( sum( (xd-theMean).^2.*f ) );
    
    pr = yd / sum(yd);  % rearrange as probability
    
    pr = pr( pr > 0 );  % ignore zero values
    entr = -sum( pr.*log2(pr) );
    
    disp(['Mean eumelanin fraction: ', num2str(theMean)]);
    disp(['Std. dev.: ', num2str(theStdDev)]);
    disp(['Information entropy: ', num2str(entr)]);
    
    % calculate 2D scatterplot of percent eu vs. intensity
    figure('Name','scatterplot');
    line( fracImageMasked(:)*100, totalImageMasked(:), 'linestyle','none','marker','.',...
        'MarkerSize',2);
    xlabel('Eumelanin percentage');
    ylabel('Pixel brightness');
    
    % ask for threshold level
    prompt = {'Threshold level:'};
    dlg_title = 'Threshold level';
    num_lines = 1;
    def = {'0.00'};
    answer = inputdlg(prompt,dlg_title,num_lines,def);
    tl = str2num(answer{1});
    if isempty(tl)
        error('invalid entry');
    end
    
    % subtract threshold level
    validMask = (totalImageMasked >= tl);
    totalImageMasked( validMask == 0 ) = 0;
    figure('Name','scatterplot');
    line( fracImageMasked(:)*100, totalImageMasked(:), 'linestyle','none','marker','.',...
        'MarkerSize',2);
    xlabel('Eumelanin percentage');
    ylabel('Pixel brightness');

     % calculate weighted histogram with threshold subtracted
    figure('name','histogram (weighted with threshold)');
    
    [n,bin] = histc(fracImageMasked(:), linspace(0,1));
    
    % now accumulate weights
    w = 0*n;
    for ii = 1:length(n)-1;
        w(ii) = sum(totalImageMasked( bin == ii) );
    end
    
    h = bar(linspace(0,1),w);
    set(h,'tag','histogramBarGraph');
    xlabel('Eumelanin fraction');
    ylabel('Sum pixel intensity');
    title('Intensity-weighted Histogram')
    axis tight
    
end

function displayStatistics( src, evt )
    % delete bins on extreme ends of percent eumelanin histogram
    h = findobj(gcbf,'tag','histogramBarGraph');
    
    xd = get(h,'xdata');
    yd = get(h,'ydata');
    
    % xd = xd(2:(end-1));
    % yd = yd(2:(end-1));
    
    % set(h,'xdata',xd);
    % set(h,'ydata',yd);
    
    f = yd / sum(yd);   % normalized probability distribution
    
    % delete(src);
    
    % title('Intensity-weighted Histogram (0% and 100% bins deleted)');
    % drawnow;
    
    theMean = sum(xd.*f);
    
    %theStdDev = sqrt( (1/length(xd)) * sum( (xd - theMean).^2) );
    theStdDev = sqrt( sum( (xd-theMean).^2.*f ) );
    
    pr = yd / sum(yd);  % rearrange as probability
    
    entr = -sum( pr.*log2(pr) );
    
    disp(['Mean eumelanin fraction: ', num2str(theMean)]);
    disp(['Std. dev.: ', num2str(theStdDev)]);
    disp(['Information entropy: ', num2str(entr)]);
    

end

function makeTransparencyControls()
    % set up controls to adjust transparency range
    
    pad = 5;    % distance between controls
    
    % label
    h1 = uicontrol('Style','text','String','alim:');
    pos1 = get(h1,'position');
    
    % minimum alpha limit edit box
    h2 = uicontrol('Style','edit','String',num2str(min(alim())));
    pos2 = get(h2,'position');
    pos2(1) = pos1(1) + pos1(3) + pad;
    set(h2,'position',pos2);
    
    % maximum alpha limit edit box
    h3 = uicontrol('Style','edit','String',num2str(max(alim())));
    pos3 = get(h3,'position');
    pos3(1) = pos2(1) + pos2(3) + pad;
    set(h3,'position',pos3);

    % set up callbacks
    set(h2, 'tag', 'edit_alphaMin', 'callback', @alimCallback);
    setappdata(h2,'ax',gca());
    set(h3, 'tag', 'edit_alphaMax', 'callback', @alimCallback);
    setappdata(h3,'ax',gca());
    
end

function alimCallback( src, evt )
    % callback function for adjusting alpha limits of an image
    
    % find alpha min
    hMin = findobj( gcbf, 'tag', 'edit_alphaMin' );
    alphaMin = str2num( get(hMin, 'string') );
    
    ax = getappdata(hMin,'ax');

    hMax = findobj( gcbf, 'tag', 'edit_alphaMax' );
    alphaMax = str2num( get(hMax,'string') );
    
    % verify that the text entered are valid numbers
    isValid = 1;
    if( isempty( alphaMin ) || isempty( alphaMax ) )
        isValid = 0;
    elseif( alphaMax <= alphaMin )
        isValid = 0;
    end
    
    if( isValid == 1)
       % set alim accordingly
       set(ax,'alim',[alphaMin, alphaMax]);
    else
        % re-set textboxes to valid numbers
        al = get(ax,'alim');
        alphaMin = al(1);
        alphaMax = al(2);
        
        set(hMin,'String',num2str(alphaMin));
        set(hMax,'String',num2str(alphaMax));
    end
   
    drawnow;    % refresh the figure
end

function gammaChanged( src, evt )
    % callback for adjusting image gamma
    strGamma = get(src, 'String' );
    numGamma = str2num( strGamma );
    
    if isempty( numGamma )
        set(src, 'String', num2str( getappdata(gcbf,'gamma') ) );
    else
        setappdata(gcbf,'gamma',numGamma);
        updatePCImage( gcbf );
    end
    
end

function mnuFourierTransform ( src, evt )

% image the fourier transform of the principal component image

    c = getappdata(gcf,'PCImageRGB');
    [nr,nc,~] = size(c);
    
    r = c(:,:,1);
    g = c(:,:,2);
    b = c(:,:,3);
    
    % get rid of negative components
    r(r < 0) = 0;
    g(g < 0) = 0;
    b(b < 0) = 0;
    
    % multiply by an fft window to prevent minimize aliasing
    wy = repmat(blackman(nr),1,nc);
    wx = repmat(blackman(nc).',nr,1);
    w = wy.*wx;
    
    
    r_f = abs(fftshift(fft2(r.*w)));
    g_f = abs(fftshift(fft2(g.*w)));
    b_f = abs(fftshift(fft2(b.*w)));
    
    imFFT = zeros(nr,nc,3);
    imFFT(:,:,1) = log(r_f);
    imFFT(:,:,2) = log(g_f);
    imFFT(:,:,3) = log(b_f);
    
    imFFT = imFFT / max(imFFT(:));
    
    figure('name','fourier transform');
    
    imshow(imFFT);
    %colormap(cmap(256,[0 0 0], [1 1 1]));

end

function mnuBipolarComp( src, evt )
    % set preferences to use a bipolar coloring scheme for a more 
    % complete preview
    chk = get(src,'Checked');
    
    switch chk
        case 'on'
            set(src,'Checked','off');
            setappdata(gcbf,'useBipolarColoring',false);
        case 'off'
            set(src,'Checked','on');
            setappdata(gcbf,'useBipolarColoring',true);
    end
    
    updatePCImage( gcbf );
end
