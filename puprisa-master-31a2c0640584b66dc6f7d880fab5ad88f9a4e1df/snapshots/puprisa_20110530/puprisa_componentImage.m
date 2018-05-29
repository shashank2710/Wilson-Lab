% NEXT: load in correct time axes!
% MAYBE: extend to more than 3PCS, give option for primary color segmenting
%       (classify pixel by strongest fit)
% maybe give option of scaling the PCs individually
% hint: use gcbf to interact with calling figure

function puprisa_componentImage( X )
    % generate a princpal component image
    % overlay components with selected colors for each component
    
    initFigure(X);
    initMenus();
    
end   

%==========================================================================
% INITIALIZATION ROUTINES
%
 function initFigure( X )
    fig = figure('Name','Component Image','Units','pixels',...
        'ResizeFcn',@resizeCallback);
    
    layoutData = initLayout( fig );
    
    % lay out axes locations
    layout(fig);
    
    setappdata(fig,'ImageStack',X);
      
    [V,S] = generatePCsFromImageStack( X );
    setappdata(fig,'eigenvectors',V);
     
    colors = [1,0,0; 0,1,0; 0,0,1]; % RGB
    setappdata(fig,'PCColors',colors);
    whichPCs = [1,2,3]; % which principal components to use
    setappdata(fig,'whichPCs', whichPCs);
    signs = [1,1,1];
    setappdata(fig,'signs',signs);
    enabledPCs = [1,1,1];   % which PC channels are enabled
    setappdata(fig,'enabledPCs',enabledPCs);
    
    % generate RGB composite image from singular values
    PCImageRGB = generatePCImageRGB( X, V, whichPCs, colors, signs );
    
    setappdata(fig, 'PCImageRGB', PCImageRGB);
    
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
 function initMenus()
    h = uimenu(gcf,'Label','PrincipalComponents');
    uimenu(h, 'Label', 'Compute from image SVD', 'Callback',@mnuComputePCFromSVD );
    uimenu(h, 'Label', 'Load Principal Components...', 'Callback',@mnuLoadPCs );
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
            [ windowWidth - layoutData.rightColumnWidth, ...
              windowHeight - ii*layoutData.PCAxesHeight, ...
              pos(3), pos(4)] );
          
         % invert checkbox
         pos = get( layoutData.PCDisplay(ii).invertCheckbox, 'Position');
           
         set( layoutData.PCDisplay(ii).invertCheckbox, 'Position',...
            [ windowWidth - layoutData.rightColumnWidth, ...
              windowHeight - ii*layoutData.PCAxesHeight+pos(4)*1.2, ...
              pos(3), pos(4)] );
             
     end
    
    drawnow;
end


function [V,S] = generatePCsFromImageStack( X )
    % reshape to a list of delay scans
    [nrows,ncols,ndelays] = size(X);
    delayScans = zeros(nrows*ncols, ndelays);
    for irow = 1:nrows
        minIndex = (irow-1)*ncols+1;
        maxIndex = (irow)*ncols;
        delayScans( minIndex:maxIndex, : ) = squeeze(X( irow, 1:ncols, : ));
    end
    
    [U,S,V] = svd( delayScans, 0 );
    
    setappdata(gcf, 'eigenvectors', V);
    setappdata(gcf, 'eigenvalues',diag(S));
    
end
 
function PCImageRGB = generatePCImageRGB( X, V, whichPCs, colors, signs,...
                                            enabledPCs )
    [nrows,ncols,~] = size(X);
    
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
            PCImage( PCImage < 0 ) = 0;

            PCcolor = colors(iPC,:);
            PCImageR = PCcolor(1)*PCImage;
            PCImageG = PCcolor(2)*PCImage;
            PCImageB = PCcolor(3)*PCImage;

            PCImageRGB(:,:,1) = PCImageRGB(:,:,1) + PCImageR;
            PCImageRGB(:,:,2) = PCImageRGB(:,:,2) + PCImageG;
            PCImageRGB(:,:,3) = PCImageRGB(:,:,3) + PCImageB;
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
        enabledPCs );
    
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
    
    % make sure there are no out-of-bound colors
    PCImageDisplayRGB(PCImageDisplayRGB < 0 ) = 0;
    PCImageDisplayRGB(PCImageDisplayRGB > 1 ) = 1;
    
    set( hImage,'CData', PCImageDisplayRGB );
    
    
    set(findobj(fig,'Tag','histColorScaleLine'),'XData',[0,0]+colorScale);

    
    updatePCDisplay( fig );
    
    drawnow;
end

function updatePCDisplay( fig )
    % redraw displayed principal components
    layoutData = getappdata( fig, 'layoutData' );
    V = getappdata(fig, 'eigenvectors' );
    colors = getappdata(fig,'PCColors');
    
    if ~isempty(V)
    
        for ii = 1:layoutData.nPCsDisplayed
            ax = layoutData.PCDisplay( ii ).PCAxes;
            h = get(ax,'Children');

            eigenvector = V(:,ii);
            t = 1:length( eigenvector );

            if isempty( h )
                axes(ax);
                hline = line(t, eigenvector,'color',colors(ii,:));
            else
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
    [V,S] = generatePCsFromImageStack( X );
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

function mnuLoadPCs( s, e )
    % recall working directory
    wd = getappdata(gcbf,'workingDirectory');
    
    % or use MATLAB's working directory if none has been set yet
    if isempty(wd)
        wd = 'C:\Users\jw295\Data';
    end
    
    [fileName, pathName] = uigetfile([wd,'\*.mat'],...
        'Select a Principal Components File' );
    
    if ~isequal(fileName,0)
        % remember the working directory for later
        setappdata(gcbf,'workingDirectory',pathName);
        
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