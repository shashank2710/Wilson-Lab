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

function puprisa( fileName )
%    close all;
    
    % set up the figure with 4 panels and a statusbar
    f = figure('Name','puprisa','Units','pixels');
    initFigure( f );
    
    initMenus( f );
    
    updateStatus('puprisa ready.');
    
    if nargin == 1
        loadData(f, fileName);
    end
    
end

%--------------------------------------------------------------------------
% INITIALIZATION ROUTINES
%
function initFigure( f )
    set(f,'Units','pixels');
    fpos = get(f,'Position');
    headerTextWidth = 200;
    fpos(3) = fpos(3)+headerTextWidth+20;
    fwidth = fpos(3);
    set(f,'Position',fpos);
    
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
        setappdata(hChannelAxes(iChannel),'Channel',iChannel);
        hChannelImages(iChannel) = image(0,'CDataMapping','scaled',...
                                    'ButtonDownFcn',@axesButtonDownFcn);
        set(hChannelAxes(iChannel), 'PlotBoxAspectRatio',[1 1 1],...
                                    'DataAspectRatioMode','auto',...
                                    'CLimMode','auto',...
                                    'YDir','normal',...
                                    'ButtonDownFcn',@axesButtonDownFcn,...
                                    'XTick',[],'YTick',[]);
                                
        ylabel(['Channel ',num2str(iChannel)]);
        
        axesContextMenu( hChannelImages(iChannel) );
    end
    
    setappdata(gcf,'hChannelAxes',hChannelAxes);
    setappdata(gcf,'hChannelImages',hChannelImages);
    
    % set up status bar
    hStatusText = uicontrol('style','text','String','ImageAnalysisJ',...
        'Units','pixels','tag','hStatusText','position',[0,0,fwidth,20]);
    
    % set up header text
    hHeaderText = uicontrol('style','edit','String','No File Opened',...
        'Units','pixels','tag','hHeaderText',...
        'position',[0,0,headerTextWidth,20],...
        'Max',2,'HorizontalAlignment','left');

    % resize everything
    setChannelAxesPos(gcf);
    % set up callbacks
    set(gcf,'ResizeFcn',@windowResizeFcn)
       
    % set up scroll wheel fn to be able to flip between slices
    set(gcf,'WindowScrollWheelFcn', @scroll_wheel);
    
    % set up keypress handler for flipping between slices
    set(gcf,'KeyPressFcn', @keyPress);
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
    uimenu(f,'Label','Calibrate scale...','Callback',@mnuCalibrateScale);
    uimenu(f,'Label','About puprisa','Callback',@mnuAbout); 
    
    f = uimenu('Label','TestImages');
    uimenu(f,'Label','Biexponential Test Set','callback',@mnuBiexpTestSet);
       
    f = uimenu('Label','Art Imaging');
    uimenu(f,'Label','Multiple files delayConst Histogram','callback',@mnuLoadMultipleFiles);
    uimenu(f,'Label','Multiple files Average delayConst Histogram','callback',@mnuLoadMultipleFilesAVG);
    uimenu(f,'Label','Average all pixels over multiple Delaystacks in a directory','callback',@mnuLoadMultipleFilesAVG2);
    uimenu(f,'Label','Average all pixels in an image for multiple files in a directory','callback',@mnuLoadMultipleFilesAVG3);
    uimenu(f,'Label','Same area average in multiple files in a directory','callback',@mnuLoadMultipleFilesSameArea);
    uimenu(f,'Label','Multiple files average, Ratio','callback',@mnuLoadMultipleFilesRatio);
    uimenu(f,'Label','Average all pixels over multiple Delaystacks in a directory (NEW)','callback',@mnuLoadMultipleFilesAVG2new);
   
    f = uimenu('Label','Scale calibration');
    uimenu(f,'Label','Scale calibration manager...','Callback',@mnuCalibrateScale);

    % populate this menu w/ available scales
    % don't forget checkmarks on current calibration
    if ispref('puprisa','scaleCalibrationSet')
        s = getpref('puprisa','scaleCalibrationSet');
        for is = 1 : length(s)
            if is == 1
                sep = 'on';
            else
                sep = 'off';
            end
            
            uimenu(f,'Label',s(is).name,'Separator',sep,'Callback',@(s,e) mnuChangeScale(s,e,is));
        end
    end
end

function mnuChangeScale(src,evt,is)
    s = getpref('puprisa','scaleCalibrationSet');

    
    setpref('puprisa','scaleX_volts_per_micron',s(is).scaleX_V_per_micron);
    setpref('puprisa','scaleY_volts_per_micron',s(is).scaleX_V_per_micron);
end

function axesContextMenu( ax )
    % make a context menu to appear on right-click of an axes object
    hcmenu = uicontextmenu();
    
    hcb1 = uimenu(hcmenu, 'Label', 'Open as delay stack...', ...
        'Callback',@cmnuOpenAsDelayStack);
    hcb1 = uimenu(hcmenu, 'Label', 'Open as z stack, en face...', ...
        'Callback',@cmnuOpenAsZStackEnFace);
    hcb1 = uimenu(hcmenu, 'Label', 'Open as z stack, cross section...', ...
        'Callback',@cmnuOpenAsZStackCrossSection);

    set(ax, 'uicontextmenu', hcmenu );
end

function cmnuOpenAsDelayStack( src, evt )
    % forces stack to be opened as delay stack, even if it is a z stack
    ax = get(gco,'parent');
    fileName = getappdata(gcbf,'fileName');

    imageStack = getappdata(ax,'imageStack');

    header = getappdata(gcbf,'fileHeader');

    chan = getappdata(ax,'Channel');
    zPos = getappdata(gcbf,'zPos');
    puprisa_viewChannel( imageStack, zPos, fileName, header, chan );
end

%--------------------------------------------------------------------------
% MISC CALLBACK FUNCTIONS
%
% Key presses, window resize, mouse button down, etc.
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

function scroll_wheel( src, evnt )
    incrementSlice(evnt.VerticalScrollCount);
end

function incrementSlice( incr )
% Change the currently displayed image slice by the given increment.
%
% This function is called whenever the user moves the scroll wheel or 
% presses the left or right arrow key.

    ik = getappdata(gcbf,'currentSlice') + incr;
    
    % Ensure the new index is within bounds.
    nSlices = getappdata(gcbf,'nSlices');
    if ik < 1
        ik = 1;
    elseif ik > nSlices
        ik = nSlices;
    end

    setappdata(gcbf,'currentSlice', ik);
    updateAll( gcbf );
end

function keyPress( src, evnt )
% handle press of arrow keys to cycle through slices in the stack
    incr = 0;
    
    if( evnt.Character == char(28) )
        % left arrow
        incr = -1;
    elseif( evnt.Character == char(29) )
        % right arrow
        incr = +1;
    end
    
    if incr ~= 0
        incrementSlice( incr );
    end
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
            
            header = getappdata(gcbf,'fileHeader');
            
            chan = getappdata(ax,'Channel');
            
            % use the appropriate viewer depending on the stack type
            switch stackType
                case 'delay stack'
                    zpos = getappdata(gcbf,'delays');
                    puprisa_viewChannel( imageStack, zpos, fileName,...
                                         header, chan );
                case 'z stack'
                    zpos = getappdata(gcbf,'zPos');
                    puprisa_viewChannelZStack(imageStack, zpos, fileName, header);
                case 'mosaic'
                    xPos = getappdata(gcbf,'xPos');
                    yPos = getappdata(gcbf,'yPos');
                    zPos = getappdata(gcbf,'delays');

                    puprisa_viewChannelMosaic(imageStack, xPos, yPos, zPos,...
                        fileName, header, chan);
            end
    end
end

%--------------------------------------------------------------------------
% MENU CALLBACK FUNCTIONS
%
function mnuAbout(src, evt)
    msgbox(['PUPRISA: PUmp PRobe Image Stack Analysis. ',...
            'Warren Lab: Duke University. ',...
            'Created 2011 by J. W. Wilson. ',...
            'Contributions by, P. Samineni, M.J. Simpson)'],...
        'About PUPRISA', 'help');
end

function mnuAutoScaleStack(src, evt)
% set clim on display to consider entire stack
% (NOT IMPLEMENTED)
    get(get(get(src,'Parent'),'Parent'))
end

function buttonDown_colorbar( src,evt,ax,  appDataTag_imageStack )
% set colormap
% (NOT IMPLEMENTED)
    
end

function mnuOpenImageStack( src, evt )
% Prompt user for an image stack and load it
    
    % recall working directory
    if ispref('puprisa','imageWorkingDirectory')
        wd = getpref('puprisa','imageWorkingDirectory');
    else
        % or use MATLAB's working directory if none has been set yet
        wd = pwd();
    end
    
    % prompt for filename
    [fileName, pathName] = uigetfile('*.dat;*.lsm;*.tif',...
        'Select an Image Stack',...
        [wd,'\']);
    
    if ~isequal(fileName,0)
        % remember the working directory for later
        setpref('puprisa','imageWorkingDirectory',pathName);
        
        % if not canceled, load the file
        loadData(gcbf, [pathName,fileName]);
    end
end

function mnuFrameShiftCorrection( src, evt )
    % prompt for alignment channel
    alignChanStr = inputdlg(...
        'Channel number for calculating alignment:','Frame shift registration.', 1, {'1'});
    alignChan = str2num(alignChanStr{1});
    
    if isempty( alignChan )
        error('invalid channel number');
    end
    
    % align slices based on cross-correlation of y-channel
    hChannelAxes = getappdata(gcbf,'hChannelAxes');

    Y = getappdata(hChannelAxes(alignChan),'imageStack');
    
    [nrows,ncols,nslices] = size(Y);
    
    nChan = getappdata(gcbf,'nChannels');
    
    % get all channel data
    for iChan = 1:nChan
        %chanData{iChan} = getappdata(hChannelAxes(1),'imageStack');
        chanData{iChan} = getappdata(hChannelAxes(iChan),'imageStack');
    end
    
    hwb = waitbar(0,'Aligning slices...');
    
    for islice = 2:nslices
        
        Y2 = Y(:,:,islice);
        %Y2 = circshift(Y1,[20,13]);
        
        

        [output, shifted] = dftregistration(...
            fft2(Y(:,:,1)), fft2(Y2),100);
        
        Y2 = real(ifft2(shifted));
        Y(:,:,islice) = Y2;
        
        % shift X accordingly
        %X2 = circshift(X2,[shifty,shiftx]);
        %
        diffphase = output(2);
        row_shift = output(3);
        col_shift = output(4);
        
        for iChan = 1:nChan
            if iChan ~= alignChan
                X2 = chanData{iChan}(:,:,islice);
                X2 = fourierShift( X2, diffphase, row_shift, col_shift );
                chanData{iChan}(:,:,islice) = X2;
            end
        end
        
        
        waitbar(islice / nslices, hwb);
    end
    
    % save the results
    setappdata(hChannelAxes(alignChan),'imageStack',Y);
    
    for iChan = 1:nChan
        if iChan ~= alignChan
            %setappdata(hChannelAxes(1),'imageStack',chanData{iChan});
            setappdata(hChannelAxes(iChan),'imageStack',chanData{iChan});
        end
    end
    
    
    close(hwb)
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

function mnuCalibrateScale( src, evt )
    [scaleX_volts_per_micron, scaleY_volts_per_micron] ...
        = puprisa_calibrateScale();
    
    % save these as application preferences
    if ~ispref('puprisa','scaleX_volts_per_micron')
        addpref('puprisa','scaleX_volts_per_micron',...
            scaleX_volts_per_micron);
        addpref('puprisa','scaleY_volts_per_micron',...
            scaleY_volts_per_micron);
    else
    
        setpref('puprisa','scaleX_volts_per_micron',scaleX_volts_per_micron);
        setpref('puprisa','scaleY_volts_per_micron',scaleY_volts_per_micron);
        
    end
    
end

function mnuLoadMultipleFiles( src, evt )
    
    %Exponential fitting of pixels in multiple files in a directory
    %Each delay stack fit parameters are saved individually
    
    clear all;
    
    fitParm = puprisa_ArtImaingDialog();
    
    % Parameters for multi-exponential fit
    minDelay = fitParm.MinDelay; %Min probe delay (ps)
    maxDelay = fitParm.MaxDelay; %Max probe delay (ps)
    nExp     = fitParm.NoOfExp; %No. of exponentials
    
    
    % check for valid inputs
    if( isempty(minDelay) || isempty(maxDelay) || isempty(nExp) )
        error('Invalid inputs to multiexponential fit');
    end
    
    %Initial conditions for multi-exponential fit
    A1 = fitParm.A1;
    T1 = fitParm.T1;
    A2 = fitParm.A2;
    T2 = fitParm.T2;
    
    %bounds for multi-exponential fit
    A1L = fitParm.A1L;
    T1L = fitParm.T1L;
    A2L = fitParm.A2L;
    T2L = fitParm.T2L;
    A1U = fitParm.A1U;
    T1U = fitParm.T1U;
    A2U = fitParm.A2U;
    T2U = fitParm.T2U;
    
    %fit options
    options = optimset('MaxIter',10000,'TolFun',1.0e-9,'MaxFunEvals',5000);
    
    %select only specific pixels with good peak Signal
    snrThresholdABS = fitParm.snrThresholdAbs;%absolute threshold
    
    %Enable binning(0/1)
    Bin = fitParm.bin;
    
    %type of fit positve/negative
    type = fitParm.type;
    
    
    %select directory
    directoryname = uigetdir('C:\Users\Prathyush\Data');
    Dir = dir(directoryname);
    numFiles = length(Dir);
    
    nameAdd = fitParm.fileName; 
    
    for i = 1:numFiles
        index = 1;
        nameLength = length(Dir(i).name);
        if(nameLength > 4)
            if(strcmp((Dir(i).name(end-3:end)),'.dat'))
                fileName = Dir(i).name;
                if ~isequal(fileName,0)
                    loadData(gcbf, [directoryname,'/',fileName]);
                    imageStack = getappdata(gcbf,'imageStackChannels');
                    if(Bin)
                        Y = binImage(imageStack{1});
                    else
                        Y = imageStack{1}; %select first channel
                    end                    
                    x = getappdata(gcbf,'delays');
                    mask = ((x >= minDelay) & (x <= maxDelay));
                    x = x(mask);
                    [m n ~] = size(Y);
                    if(type > 0)
                        peakYimage = (max(max(max(Y))));
                    elseif(type < 0)
                        peakYimage = (-1)*(min(min(min(Y))));
                        snrThresholdABS = snrThresholdABS*(-1);
                    end
                    snrThresholdLow = (fitParm.snrThresholdLow/100)*peakYimage;
                    snrThresholdHigh = (fitParm.snrThresholdHigh/100)*peakYimage;
                    
                    if(snrThresholdLow < snrThresholdABS)
                        snrThresholdLow = snrThresholdABS;
                    end
                    
                    for i2 = 1:m
                        for j = 1:n
                            y = (squeeze(Y(i2,j,:)))';
                            y = y(mask);
                            if (nExp == 1) 
                                fn = @(A,xx) A(1)*exp(-xx./A(2))+ A(3);
                                if(type > 0)
                                    peakI = max(y);
                                elseif(type < 0)
                                    peakI = ((-1)*min(y));
                                end
                                if (peakI > snrThresholdLow && peakI < snrThresholdHigh) %select only specific pixels with good Signal
                                    [AA,~,Residual,~,~,~,J] = lsqcurvefit(fn, [A1,T1,0.01],...
                                        x, finalDecayCurve,[A1L,T1L,-Inf],[A1U,T1U,Inf],options);
                                    ci = nlparci(AA,Residual,'jacobian',J);%95% confidence intervals
                                    AAerror = (ci(:,2)-ci(:,1))/(2*1.96);%standard error
                                    includePixel = 0;
                                    if(type > 0 && AA(1) > 0)%select only fits with positive signal
                                        includePixel = 1;
                                    elseif(type < 0 && AA(1) < 0)%select only fits with negative signal
                                        includePixel = 1;
                                    end
                                    if(includePixel)
                                        tauShort(index) = AA(2);
                                        ampShort(index) = AA(1);
                                        tauShortError(index) = AAerror(2);
                                        ampShortError(index) = AAerror(1);
                                        index= index+1;
                                    end
                                end
                            else
                                fn = @(A,xx) A(1)*exp(-xx./A(2)) + A(3)*exp(-xx./A(4))+A(5);
                                if(type > 0)
                                    peakI = max(y);
                                elseif(type < 0)
                                    peakI = ((-1)*min(y));
                                end
                                if (peakI > snrThresholdLow && peakI < snrThresholdHigh) %select only specific pixels with good Signal
                                    [AA,~,Residual,~,~,~,J] = lsqcurvefit(fn, [A1,T1,A2,T2,0.01],...
                                    x, y,[A1L,T1L,A2L,T2L,-Inf],[A1U,T1U,A2U,T2U,Inf],options);
                                    ci = nlparci(AA,Residual,'jacobian',J);%95% confidence intervals
                                    AAerror = (ci(:,2)-ci(:,1))/(2*1.96);%standard error
                                    includePixel = 0;
                                    if(type > 0 && AA(1) > 0 && AA(3) > 0)%only fits with positive signal
                                        includePixel = 1;
                                    elseif(type < 0 && AA(1) < 0 && AA(3) < 0) %select only fits with negative signal
                                        includePixel = 1;
                                    end
                                    if (includePixel)  
                                        if(AA(2) < AA(4))%assign to short or long component group
                                            tauLong(index) = AA(4);
                                            tauShort(index) = AA(2);
                                            ampLong(index) = AA(3);
                                            ampShort(index) = AA(1);
                                            tauLongError(index) = AAerror(4);
                                            tauShortError(index) = AAerror(2);
                                            ampLongError(index) = AAerror(3);
                                            ampShortError(index) = AAerror(1);
                                        else
                                            tauLong(index) = AA(2);
                                            tauShort(index) = AA(4);
                                            ampLong(index) = AA(1);
                                            ampShort(index) = AA(3);
                                            tauLongError(index) = AAerror(2);
                                            tauShortError(index) = AAerror(4);
                                            ampLongError(index) = AAerror(1);
                                            ampShortError(index) = AAerror(3);
                                        end
                                        index= index+1;
                                    end
                                end
                            end
                        end
                    end
                    
                    
                    FileName_S = [fileName(1:(end-4)),nameAdd,'_tauShort.txt'];
                    FileName_AS = [fileName(1:(end-4)),nameAdd,'_ampShort.txt'];
                    FileName_Serr = [fileName(1:(end-4)),nameAdd,'_tauShortError.txt'];
                    FileName_ASerr = [fileName(1:(end-4)),nameAdd,'_ampShortError.txt'];
                    
                    %write to file
                    fid = fopen(FileName_S, 'w');
                    fprintf(fid, '%12.6f\n', tauShort);
                    fclose(fid);
                    fid = fopen(FileName_AS, 'w');
                    fprintf(fid, '%12.6f\n', ampShort);
                    fclose(fid);
                    fid = fopen(FileName_Serr, 'w');
                    fprintf(fid, '%12.6f\n', tauShortError);
                    fclose(fid);
                    fid = fopen(FileName_ASerr, 'w');
                    fprintf(fid, '%12.6f\n', ampShortError);
                    fclose(fid);
                    
                    if (nExp > 1)
                        FileName_T = [fileName(1:(end-4)),nameAdd,'_tauLong.txt'];
                        FileName_AT = [fileName(1:(end-4)),nameAdd,'_ampLong.txt'];
                        FileName_Terr = [fileName(1:(end-4)),nameAdd,'_tauLongError.txt'];
                        FileName_ATerr = [fileName(1:(end-4)),nameAdd,'_ampLongError.txt'];
                        fid = fopen(FileName_T, 'w');
                        fprintf(fid, '%12.6f\n', tauLong);
                        fclose(fid);
                        fid = fopen(FileName_AT, 'w');
                        fprintf(fid, '%12.6f\n', ampLong);
                        fclose(fid);
                        fid = fopen(FileName_Terr, 'w');
                        fprintf(fid, '%12.6f\n', tauLongError);
                        fclose(fid);
                        fid = fopen(FileName_ATerr, 'w');
                        fprintf(fid, '%12.6f\n', ampLongError);
                        fclose(fid);
                        clear tauLongError tauLong ampLongError ampLong;
                    end
                    clear tauShort ampShort tauShortError ampShortError;
                end
            end
        end
    end
    
end

function mnuLoadMultipleFilesAVG( src, evt )
    

    %Exponentail fitting of all pixels in multiple files in a directory
    %Saved into one file to average over multiple delay stacks
    
    clear all;
    
    fitParm = puprisa_ArtImaingDialog();
    
    % Parameters for multi-exponential fit
    minDelay = fitParm.MinDelay; %Min probe delay (ps)
    maxDelay = fitParm.MaxDelay; %Max probe delay (ps)
    nExp     = fitParm.NoOfExp; %No. of exponentials
    
    
    % check for valid inputs
    if( isempty(minDelay) || isempty(maxDelay) || isempty(nExp) )
        error('Invalid inputs to multiexponential fit');
    end
    
    %Initial conditions for multi-exponential fit
    A1 = fitParm.A1;
    T1 = fitParm.T1;
    A2 = fitParm.A2;
    T2 = fitParm.T2;
    
    %bounds for multi-exponential fit
    A1L = fitParm.A1L;
    T1L = fitParm.T1L;
    A2L = fitParm.A2L;
    T2L = fitParm.T2L;
    A1U = fitParm.A1U;
    T1U = fitParm.T1U;
    A2U = fitParm.A2U;
    T2U = fitParm.T2U;
    
    %fit options
% changed MCF 2012-05-25: Added the Display->off option
		options = optimset('MaxIter',10000,'TolFun',1.0e-9,'MaxFunEvals',6000,'Display','off');
    
    %select only specific pixels with good peak Signal
    snrThresholdABS = fitParm.snrThresholdAbs;%absolute threshold
    
    %Enable binning(0/1)
    Bin = fitParm.bin;
    
    %type of fit positve/negative
    type = fitParm.type;
    
    %select directory
    directoryname = uigetdir('C:\Users\Prathyush\Data');
    Dir = dir(directoryname);
    numFiles = length(Dir);
    
    %filename to write data
    name = fitParm.fileName;

    index = 1;
    for i = 1:numFiles
        nameLength = length(Dir(i).name);
        if(nameLength > 4)
            if(strcmp((Dir(i).name(end-3:end)),'.dat'))
                fileName = Dir(i).name;
                if ~isequal(fileName,0)
% changed MCF 2012-05-25: display file name to monitor progress
                    disp(fileName);
                    loadData(gcbf, [directoryname,'/',fileName]);
                    imageStack = getappdata(gcbf,'imageStackChannels');
                    if(Bin)
                        Y = binImage(imageStack{1});
                    else 
                        Y = imageStack{1}; %select first channel
                    end
                    
                    x = getappdata(gcbf,'delays');
                    mask = ((x >= minDelay) & (x <= maxDelay));
                    x = x(mask);
                    [m n ~] = size(Y);
                    if(type > 0)
                        peakYimage = (max(max(max(Y))));
                    elseif(type < 0)
                        peakYimage = (-1)*(min(min(min(Y))));
                        snrThresholdABS = snrThresholdABS*(-1);
                    end
                    snrThresholdLow = (fitParm.snrThresholdLow/100)*peakYimage;
                    snrThresholdHigh = (fitParm.snrThresholdHigh/100)*peakYimage;
                    
                    if(snrThresholdLow < snrThresholdABS)
                        snrThresholdLow = snrThresholdABS;
                    end
                    
                    for i2 = 1:m
                        for j = 1:n
                            y = (squeeze(Y(i2,j,:)))';
                            y = y(mask);
                            if (nExp == 1)
                                fn = @(A,xx) A(1)*exp(-xx./A(2))+ A(3);
                                if(type > 0)
                                    peakI = max(y);
                                elseif(type < 0)
                                    peakI = ((-1)*min(y));
                                end
                                if (peakI > snrThresholdLow && peakI < snrThresholdHigh) %select only specific pixels with good Signal
                                    [AA,~,Residual,~,~,~,J] = lsqcurvefit(fn, [A1,T1,0.01],...
                                        x, finalDecayCurve,[A1L,T1L,-Inf],[A1U,T1U,Inf],options);
                                    ci = nlparci(AA,Residual,'jacobian',J);%95% confidence intervals
                                    AAerror = (ci(:,2)-ci(:,1))/(2*1.96);%standard error
                                    includePixel = 0;
                                    if(type > 0 && AA(1) > 0)%select only fits with positive signal
                                        includePixel = 1;
                                    elseif(type < 0 && AA(1) < 0)%select only fits with negative signal
                                        includePixel = 1;
                                    end
                                    if(includePixel)
                                        tauShort(index) = AA(2);
                                        ampShort(index) = AA(1);
                                        tauShortError(index) = AAerror(2);
                                        ampShortError(index) = AAerror(1);
                                        index= index+1;
                                    end
                                end
                            else
                                fn = @(A,xx) A(1)*exp(-xx./A(2)) + A(3)*exp(-xx./A(4))+A(5);
                                if(type > 0)
                                    peakI = max(y);
                                elseif(type < 0)
                                    peakI = ((-1)*min(y));
                                end
                                if (peakI > snrThresholdLow && peakI < snrThresholdHigh) %select only specific pixels with good Signal
                                    [AA,~,Residual,~,~,~,J] = lsqcurvefit(fn, [A1,T1,A2,T2,0.01],...
                                        x, y,[A1L,T1L,A2L,T2L,-Inf],[A1U,T1U,A2U,T2U,Inf],options);
                                    ci = nlparci(AA,Residual,'jacobian',J);%95% confidence intervals
                                    AAerror = (ci(:,2)-ci(:,1))/(2*1.96);%standard error
                                    includePixel = 0;
                                    if(type > 0 && AA(1) > 0 && AA(3) > 0)%only fits with positive signal
                                        includePixel = 1;
                                    elseif(type < 0 && AA(1) < 0 && AA(3) < 0) %select only fits with negative signal
                                        includePixel = 1;
                                    end
                                    if(includePixel)
                                        if(AA(2) < AA(4))%assign to short or long component group
                                            tauLong(index) = AA(4);
                                            tauShort(index) = AA(2);
                                            ampLong(index) = AA(3);
                                            ampShort(index) = AA(1);
                                            tauLongError(index) = AAerror(4);
                                            tauShortError(index) = AAerror(2);
                                            ampLongError(index) = AAerror(3);
                                            ampShortError(index) = AAerror(1);
                                        else
                                            tauLong(index) = AA(2);
                                            tauShort(index) = AA(4);
                                            ampLong(index) = AA(1);
                                            ampShort(index) = AA(3);
                                            tauLongError(index) = AAerror(2);
                                            tauShortError(index) = AAerror(4);
                                            ampLongError(index) = AAerror(1);
                                            ampShortError(index) = AAerror(3);
                                        end
% changed MCF 2012-05-25: log intensity
																				intens(index) = sqrt(sum(y.^2));
                                        index= index+1;
                                    end
                                end
                            end
                            
                        end
                    end
                end
            end
        end
    end
    
    FileName_S = [name,'_tauShort.txt'];
    FileName_AS = [name,'_ampShort.txt'];
    FileName_Serr = [name,'_tauShortError.txt'];
    FileName_ASerr = [name,'_ampShortError.txt'];

		
    %write to file
    fid = fopen(FileName_S, 'w');
    fprintf(fid, '%12.6f\n', tauShort);
    fclose(fid);
    fid = fopen(FileName_AS, 'w');
    fprintf(fid, '%12.6f\n', ampShort);
    fclose(fid);
    fid = fopen(FileName_Serr, 'w');
    fprintf(fid, '%12.6f\n', tauShortError);
    fclose(fid); 
    fid = fopen(FileName_ASerr, 'w');
    fprintf(fid, '%12.6f\n', ampShortError);
    fclose(fid);
    
    if (nExp > 1)
        FileName_L = [name,'_tauLong.txt'];
        FileName_AL = [name,'_ampLong.txt'];
        FileName_Lerr = [name,'_tauLongError.txt'];
        FileName_ALerr = [name,'_ampLongError.txt'];
        fid = fopen(FileName_L, 'w');
        fprintf(fid, '%12.6f\n', tauLong);
        fclose(fid);
        fid = fopen(FileName_AL, 'w');
        fprintf(fid, '%12.6f\n', ampLong);
        fclose(fid);
        fid = fopen(FileName_Lerr, 'w');
        fprintf(fid, '%12.6f\n', tauLongError);
        fclose(fid);
        fid = fopen(FileName_ALerr, 'w');
        fprintf(fid, '%12.6f\n', ampLongError);
        fclose(fid);

% changed MCF 2012-05-25: write combo file (including intensity)
				FileName_combo = [name,'_combo.txt'];
				temp=[intens' tauShort' ampShort' tauLong' ampLong'];
				dlmwrite(FileName_combo, temp, '\t');	
		end
end


function mnuLoadMultipleFilesAVG2( src, evt )

%Average all pixels in multiple files in a directory
%And save the fit parameters into an excel file
%Used for averaging over multiple delay stacks

clear all;

fitParm = puprisa_ArtImaingDialog();

% Parameters for multi-exponential fit
minDelay = fitParm.MinDelay; %Min probe delay (ps)
maxDelay = fitParm.MaxDelay; %Max probe delay (ps)
nExp     = fitParm.NoOfExp; %No. of exponentials


% check for valid inputs
if( isempty(minDelay) || isempty(maxDelay) || isempty(nExp) )
    error('Invalid inputs to multiexponential fit');
end

%Initial conditions for multi-exponential fit
A1 = fitParm.A1;
T1 = fitParm.T1;
A2 = fitParm.A2;
T2 = fitParm.T2;

%bounds for multi-exponential fit
A1L = fitParm.A1L;
T1L = fitParm.T1L;
A2L = fitParm.A2L;
T2L = fitParm.T2L;
A1U = fitParm.A1U;
T1U = fitParm.T1U;
A2U = fitParm.A2U;
T2U = fitParm.T2U;

%fit options
options = optimset('MaxIter',20000,'TolFun',1.0e-9,'MaxFunEvals',20000);

%select only specific pixels with good peak Signal
snrThresholdABS = fitParm.snrThresholdAbs;%absolute threshold

%Enable binning(0/1)
Bin = fitParm.bin;

%type of fit positve/negative
type = fitParm.type;

%select directory
directoryname = uigetdir('C:\Users\Prathyush\Data');
Dir = dir(directoryname);
numFiles = length(Dir);

%filename to write data
name = fitParm.fileName;
avgY = 0;
index = 0;

for i = 1:numFiles
    nameLength = length(Dir(i).name);
    if(nameLength > 4)
        if(strcmp((Dir(i).name(end-3:end)),'.dat'))
            fileName = Dir(i).name;
            if ~isequal(fileName,0)
                loadData(gcbf, [directoryname,'/',fileName]);
                imageStack = getappdata(gcbf,'imageStackChannels');
                if(Bin)
                    Y = binImage(imageStack{1});
                else
                    Y = imageStack{1}; %select first channel
                end
                x = getappdata(gcbf,'delays');
                mask = ((x >= minDelay) & (x <= maxDelay));
                x = x(mask);
                [m n ~] = size(Y);
                if(type > 0)
                    peakYimage = (max(max(max(Y))));
                elseif(type < 0)
                    peakYimage = (-1)*(min(min(min(Y))));
                    snrThresholdABS = snrThresholdABS*(-1);
                end
                snrThresholdLow = (fitParm.snrThresholdLow/100)*peakYimage;
                snrThresholdHigh = (fitParm.snrThresholdHigh/100)*peakYimage;
                
                if(snrThresholdLow < snrThresholdABS)
                    snrThresholdLow = snrThresholdABS;
                end
                
                for i2 = 1:m
                    for j = 1:n
                        y = (squeeze(Y(i2,j,:)))';
                        y = y(mask);
                        if(type > 0)
                            peakI = max(y);
                        elseif(type < 0)
                            peakI = ((-1)*min(y));
                        end
                       
                        %add up all the decays
                        if (peakI > snrThresholdLow && peakI < snrThresholdHigh)%select only specific pixels with restricted Signal
                            avgY = avgY + y;
                            index= index+1;
                        end
                    end
                end
            end
        end
    end
end

%averaged decay curve
finalDecayCurve = avgY/index;
               
if (nExp == 1)
    fn = @(A,xx) A(1)*exp(-xx./A(2))+ A(3);
    [AA,~,Residual,~,~,~,J] = lsqcurvefit(fn, [A1,T1,0.01],...
        x, finalDecayCurve,[A1L,T1L,-Inf],[A1U,T1U,Inf],options);
    ci = nlparci(AA,Residual,'jacobian',J);%95% confidence intervals
    AAerror = (ci(:,2)-ci(:,1))/(2*1.96);
    tauShort = AA(2);
    ampShort = AA(1);
    tauShortError = AAerror(2);
    ampShortError = AAerror(1);
    baseline = AA(3);
    excelWrite = {'TimeConstShort','AmpShort','Baseline','TimeConstShortError'...
    ,'AmpShortError';tauShort,ampShort,baseline,tauShortError,ampShortError};
else
    fn = @(A,xx) A(1)*exp(-xx./A(2)) + A(3)*exp(-xx./A(4))+A(5);
    [AA,~,Residual,~,~,~,J] = lsqcurvefit(fn, [A1,T1,A2,T2,0.01],...
        x, finalDecayCurve,[A1L,T1L,A2L,T2L,-Inf],[A1U,T1U,A2U,T2U,Inf],options);
    ci = nlparci(AA,Residual,'jacobian',J);%95% confidence intervals
    AAerror = (ci(:,2)-ci(:,1))/(2*1.96);
    
    if(AA(2) < AA(4))%assign to short or long component group
        tauLong = AA(4);
        tauShort = AA(2);
        ampLong = AA(3);
        ampShort = AA(1);
        tauLongError = AAerror(4);
        tauShortError = AAerror(2);
        ampLongError = AAerror(3);
        ampShortError = AAerror(1);
    else
        tauLong = AA(2);
        tauShort = AA(4);
        ampLong = AA(1);
        ampShort = AA(3);
        tauLongError = AAerror(2);
        tauShortError = AAerror(4);
        ampLongError = AAerror(1);
        ampShortError = AAerror(3);
    end
    baseline = AA(5);
    excelWrite = {'TimeConstShort','AmpShort','TimeConstLong','AmpLong','Baseline','TimeConstShortError'...
    ,'AmpShortError','TimeConstLongError','AmpLongError';tauShort,ampShort,tauLong...
    ,ampLong,baseline,tauShortError,ampShortError,tauLongError,ampLongError};
end

FileName = [name,'.xlsx'];
xlswrite(FileName,excelWrite);
 
end

function mnuLoadMultipleFilesAVG3( src, evt )

%Average all pixels in an image and do exponential fitting
%For multiple files in a directory and save fit parameters into a
%single excel file

clear all;
fitParm = puprisa_ArtImaingDialog();

% Parameters for multi-exponential fit
minDelay = fitParm.MinDelay; %Min probe delay (ps)
maxDelay = fitParm.MaxDelay; %Max probe delay (ps)
nExp     = fitParm.NoOfExp; %No. of exponentials


% check for valid inputs
if( isempty(minDelay) || isempty(maxDelay) || isempty(nExp) )
    error('Invalid inputs to multiexponential fit');
end

%Initial conditions for multi-exponential fit
A1 = fitParm.A1;
T1 = fitParm.T1;
A2 = fitParm.A2;
T2 = fitParm.T2;

%bounds for multi-exponential fit
A1L = fitParm.A1L;
T1L = fitParm.T1L;
A2L = fitParm.A2L;
T2L = fitParm.T2L;
A1U = fitParm.A1U;
T1U = fitParm.T1U;
A2U = fitParm.A2U;
T2U = fitParm.T2U;


%fit options
options = optimset('MaxIter',5000,'TolFun',1.0e-8,'MaxFunEvals',5000);

%select only specific pixels with good peak Signal
snrThresholdABS = fitParm.snrThresholdAbs;%absolute threshold

%Enable binning(0/1)
Bin = fitParm.bin;

%type of fit positve/negative
type = fitParm.type;

%select directory
directoryname = uigetdir('C:\Users\Prathyush\Data');
Dir = dir(directoryname);
numFiles = length(Dir);

%filename to write data
name = fitParm.fileName;

if (nExp == 1)
    excelWrite = {'File Name','TimeConstShort','AmpShort','Baseline'...
        ,'TimeConstShortError','AmpShortError'};
else
    excelWrite = {'File Name','TimeConstShort','AmpShort','TimeConstLong','AmpLong','Baseline'...
        ,'TimeConstShortError','AmpShortError','TimeConstLongError','AmpLongError'};
end
fileindex = 2;

for i = 1:numFiles
    avgY = 0;
    index = 0;
    nameLength = length(Dir(i).name);
    if(nameLength > 4)
        if(strcmp((Dir(i).name(end-3:end)),'.dat'))
            fileName = Dir(i).name;
            if ~isequal(fileName,0)
                loadData(gcbf, [directoryname,'/',fileName]);
                imageStack = getappdata(gcbf,'imageStackChannels');
                if(Bin)
                    Y = binImage(imageStack{1});
                else
                    Y = imageStack{1}; %select first channel
                end
                x = getappdata(gcbf,'delays');
                mask = ((x >= minDelay) & (x <= maxDelay));
                x = x(mask);
                
                [m n ~] = size(Y);
                if(type > 0)
                    peakYimage = (max(max(max(Y))));
                elseif(type < 0)
                    peakYimage = (-1)*(min(min(min(Y))));
                    snrThresholdABS = snrThresholdABS*(-1);
                end
                snrThresholdLow = (fitParm.snrThresholdLow/100)*peakYimage;
                snrThresholdHigh = (fitParm.snrThresholdHigh/100)*peakYimage;
                
                if(snrThresholdLow < snrThresholdABS)
                    snrThresholdLow = snrThresholdABS;
                end
                
                for i2 = 1:m
                    for j = 1:n
                        y = (squeeze(Y(i2,j,:)))';
                        y = y(mask);
                        if(type > 0)
                            peakI = max(y);
                        elseif(type < 0)
                            peakI = ((-1)*min(y));
                        end
                        
                        %add up all the decays
                        if (peakI > snrThresholdLow && peakI < snrThresholdHigh) %select only specific pixels with good Signal
                            avgY = avgY + y;
                            index= index+1;
                        end
                    end
                end
                finalDecayCurve = avgY/index;
                if (nExp == 1)
                    fn = @(A,xx) A(1)*exp(-xx./A(2))+ A(3);
                    [AA,~,Residual,~,~,~,J] = lsqcurvefit(fn, [A1,T1,0.01],...
                        x, finalDecayCurve,[A1L,T1L,-Inf],[A1U,T1U,Inf],options);
                    ci = nlparci(AA,Residual,'jacobian',J);%95% confidence intervals
                    AAerror = (ci(:,2)-ci(:,1))/(2*1.96);
                    tauShort = AA(2);
                    ampShort = AA(1);
                    tauShortError = AAerror(2);
                    ampShortError = AAerror(1);
                    baseline = AA(3);
                    excelWrite(fileindex,:) = {fileName,tauShort,ampShort...
                        ,baseline,tauShortError,ampShortError};
                else
                    fn = @(A,xx) A(1)*exp(-xx./A(2)) + A(3)*exp(-xx./A(4))+A(5);
                    [AA,~,Residual,~,~,~,J] = lsqcurvefit(fn, [A1,T1,A2,T2,0.01],...
                        x, finalDecayCurve,[A1L,T1L,A2L,T2L,-Inf],[A1U,T1U,A2U,T2U,Inf],options);
                    ci = nlparci(AA,Residual,'jacobian',J);%95% confidence intervals
                    AAerror = (ci(:,2)-ci(:,1))/(2*1.96);%standard error
                    
                    if(AA(2) < AA(4))%assign to short or long component group
                        tauLong = AA(4);
                        tauShort = AA(2);
                        ampLong = AA(3);
                        ampShort = AA(1);
                        tauLongError = AAerror(4);
                        tauShortError = AAerror(2);
                        ampLongError = AAerror(3);
                        ampShortError = AAerror(1);
                    else
                        tauLong = AA(2);
                        tauShort = AA(4);
                        ampLong = AA(1);
                        ampShort = AA(3);
                        tauLongError = AAerror(2);
                        tauShortError = AAerror(4);
                        ampLongError = AAerror(1);
                        ampShortError = AAerror(3);
                    end
                    baseline = AA(5);
                    excelWrite(fileindex,:) = {fileName,tauShort,ampShort,tauLong...
                        ,ampLong,baseline,tauShortError,ampShortError,tauLongError,ampLongError};
                end
                fileindex = fileindex+1;
            end
        end
    end
end
writefileName = [name,'.xlsx'];
xlswrite(writefileName,excelWrite);

end

function mnuLoadMultipleFilesSameArea( src, evt )

%Average all pixels in selected areas of the image for each file in a
%directory
%Areas are selected by the thresholding the selected file
%single excel file to store data of all  files in the directory

clear all;

% Parameters for multi-exponential fit
minDelay = 0.21; %Min probe delay (ps)
maxDelay = Inf; %Max probe delay (ps)
nExp     = 2; %No. of exponentials

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

%fit options
options = optimset('MaxIter',5000,'TolFun',1.0e-8,'MaxFunEvals',5000);

%select only specific pixels with good peak Signal
snrThresholdABS = 0;

%Enable binning(0/1)
Bin = 0;

%select directory
[fileName, pathName] = uigetfile(['C:\Users\Prathyush\Data','/*.dat'],...
    'Select a reference imageStack' );
Dir = dir(pathName);
numFiles = length(Dir);

%filename to write data
name = inputdlg({'Save file name'},'File Name',1,{'test'});

if ~isequal(fileName,0)
    loadData(gcbf, [pathName,'/',fileName]);
    imageStack = getappdata(gcbf,'imageStackChannels');
    if(Bin)
        Y = binImage(imageStack{1});
    else
        Y = imageStack{1}; %select first channel
    end
    x = getappdata(gcbf,'delays');
    mask = ((x >= minDelay) & (x <= maxDelay));
    [m n ~] = size(Y);
    MinYimage = (min(min(min(Y))));
    snrThresholdLow = 0.03*MinYimage; %threshold as a % of peak value
    snrThresholdHigh = 1*MinYimage; %
    if(snrThresholdLow > snrThresholdABS)
        snrThresholdLow = snrThresholdABS;
    end
    for i2 = 1:m
        for j = 1:n
            y = (squeeze(Y(i2,j,:)))';
            y = y(mask);
            [minI ~] = min(y);
            if (minI < snrThresholdLow && minI > snrThresholdHigh) %select only specific pixels with good Signal
                RefImageMask(i2,j) = 1;
            else
                RefImageMask(i2,j) = 0;% reject the pixel
            end
        end
    end
end

excelWrite = {'File Name','Decay Const(Short)','AmpShort','Decay Const(Long)','AmpLong','Fit Error'};
fileindex = 2;

for i = 1:numFiles
    avgY = 0;
    index = 0;
    nameLength = length(Dir(i).name);
    if(nameLength > 4)
        if(strcmp((Dir(i).name(end-3:end)),'.dat'))
            fileName = Dir(i).name;
            if ~isequal(fileName,0)
                loadData(gcbf, [pathName,'/',fileName]);
                imageStack = getappdata(gcbf,'imageStackChannels');
                if(Bin)
                    Y = binImage(imageStack{1});
                else
                    Y = imageStack{1}; %select first channel
                end
                x = getappdata(gcbf,'delays');
                mask = ((x >= minDelay) & (x <= maxDelay));
                x = x(mask);
                [dummy sizeOFx] = size(x);
                [m n ~] = size(Y);
                for i2 = 1:m
                    for j = 1:n
                        if(RefImageMask(i2,j))
                            y = (squeeze(Y(i2,j,:)))';
                            y = y(mask);
                            [minI ~] = min(y);
                            %add up all the decays
                            if (minI < 0) %select only negative pixels
                                avgY = avgY + y;
                                index= index+1;
                            end
                        end
                    end
                end
                finalDecayCurve = avgY/index;
                if (nExp == 1)
                    fn = @(A,xx) A(1)*exp(-xx./A(2))+ A(3);
                    [AA,Residual] = lsqcurvefit(fn, [A1,T1,0.01], x, finalDecayCurve,[],[],options);
                    tauShort = AA(2);
                    ampShort = AA(1);
                    fitError = sqrt(Residual/sizeOFx);
                else
                    fn = @(A,xx) A(1)*exp(-xx./A(2)) + A(3)*exp(-xx./A(4))+A(5);
                    [AA,Residual] = lsqcurvefit(fn, [A1,T1,A2,T2,0.01], x, finalDecayCurve,[],[],options);
                    if(AA(2) < AA(4))%assign to short or long component group
                        tauLong = AA(4);
                        tauShort = AA(2);
                        ampLong = AA(3);
                        ampShort = AA(1);
                    else
                        tauLong = AA(2);
                        tauShort = AA(4);
                        ampLong = AA(1);
                        ampShort = AA(3);
                    end
                    fitError = sqrt(Residual/sizeOFx);
                end
                excelWrite(fileindex,:) = {fileName,tauShort,ampShort,tauLong,ampLong,fitError};
                fileindex = fileindex+1;
            end
        end
    end
end
writefileName = [name{1},'.xlsx'];
xlswrite(writefileName,excelWrite);
end

function mnuLoadMultipleFilesRatio( src, evt )

%calculate the number of +ve and -ve pixels in an image
%over multiple files ina directory

clear all;

%Enable binning(0/1)
Bin = 0;

%select directory
directoryname = uigetdir('C:\Users\Prathyush\Data');
Dir = dir(directoryname);
numFiles = length(Dir);

%filename to write data

NegativePix = 0;
PositivePix = 0;
TotalPix = 0;
for i = 1:numFiles
    nameLength = length(Dir(i).name);
    if(nameLength > 4)
        if(strcmp((Dir(i).name(end-3:end)),'.dat'))
            fileName = Dir(i).name;
            if ~isequal(fileName,0)
                loadData(gcbf, [directoryname,'/',fileName]);
                imageStack = getappdata(gcbf,'imageStackChannels');
                if(Bin)
                    Y = binImage(imageStack{1});
                else
                    Y = imageStack{1}; %select first channel
                end
                x = getappdata(gcbf,'delays');
                mask = ((x >= 0) & (x <= 1));
                [m n ~] = size(Y);
                ThresholdPositive = 0.05;
                ThresholdNegative = -0.08;
               
                for i2 = 1:m
                    for j = 1:n
                        y = (squeeze(Y(i2,j,:)))';
                        y = y(mask);
                        [minI ~] = min(y);
                        [maxI ~] = max(y);
                        TotalPix = TotalPix+1;
                        if (minI < ThresholdNegative)%select only specific pixels with good Signal
                            NegativePix = NegativePix+1;
                        elseif(maxI > ThresholdPositive)
                            PositivePix = PositivePix+1;
                        end
                    end
                end
            end
        end
    end
end

Pigment = (NegativePix*100)/TotalPix
Impurity = (PositivePix*100)/TotalPix

end

function mnuLoadMultipleFilesAVG2new( src, evt )

%Average all pixels in multiple files in a directory
%And save the fit parameters into an excel file
%Used for averaging over multiple delay stacks

%NEW indicates use of fit function with gaussian convolution

clear all;

fitParm = puprisa_ArtImaingDialogNEW();

nExp  = fitParm.NoOfExp; %No. of exponentials
xfwhm = fitParm.CCpulsewidth/1000;%cross-correlation width in ps
x0 = fitParm.DelayOffset/1000;%pump-probe delay offset in ps

if (nExp == 1)
    LB = [fitParm.InstRespL,fitParm.A1L,fitParm.T1L,-Inf];
    UB = [fitParm.InstRespU,fitParm.A1U,fitParm.T1U,Inf];
    InCon = [fitParm.InstResponse,fitParm.A1,fitParm.T1,0];
else
    %bounds for multi-exponential fit
    LB = [fitParm.InstRespL,fitParm.A1L,fitParm.T1L,fitParm.A2L,fitParm.T2L,-Inf];
    UB = [fitParm.InstRespU,fitParm.A1U,fitParm.T1U,fitParm.A2U,fitParm.T2U,Inf];
    %Initial conditions for multi-exponential fit
    InCon = [fitParm.InstResponse,fitParm.A1,fitParm.T1,fitParm.A2,fitParm.T2,0];
end


%fit options
options = optimset('MaxIter',20000,'TolFun',1.0e-9,'MaxFunEvals',20000);

%select only specific pixels with good peak Signal
snrThresholdABS = fitParm.snrThresholdAbs;%absolute threshold

%Enable binning(0/1)
Bin = fitParm.bin;

%type of fit positve/negative
type = fitParm.type;

%select directory
directoryname = uigetdir('C:\Users\Prathyush\Data');
Dir = dir(directoryname);
numFiles = length(Dir);

%filename to write data
name = fitParm.fileName;
avgY = 0;
index = 0;

for i = 1:numFiles
    nameLength = length(Dir(i).name);
    if(nameLength > 4)
        if(strcmp((Dir(i).name(end-3:end)),'.dat'))
            fileName = Dir(i).name;
            if ~isequal(fileName,0)
                loadData(gcbf, [directoryname,'/',fileName]);
                imageStack = getappdata(gcbf,'imageStackChannels');
                if(Bin)
                    Y = binImage(imageStack{1});
                else
                    Y = imageStack{1}; %select first channel
                end
                x = getappdata(gcbf,'delays');
                [m n ~] = size(Y);
                if(type > 0)
                    peakYimage = (max(max(max(Y))));
                elseif(type < 0)
                    peakYimage = (-1)*(min(min(min(Y))));
                    snrThresholdABS = snrThresholdABS*(-1);
                end
                snrThresholdLow = (fitParm.snrThresholdLow/100)*peakYimage;
                snrThresholdHigh = (fitParm.snrThresholdHigh/100)*peakYimage;
                
                if(snrThresholdLow < snrThresholdABS)
                    snrThresholdLow = snrThresholdABS;
                end
                
                for i2 = 1:m
                    for j = 1:n
                        y = (squeeze(Y(i2,j,:)))';
                        if(type > 0)
                            peakI = max(y);
                        elseif(type < 0)
                            peakI = ((-1)*min(y));
                        end
                        
                        %add up all the decays
                        if (peakI > snrThresholdLow && peakI < snrThresholdHigh)%select only specific pixels with restricted Signal
                            avgY = avgY + y;
                            index= index+1;
                        end
                    end
                end
            end
        end
    end
end

%averaged decay curve
finalDecayCurve = avgY/index;

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
    [AA,~,Residual,~,~,~,J] = lsqcurvefit(fn,InCon,x,finalDecayCurve,LB,UB,options);
    ci = nlparci(AA,Residual,'jacobian',J);%95% confidence intervals
    AAerror = (ci(:,2)-ci(:,1))/(2*1.96); %standard error
    tau1 = AA(3);
    amp1 = AA(2);
    tau1Error = AAerror(3);
    amp1Error = AAerror(2);
    InstResponse = AA(1);
    excelWrite = {'TimeConst','Amplitude','InstResponse','TimeConstError'...
        ,'AmpError';tau1,amp1,InstResponse,tau1Error,amp1Error};
else
    fn = @(A,x) multiExpFitFunc(A,x,xNew,I,delFunc,stepFunc);
    [AA,~,Residual,~,~,~,J] = lsqcurvefit(fn,InCon,x,finalDecayCurve,LB,UB,options);
    ci = nlparci(AA,Residual,'jacobian',J);%95% confidence intervals
    AAerror = (ci(:,2)-ci(:,1))/(2*1.96); %standard error
        
    tau2 = AA(5);
    tau1 = AA(3);
    amp2 = AA(4);
    amp1 = AA(2);
    tau2Error = AAerror(5);
    tau1Error = AAerror(3);
    amp2Error = AAerror(4);
    amp1Error = AAerror(2);
    InstResponse = AA(1);
    
    excelWrite = {'TimeConst1','Amp1','TimeConst2','Amp2','InstResponse','TimeConst1Error'...
    ,'Amp1Error','TimeConst2Error','Amp2Error';tau1,amp1,tau2...
    ,amp2,InstResponse,tau1Error,amp1Error,tau2Error,amp2Error};
end

FileName = [name,'.xlsx'];
xlswrite(FileName,excelWrite');

figure;
ee = fn(AA,x);
plot(x,finalDecayCurve,x,ee,'-r');

end

function y = SingleExpFitFunc(A,t,tNew,I,delFunc,stepFunc)

yy = A(1)*delFunc...
    + (A(2)*exp(-(tNew.*stepFunc)/A(3))).*stepFunc;

yyConv = (yy*I);

y = interp1(tNew,yyConv,t,'spline')+ A(4);

end

function y = multiExpFitFunc(A,t,tNew,I,delFunc,stepFunc)

yy = A(1)*delFunc...
    + (A(2)*exp(-(tNew.*stepFunc)/A(3))...
    + A(4)*exp(-(tNew.*stepFunc)/A(5))).*stepFunc;

yyConv = (yy*I);

y = interp1(tNew,yyConv,t,'spline')+ A(6);

end

%--------------------------------------------------------------------------
% CALLBACK SUPPORT FUNCTIONS
%


function X = fourierShift( X, diffphase, row_shift, col_shift )
    % copied from dft registration code
    [nr,nc]=size(X);
    Nr = ifftshift([-fix(nr/2):ceil(nr/2)-1]);
    Nc = ifftshift([-fix(nc/2):ceil(nc/2)-1]);
    [Nc,Nr] = meshgrid(Nc,Nr);
    Greg = fft2(X).*exp(i*2*pi*(-row_shift*Nr/nr-col_shift*Nc/nc));
    Greg = Greg*exp(i*diffphase);
    X = real(ifft2(Greg));
end

function setChannelAxesPos( f )
    % get relevant objects
    hChannelAxes = getappdata(f,'hChannelAxes');
    hStatusText = findobj(f,'tag','hStatusText');
    hHeaderText= findobj(f,'tag','hHeaderText');
    
    % set statusbar position
    statusTextPos = get(hStatusText,'position');
    statusTextHeight = statusTextPos(4);
    headerTextWidth = 200;
    
    % assume figure units are pixels
    pos = get(f,'Position');
    figWidth = pos(3) - headerTextWidth;
    figHeight = pos(4) - statusTextHeight;
    bottom = statusTextHeight;
    
    % set header text position
    headerTextPos = get(hHeaderText,'position');
    headerTextPos(4) = figHeight;
    headerTextPos(2) = statusTextHeight;
    headerTextPos(1) = pos(3) - headerTextWidth;
    set(hHeaderText,'position',headerTextPos);
    
    
    nChanRow = getappdata(f,'nChanRow');
    nChanCol = getappdata(f,'nChanRow');

  
    
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
% Load data, then parse into figure appdata structures, and call update
% routines to display the newly loaded stack
%
    setappdata(f,'fileName',fileName);

    % set cursor to 'watch' to indicate to the user we're busy
    set(f,'pointer','watch');
    drawnow;

    % load the image stack data
    [slices, header] = puprisa_readImageStack( fileName);
    nSlices = length(slices);
    
    % exit this function if the stack failed to load
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
    
    % parse loaded data into one image stack per channel
    imageStackThisChannel = zeros(nRows, nCols);
    for iChannel = 1:nChannels
        % pull out each channel
        for iSlice = 1:nSlices
            imageStackThisChannel(:,:, iSlice) ...
                = slices(iSlice).imageData{iChannel};
        end
        imageStackChannels{iChannel} = imageStackThisChannel;
    end
        
    % store the image stacks in this figure's appdata
    setappdata(f,'imageStackChannels',imageStackChannels);
    
    % store delay axes in this figure's appdata
    setappdata(f, 'delays', cell2mat({slices.delays}));
    setappdata(f, 'xPos', cell2mat({slices.posX}));
    setappdata(f, 'yPos', cell2mat({slices.posY}));
    setappdata(f, 'zPos', cell2mat({slices.posZ}));
    
    % determine whether this was a delay stack, z stack, or mosaic
    stackType = '';
    
    % first see if it was a mosaic
    if isfield( header, 'mosaic' )
        if header.mosaic == 1
            stackType = 'mosaic';
        end
    end
    
    % if not a mosaic, then see if it's a z stack or a delay stack
    if isempty(stackType)
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
    updateAll( f );
    set(f,'pointer','arrow');
end

function mnuBiexpTestSet( src, evt )
    % come up with a 32x32 test set
    
    nrow = 32;
    ncol =32;
    nt = 100;
    t = linspace(0,200,nt);
    
%     % generate biexp. decays for the 4 quadrants
%     %[A10,T10,A20,T20,off0]
%     A_ul = [.5,3,.5,0.2,0];
%     A_ll = [0.5,3,-0.5,0.2,0];
%     A_ur = [-2.0,5, 1.0, 0.5, 0];
%     A_lr = [-0.1, 0.4, 1.0, 6, 0];
%     
%     fn = @(A,xx) A(1)*exp(-xx./A(2)) + A(3)*exp(-xx./A(4))+A(5);
%     
%     f_ul = fn(A_ul, t);
%     f_ll = fn(A_ll, t);
%     f_ur = fn(A_ur, t);
%     f_lr = fn(A_lr, t);
%     
%     % make an image stack
%     X(1:128,1:128,:) = repmat(reshape(f_ul,[1,1,nt]),[128,128,1]);
%     X(129:256,1:128,:) = repmat(reshape(f_ll,[1,1,nt]),[128,128,1]);
%     X(1:128,129:256,:) = repmat(reshape(f_ur,[1,1,nt]),[128,128,1]);
%     X(129:256,129:256,:) = repmat(reshape(f_lr,[1,1,nt]),[128,128,1]);


T1=[0.6 1 5 35];
T2=[4 6 10 60];
A1=[-0.3 -0.3 0.4 0.8];
A2=[-0.1 -0.05 0.1 0.2];
grid=32;

for i= 1:2
    for j = 1:2
            fn=A1(2*(i-1)+j)*exp(-t/T1(2*(i-1)+j))+A2(2*(i-1)+j)*exp(-t/T2(2*(i-1)+j));
            R=reshape(fn,[1,1,nt]);
            blk=repmat(R, [grid,grid,1]);
            X(1+(i-1)*grid:i*grid,1+(j-1)*grid:j*grid,:)=blk;
    end
end



          X=repmat(X, [256/(grid*2),256/(grid*2),1])+0.01*randn(256,256,nt);

    
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
    setappdata(gcbf, 'fileName', 'biexptestset');
    
    hChannelImages = getappdata(gcbf,'hChannelImages');
    hChannelAxes = getappdata(gcbf,'hChannelAxes');
    
    setappdata(hChannelAxes(1),'imageStack', X);
    
    set(hChannelImages(1),'XData',[1,ncol],'YData',[1,nrow],...
            'CData',squeeze(X(:,:,1)));
    set(hChannelAxes(1),'XLim',[1,ncol],'YLim',[1,nrow]);

    
    % display all channels for current slice
    updateAll( gcbf );
end

function Xnew = binImage(X)
    % bin image, lowering resolution to improve SNR
    % X is the image stack
    binFactor = 2;
      
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

    nChannels = getappdata(f,'nChannels');
    for iChannel = 1:nChannels
        imageStack = getappdata(hChannelAxes(iChannel),'imageStack');
        imageSlice = squeeze(imageStack(:,:,currentSlice));
        
        set(hChannelImages(iChannel),'CData',imageSlice);
    end
    
    fileName = getappdata(f,'fileName');
    updateStatus([fileName,': slice ', num2str(currentSlice), ' / ', ...
        num2str(nSlices)]);
    
    % update header display
    % EASIEST WAY TO FIX THIS:
    % change readImageStack to retain a header string,
    % in addition to the struct
    hHeaderText = findobj(f,'tag','hHeaderText');
    header = getappdata(f,'fileHeader');
    set(hHeaderText,'string',header.fullHeaderText);
    
    drawnow;
end