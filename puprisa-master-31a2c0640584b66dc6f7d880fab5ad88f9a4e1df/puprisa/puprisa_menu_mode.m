function puprisa_menu_mode( )
    % create a menu to change mode from pumpProbe to fluor
    % input args: figure, and relevant axes object(s)
    
    f = uimenu('Label','ImageMode');
    
    uimenu(f,'Label','Pump-probe','Callback',@mnuPumpProbe,...
        'Checked','on','Tag','mnuModePumpProbe' );
    uimenu(f,'Label','Fluorescence','Callback',@mnuFluor,...
        'Checked','off','Tag','mnuModeFluor' );
    uimenu(f,'Label','FluorescenceInv','Callback',@mnuFluorInv,...
        'Checked','off','Tag','mnuModeFluorInv' );
    uimenu(f,'Label','Confocal','Callback',@mnuConfocal,...
        'Checked','off','Tag','mnuModeConfocal' );
    
    uimenu(f,'Label','Stretch Contrast','Separator','on',...
        'Callback',@mnuStretchContrast );
end

function mnuPumpProbe(src, evt)
    h = guihandles(gcbf);
    
    % change colormap
    colormap(puprisa_colorMap('pumpProbe'));
    
    % change to symmetric scaling
    cl = get(h.imageAxes,'clim');
    mxCl = max(abs(cl));
    set(h.imageAxes,'clim', [-1,1]*mxCl );
    
    % change checkboxes
    checkOnly(h.mnuModePumpProbe);
    
    setappdata(gcbf,'imagePolarity','bipolar');
end

function mnuFluor(src,evt)
    h = guihandles(gcbf);
    
    % change colormap
    colormap(puprisa_colorMap('fluor'));
    
     % change to non-symmetric scaling
    cl = get(h.imageAxes,'clim');
    set(h.imageAxes,'clim', [cl(1),0] );
    
    % change checkboxes
    checkOnly(h.mnuModeFluor);
    
    setappdata(gcbf,'imagePolarity','negative');
end

function mnuFluorInv(src,evt)
    h = guihandles(gcbf);
    
    % change colormap
    colormap(flipud(puprisa_colorMap('fluor')));
    
     % change to non-symmetric scaling
    cl = get(h.imageAxes,'clim');
    set(h.imageAxes,'clim', [0,cl(2)] );
    
    % change checkboxes
    checkOnly(h.mnuModeFluorInv);
    
    setappdata(gcbf,'imagePolarity','positive');
end

function mnuConfocal(src,evt)
    h = guihandles(gcbf);
    
    % change colormap
    colormap(gray);
    
     % change to non-symmetric scaling
    cl = get(h.imageAxes,'clim');
    set(h.imageAxes,'clim', [0,cl(2)] );
    
    % change checkboxes
    checkOnly(h.mnuModeConfocal);
    
    setappdata(gcbf,'imagePolarity','positive');
end

function checkOnly( h )
    parent = get(h,'parent');
    allMenuItems = get(parent,'children');
    set(allMenuItems,'checked','off');
    set(h,'checked','on');
end

function mnuStretchContrast( src, evt )
    % enable/disable contrast stretching
    
    if strcmp( get( src, 'Checked' ), 'on' )
        % item was previously checked--disable this feature
        setappdata( gcbf, 'stretchContrast', 0 );
        set( src, 'Checked', 'off' );
    else
        % item was previously unchecked--enable this feature
        setappdata( gcbf, 'stretchContrast', 1 );
        set( src, 'Checked', 'on' );
    end
end