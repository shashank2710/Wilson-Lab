function hg = puprisa_scaleBar( length_in_microns, scale_px_per_micron )
    % todo:
    % set up so that zoom events cause the scalebar to redraw appropriately
     
    oldAxes = gca();
    
    % delete old scalebar, if present
    
    % set up patch object for scalebar
    width_px = length_in_microns * scale_px_per_micron;
    
    height = 5;
    
    x = [0,      0, width_px, width_px, 0];
    y = [0, height,   height,        0, 0];
    
    x0 = 20;
    y0 = 20;
    
    x = x + x0;
    y = y + y0;
    
    hpatch = patch(x,y,'white','tag','scaleBarPatch');
    
    
    % set up text label for scalebar
    ymargin = 2;
    tx0 = (2*x0 +  width_px)/2;
    ty0 = y0 + height + ymargin;
    
    htxt = text(tx0, ty0, [num2str(length_in_microns),' \mum'],...
        'color','white','EdgeColor','none',...
        'HorizontalAlignment','center',...
        'VerticalAlignment','bottom',...
        'FontWeight','bold','Tag','scaleBarText');
    
    % gather scalebar objects into a group
    hg = hggroup;
    set(hg,'tag','scaleBarGroup');
    setappdata(hg,'scale_px_per_micron',scale_px_per_micron);
    setappdata(hg,'length_in_microns',length_in_microns);
    set(hpatch,'parent',hg);
    set(htxt,'parent',hg);
    
    % set mouseDown events so that double-clicking the scalebar brings
    % up a dialog to change
    set([hpatch,htxt],'ButtonDownFcn',@scaleBarButtonDown);
    
    % switch back to old axes
    
    axes(oldAxes);

end

function scaleBarButtonDown(src, evt)
    % When the user double-clicks the scalebar, raise a dialog box to
    % prompt for the desired scale in microns

    % examine the SelectionType to determine whether the user
    % double-clicked on the scalebar
    if( strcmp( 'open', get(gcbf,'SelectionType') ) )
        % first, look up the old scale to use as a default response in the
        % prompt.
        
        % find parent axes
        axParent = get(src,'parent');

        % find scalebar group
        hg = findobj( axParent,'tag','scaleBarGroup');
        oldScale = getappdata(hg,'length_in_microns');
            
        % prompt user for new scale
        a = inputdlg('Enter new scale (microns)','Scale',1,{num2str(oldScale)});
        
        % if they hit 'cancel', then quit
        if isempty(a)
            return;
        end

        % convert the entry to a number
        l = str2num(a{1});
        
        % if not a valid number, then quit
        if isempty(l)
            return;
        end

        changeScale(src,l);
    end
end

function changeScale( h, length_in_microns )
    % change a scalebar to a new length
    % h is a handle to any GUI element in the scalebar
    
    % find parent axes
    axParent = get(h,'parent');
    
    % find scalebar group
    hg = findobj( axParent,'tag','scaleBarGroup');
    
    % retrieve the scale in pixels per micron
    scale_px_per_micron = getappdata(hg,'scale_px_per_micron');

    % set the patch object's width accordingly
    
    width_px = length_in_microns * scale_px_per_micron
    
    hpatch = findobj(hg,'tag','scaleBarPatch');
    x = get(hpatch,'xdata');
    x(3) = x(1) + width_px;
    x(4) = x(1) + width_px;
    set(hpatch,'xdata',x);
    
    % adjust the text label
    htxt = findobj(hg,'tag','scaleBarText');
    ymargin = 2;
    y0 = 20;    
    height = 5;
    tx0 = (2*x(1) +  width_px)/2;
    ty0 = y0 + height + ymargin;
    set(htxt,'Position',[tx0,ty0]);
    set(htxt,'String', [num2str(length_in_microns),' \mum']);
    
    % remember the new scale 
    setappdata(hg,'length_in_microns',length_in_microns);

end