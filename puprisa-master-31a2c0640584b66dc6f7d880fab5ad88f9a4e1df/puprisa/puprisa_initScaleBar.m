function puprisa_initScaleBar( ax, header )
% calculates scale in pixels per micron,
% then draws the scalebar
    
    if ispref('puprisa','scaleX_volts_per_micron')
        scaleX_volts_per_micron = ...
            getpref('puprisa','scaleX_volts_per_micron');
        
        % make sure the scale actually contains a number before proceeding
        if ~isempty( scaleX_volts_per_micron )
            
            % calculate scale in terms of pixels per micron
            nx = header.pixelsperline;
            scanRangeX = header.scanrangex;

            scaleX_px_per_micron = scaleX_volts_per_micron * nx / scanRangeX

            setappdata(gcbf,'scaleX_px_per_micron',scaleX_px_per_micron);

            % draw the scalebar
            puprisa_scaleBar(50,scaleX_px_per_micron);
        end
    end
end