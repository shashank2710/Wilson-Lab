function [t,X, z] = puprisa_getChannelFromSlices(S, chan, header)
    % parses puprisa imageStack structure to get a single channel
    %
    % Example: load an imagestack, and extract channel 1, then display the
    % first frame in that stack
    %
    %   S = puprisa_readImageStack( fileName );
    %   [t,X] = getChannelFromSlices(S,1);
    %   imagesc(X(:,:,1));

    t = [S.delays];
    z = [S.posZ];
    nt = length(t);

    [nr,nc] = size(S(1).imageData{chan});

    X = zeros(nr,nc,nt);

    for it = 1:nt
        X(:,:,it) = S(it).imageData{chan};
    end
    
    % is it a mosaic?
    isMosaic = 0;
    if isfield( header, 'mosaic' )
        if header.mosaic == 1
            isMosaic = 1;
        end
    end
    
    if isMosaic
        % if so, then stitch it together first
        
        % convert positions to pixels
        if ~ispref('puprisa','scaleX_volts_per_micron')
            error('Mosaics require field of view calibration');
        end

        scaleX_volts_per_micron = ...
            getpref('puprisa','scaleX_volts_per_micron');
        scaleY_volts_per_micron = ...
            getpref('puprisa','scaleY_volts_per_micron');

        [nx,ny,nt] = size(X);

        scanRangeX_volts = header.scanrangex;
        scaleX_px_per_micron = scaleX_volts_per_micron * nx / scanRangeX_volts;
    
        scanRangeY_volts = header.scanrangey;
        scaleY_px_per_micron = scaleY_volts_per_micron * ny / scanRangeY_volts;
        
        setappdata(gcbf,'scaleX_px_per_micron',scaleX_px_per_micron);
        setappdata(gcbf,'scaleY_px_per_micron',scaleY_px_per_micron);
    
        
        xPos = cell2mat({S.posX});
        yPos = cell2mat({S.posY});
        t = cell2mat({S.delays});
        
        xPos_px = xPos * scaleX_px_per_micron;
        yPos_px = yPos * scaleY_px_per_micron;

        X = flipdim(X,1);
        X = puprisa_stitchMosaic(X,xPos_px,yPos_px, t);
        
        t = unique(t);  % reduce delay vector to those used for a single tile
        z = unique(z);
    end
end