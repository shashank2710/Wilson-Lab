function rgb = getimagergb(ax)
    % just like getimage, but automatically converts indexed images to RGB,
    % to match the displayed appearance
    
    if nargin == 0
        ax = gca();
    end
    
    [cdata, flag] = getimage(ax);
    
    switch flag
        case {2,3}
            % rescale the data to match the indexed colors
            cm = colormap(ax);
            cl = get(ax,'clim');
            
            [nColors,~] = size(cm);
            scaled = round(scaleImage(cdata, cl, [1, nColors], 1));

            rgb = ind2rgb( uint8(scaled), cm );
        case 4
            % already an RGB image
            rgb = cdata;
        otherwise
            error('getimagergb works with RGB or scaled images only');
    end
    
    % apply alphamap if specified
    him = imhandles(ax);
    A = get(him,'alphadata');
    amapping = get(him,'AlphaDataMapping');

    if length(A) > 1
        if ~strcmp(amapping, 'scaled');
            error('Only scaled alpha data mapping is supported');
        end
        
        % scale the alpha map to range from 0 to 1
        al = get(ax,'ALim');
        A = (A - al(1)) / (al(2)-al(1));
        A(A<0) = 0;
        A(A>1) = 1;
        
        % repeat alpha map to cover R, G, and B channels
        AA = repmat(A,[1,1,3]);
        
        % apply alpha map to the image
        rgb = rgb .* AA;
    end
    
end


function Y = scaleImage(X, oldRange, newRange, truncate)
    % stretch image values from old range to new range
    % for example, to map an image with values 0--1 to a 256-color index,
    % Y = stretchValues(X, [0,1], [1,256]);
    % 
    % or to stretch and truncate a wide range to 0--1:
    % Y = stretchValues(X, [0.1, 502], [0,1], 1);
    
    factor = (newRange(2) - newRange(1)) / (oldRange(2) - oldRange(1));
    
    % offset
    Y = X - oldRange(1);
    
    % now scale
    Y = Y * factor;
    
    % then add new offset
    Y = Y + newRange(1);
    
    if( nargin == 4 )
        if( truncate )
            Y(Y<newRange(1)) = newRange(1);
            Y(Y>newRange(2)) = newRange(2);
        end
    end
end
    