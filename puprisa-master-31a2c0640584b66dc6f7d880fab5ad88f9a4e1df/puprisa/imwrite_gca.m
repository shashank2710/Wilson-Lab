function imwrite_gca( varargin )
    % similar to imwrite
    % but grabs an indexed image and colormap from current axes
    % need to modify this to write out the entire image, with graphics
    % objects overlaid... or see if we can make something to export a
    % figure
    cm = colormap();
    cl = caxis();
    im = findobj(gca,'Type','Image');
    
    if isempty(im)
        error('Must have an image-containing axes as current axes.');
    end
    
    
    cdata = get(im,'cdata');
    
    [nx,ny,nc] = size(cdata);
    
    if nc == 1
        % rescale the data to match the indexed colors
        [nColors,~] = size(cm);
        scaled = round(scaleImage(cdata, cl, [1, nColors], 1));


        rgb = ind2rgb( scaled, cm );
    elseif nc == 3
        rgb = cdata;
    end
    
    imwrite( rgb, varargin{:} );
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
    