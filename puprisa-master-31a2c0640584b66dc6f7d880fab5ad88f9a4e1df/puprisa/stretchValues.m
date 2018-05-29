function Y = stretchValues(X, oldRange, newRange, truncate)
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
    