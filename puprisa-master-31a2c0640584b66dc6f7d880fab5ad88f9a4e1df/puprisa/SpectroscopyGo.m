    % Exports spectroscopy curves from regions of interest
    tlc = findobj(gcf,'tag','lTransient');
    x1 = get(tlc(1),'xdata');
    y1 = get(tlc(1),'ydata');
    %x2 = get(tlc(2),'xdata');
    %y2 = get(tlc(2),'ydata');
    %Normalizes minimum intensity to 1
    %if abs(max(y1))>= abs(min(y1))
    %    y1 = y1/max(y1);
    %else 
    %    y1 = -y1/min(y1);    
    %end
    %if abs(max(y1))>= abs(min(y1))
    %    y2 = y2/max(y2);
    %else
    %    y2 = -y2/min(y2);
    %end
    max1 = max(y1)
    min1 = min(y1)
    %max2 = max(y2)
    %min2 = min(y2)
    %ratio = max1/(-min1)