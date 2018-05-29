function cm = puprisa_colorMap(cmName)
% returns a colormap for different purposes
    switch cmName
        case 'fluor'
            % looks great on screen, but terrible in print
            cm = cmap(256,[1,1,1],[0,1,0],[0,0,0]);
        case 'fluorPrint'
            % looks great on a laser printer
            cm = cmap(256,[1,1,1],[.8,1,.8],[.5,1,.5],[0,1,0],[0,0.1,0]);
        case 'fluorHDR1'
            cm = cmap(256,[1,1,1],[0,1,0],[1,0,1],[0,0,0]);
        case 'fluorHDR2'
            cm = cmap(256,[1,1,1],[0,1,1],[0,1,0],[0,0,0]);
        case 'fluorHDR3'
            cm = cmap(256,[1,1,1],[1,1,0],[1,0,0],[0,0,0]);
        case 'fluorHDR4'
            cm = cmap(256,[1,1,1],[1,0,1],[0,1,1],[0,0,0]);
        case 'fluorHDR5'
            cm = cmap(256,[1,1,1],[1,0,1],[0,1,1],[0,0,0]);
        case 'fluorHDR6'
            % my favorite so far
            cm = cmap(256,[1,1,1],[0,1,0],[1,1,0],[0,0,0]);
        case 'fluorHDR7'
            cm = cmap(256,[1,1,1],[0,1,1],[0,1,0],[1,1,0],[0,0,0]);
        case 'fluorHDR8'
            cm = cmap(256,[1,1,1],[1,0,1],[0,1,0],[1,1,0],[0,0,0]);
        case 'fluorHDR9'
            cm = cmap(256,[1,1,1],[1,0,1],[0,1,1],[0,1,0],[1,1,0],[0,0,0]);
        case 'fluorMasters'
            cm = cmap(256,[1,1,1],[1,1,0],[0,1,0],[1,0,1],[0,0,0]);
        case 'fluorSimple',
            cm = cmap(256, [0 1 0],[0 0 0]);
        case 'SHG'
            cm = cmap(256,[0.6,0.5,1],[0,0,1],[0,0.2,0.3],[0,0,0]);
        case 'pumpProbe'
            cm = cmap(256,[0,1,1],[0,0,1],[0,0,0],[1,0,0],[1,1,0]);
        case 'pumpProbeSimple'
            cm = cmap(256,[0 0 1],[0 0 0],[1 0 0]);
        case 'euPheo'
            cm = cmap(256,[0 1 0],[1 1 0],[1 0 0]);
    end
    
if strfind(cmName,'fluor')
    brightScale = linspace(1,0,256);
    for ii = 1:256
        c = cm(ii,:);
        n = sqrt(sum(c.^2));
        c = c / n * sqrt(3);
        c = c*brightScale(ii);

        if( n == 0 )
            cm(ii,:) = [0 0 0];
        else

            cm(ii,:) = c;
        end
    end

    cm( cm > 1 ) = 1;
end

    if nargout == 0
        colormap(cm);
        cm = [];
    end

end