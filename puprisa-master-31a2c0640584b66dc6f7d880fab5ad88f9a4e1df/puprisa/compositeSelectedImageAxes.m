function im = compositeSelectedImageAxes(mode,varargin)
% make RGB composite of selected image axes
%
% USAGE:
% To additively combine two axes, select them, and call
% compositeSelectedImageAxes() with no arguments.
%
% To combine with a threshold, select 2 axes, and then call, e.g.
% compositeSelectedImageAxes( 'thresholdedReplace', 10 );
%
% to reverse the order,
% compositeSelectedImageAxes( 'thresholdedReplace', 10, 'reverse' );


    if nargin == 0
        mode = 'sum';
    end

    h = findobj(0,'selected','on','type','image');

    if nargin == 3
        if( strcmp(varargin{2}, 'reverse') )
            h=flipud(h);
        end
    end
    
    im = [];
        figure();
    
    for ih = 1:length( h );
        ax = get(h(ih),'parent');
        thisImageRGB = getimagergb( ax );
        
        if isempty(im)
            im = thisImageRGB;
        else
            switch mode
                case 'sum'
                    im = im + thisImageRGB;
                case 'lightenOnly'
                    mask = sum(thisImageRGB,3) > sum(im,3);
                    im = im + thisImageRGB.*repmat(mask,[1,1,3]);
                case 'thresholdedReplace'
                    mask = sum(thisImageRGB,3) > varargin{1};
                    mask = repmat(mask,[1,1,3]);
                    im = im.*(1-mask) + thisImageRGB.*mask;
            end
        end
        
        imshow( im );
        pause(1)
        
    end
    
    if nargout == 0
        figure();
        imshow( im );
        set(gca,'ydir','normal')
        im = [];
    end
end
        