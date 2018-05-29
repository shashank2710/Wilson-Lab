% principal component analysis heterogeneity map
function puprisa_PCA_hetMap( X )

    
   % reshape to a list of delay scans
   [nrows,ncols,ndelays] = size(X);
   delayScans = zeros(nrows*ncols, ndelays);
   for irow = 1:nrows
       minIndex = (irow-1)*ncols+1;
       maxIndex = (irow)*ncols;
       delayScans( minIndex:maxIndex, : ) = squeeze(X( irow, 1:ncols, : ));
   end
   
   O = zeros(nrows,ncols);
   
   hwb = waitbar(0,'calculating overlap image...')
   
    for irow=2:(nrows-1)
        for icol=2:(ncols-1)
            % select 3x3 pixels in region
           
            Y = X( (irow-1):(irow+1), (icol-1):(icol+1), : );
           
            yy = [squeeze(Y(1,1,:)), squeeze(Y(1,2,:)), squeeze(Y(1,3,:)), ...
                  squeeze(Y(2,1,:)), squeeze(Y(2,2,:)), squeeze(Y(2,3,:)), ...
                  squeeze(Y(3,1,:)), squeeze(Y(3,2,:)), squeeze(Y(3,3,:))].';
             
            %y0 = squeeze(Y(2,2,:)).';
            
            % normalize
            %n = sqrt(sum(yy.*yy,2));
            %yy = yy ./ repmat(n,[1,ndelays]);
            
            %n0 = sqrt(sum(y0.*y0,2));
            %y0 = y0 / n0;
       
            % permute overlap dot product
            %O(irow,icol) = sum(sum(triu(yy*yy.',1)));
            
            % dot product with center pixel
            %O(irow,icol) = sum(y0*yy.') / 8;
            
            % subtraction similarity
            %q = sum( repmat(y0,[8,1])- yy, 2 );
            %O(irow, icol) = sqrt( sum(q.^2) );
            
            %PCA similarity measure
            % subtract means
            %chanMeans = mean( yy );
            %yy = yy - repmat( chanMeans, 3*3, 1 );

            [~,S,~] = svd(single(yy),'econ');
            s = diag(S);
            O(irow,icol) = s(2) / s(1);
           
        end
        
        waitbar( irow / nrows, hwb );
    end
    
    close(hwb);
    
   figure();
   img = imagesc(O);
   set(gca,'ydir','normal','color','k');
   colormap(jet);caxis([0,1]);
   
   Xbg = sqrt(sum((X.^2),3));
   set(img,'alphadata',Xbg,'alphadatamapping','scaled');
   alim([0, mean(Xbg(:))+3*std(Xbg(:))])
   axis image;
   
    colorbar;
end