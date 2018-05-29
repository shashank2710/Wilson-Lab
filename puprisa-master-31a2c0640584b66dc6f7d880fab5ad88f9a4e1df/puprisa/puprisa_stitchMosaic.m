% stitch a mosiac from tiles with x and y coordinates specified in pixels
% overlapped regions are rewritten with latter tile in the list
% tiles is an n by m by p array,
% where n is n. rows, m is n. cols, and p is n tiles
%
% images are lined up with coordinates rounded to the nearest pixel
% 
% xpos and ypos are lists of x and y coords for each tile,
% in units of pixels
%
% future improvements: subpixel registration
function M = puprisa_stitchMosaic( tiles, xPos, yPos, zPos )
    % since we can't do subpixel registration yet, round coordinates
    xPos = round(xPos);
    yPos = round(yPos);
    figure;
    colormap(gray);
    
    zPosMosaic = unique(zPos);
    nZSlices = length(zPosMosaic);
    
    % get the size and number of tiles
    [nRowPerTile, nColPerTile, nTiles] = size(tiles);
    
    % figure out the total dimensions of the final mosaic
    minX = min(xPos) - nColPerTile/2;
    maxX = max(xPos) + nColPerTile/2;
    minY = min(yPos) - nRowPerTile/2;
    maxY = max(yPos) + nRowPerTile/2;
    
    % make an array large enough to hold the mosaic
    nCol = maxX - minX;
    nRow = maxY - minY;
    M = zeros( nRow, nCol, nZSlices );
    
    % fill the array one tile at a time
    for iTile = 1:nTiles
        % determine which indices to fill with the tile
        i2 = ((xPos(iTile) - nColPerTile/2) ...
                : (xPos(iTile) + nColPerTile/2-1)) - minX + 1;
        i1 = ((yPos(iTile) - nRowPerTile/2) ...
                : (yPos(iTile) + nRowPerTile/2-1)) - minY + 1;
            
        iz = find( zPosMosaic == zPos(iTile) );
            
        % fill that patch
        M(i1, i2, iz) = tiles(:,:,iTile);
        imagesc(squeeze(M(:,:,iz)));
        drawnow;
    end
    close;
end