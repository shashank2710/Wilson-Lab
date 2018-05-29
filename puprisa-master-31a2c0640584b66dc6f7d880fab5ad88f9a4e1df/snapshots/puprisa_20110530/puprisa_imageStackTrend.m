function trend = puprisa_imageStackTrend( imageStack )
    [nRow, nCol, ~] = size(imageStack);
    nPix = nRow * nCol;
    trend = squeeze(sum(sum(abs(imageStack),1),2)) / nPix;
end