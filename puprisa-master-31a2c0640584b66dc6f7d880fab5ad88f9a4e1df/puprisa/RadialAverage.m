% Warren lab's finest code for taking a fourier transform of a melanin
% contrast image and calculating the radial average with variable bin size
% (change increm to change bin size or switch to not binnin option if no binning
% is preferred.

% Written by Mary Jane with the help of Jesse, Kevin, Prathyush and Michael
% on July 26, 2012.

% must first click on the spectral decomp image
h = gco;
Im = get(h,'CData');
[nr,nc,~] = size(Im);

% melanin needs to be red and green
r = Im(:,:,1);
g = Im(:,:,2);

% apply blackman window to time domain data before taking FFT
wy = repmat(blackman(nr),1,nc);
wx = repmat(blackman(nc).',nr,1);
w = wy.*wx;

% shift the zero frequency to the center
r_f = abs(fftshift(fft2(r.*w)));
g_f = abs(fftshift(fft2(g.*w)));

% make the fourier transform matrix
imFFT = zeros(nr,nc,3);
imFFT(:,:,1) = log(r_f);
imFFT(:,:,2) = log(g_f);

% normalize it
imFFT = imFFT / max(imFFT(:));

% maybe someday might need to use just one type of melanin
%figure('name','fourier transform eumelanin');imshow(imFFT(:,:,1));
%figure('name','fourier transform pheomelanin');imshow(imFFT(:,:,2));

% make the color image gray scale (get rid of eu/pheo information)
imFFTgray = rgb2gray(imFFT);

% show the fourier transform if you want to see it
figure('name','fourier transform melanin');
imshow(imFFTgray);

% make a matrix of the distance of each pixel from the center
d = zeros(nr,nc);
x1 = 1;

% calculate the distance from the center varying x and y values up to the
% maximum nr and nc
% calculate the distance from the center varying the x value
while x1 < nc+1
    y1 = 1;
    % calculate the distance from the center varying the y value
    while y1 < nr+1
        % thank you Pythagoras
        d(x1,y1) = sqrt((x1-nc/2)^2+(y1-nr/2)^2);
        y1 = y1+1;
    end
    x1 = x1+1;
end

% make the distance matrix a list
d = reshape(d,nr*nc,1);

% make the FFT value matrix a list
f = reshape(imFFTgray,nr*nc,1);

% combine lists into a single matrix
dfmatrix = zeros(nr*nc,2);
dfmatrix(:,1) = d;
dfmatrix(:,2) = f;

% sort data in ascending order of distance value
dfmatsorted = sortrows(dfmatrix,1);

% average the FFT values which have the same distance value
dfAve = zeros(10, 2);

% start with averaged matrix (matrix being written) row 1
RowAveMat = 1;

% start with sorted matrix row 1
RowSortMat = 1;

% have bins start with zero
bin1 = 0;

% increment by any value desired - may add this as a prompt
increm = .5;

% make a temporary matrix to hold values with the same distance values
tempMat = zeros(1,2);
while RowSortMat < nr*nc;
    % make the first column values of the temporary matrix all the same 
    % distance values
    tempMat(:,1) = dfmatsorted(RowSortMat,1);
    
    % make the second column values be all of the FFT values that 
    % correspond to that single distance value
    tempRow = 1;
    tempMat(tempRow,2) = dfmatsorted(RowSortMat,2);
    
    % if binning
    % see if the distance values are between values of bin1 and bin1+increm
        while dfmatsorted(RowSortMat,1)>=bin1 && dfmatsorted(RowSortMat,1)<bin1+increm && RowSortMat < nr*nc
            
            % if it is, add the FFT value to the temporary matrix
            tempMat(tempRow,2) = dfmatsorted(RowSortMat,2);
            
            % then go on to the next value to see if its below i+3
            RowSortMat=RowSortMat+1;
            
            % and go on to the next row in the temporary matrix
            tempRow=tempRow+1;
           
            % end when it finds a distance value that is more than the
            % increment value
        end
        
    % if not binning
        %while dfmatsorted(RowSortMat,1) == dfmatsorted(RowSortMat+1,1)
            % if they are, add the FFT value to the temporary matrix
            %tempMat(tempRow+1,2) = dfmatsorted(RowSortMat+1,2);
            % then go on to the next value to see if its the same
            %RowSortMat=RowSortMat+1;
            % and go on to the next row in the temporary matrix
            %tempRow=tempRow+1;
        %end when it finds a distance value that is different from the 
        %previous
        %end
    
    % increment bins
    bin1 = bin1+increm;
    
    % average the temporary matrix values
    dAve = tempMat(1,1);
    fAve = mean(tempMat(:,2));
    
    % put them into the averaged distance and fft value matrix one at a 
    % time
    dfAve(RowAveMat,1) = dAve;
    dfAve(RowAveMat,2) = fAve;
        
    % increment through the sorted df value matrix starting with the
    % the next row in the sorted matrix
    RowSortMat = RowSortMat+1;
        
    % prepare to add the next averaged values to the next row of the 
    % averaged matrix
    RowAveMat = RowAveMat+1;
        
    % clear temporary matrix and remake
    clear tempMat;
    tempMat = zeros(1,2);
end

% see fabulous plot
figure('name','radially averaged fourier transform'); 
plot(dfAve(:,1),dfAve(:,2));
% width as soon as it starts dropping times 2.4048 equals the FW radius of the
% disk in pixel units