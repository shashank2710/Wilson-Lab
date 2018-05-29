% generate common principal components basis for heart and BAT tissue

% files that had 8 delay points
files = {'bat_01\062416c_CH2.tif', ...
    'bat_01\062416d_CH2.tif', ...
    'heart_01\070116b_CH2.tif', ...
    'heart_01\070116c_CH2.tif'};
     
% load each file

do_denoise = 0;
binfactor = 2;
for ii = 1:length(files)
    dat = puprisa_readImageStack(files{ii});
    [t,X,~] = puprisa_getChannelFromSlices(dat, 2, '');
    
    delays{ii} = t;
    
    if do_denoise
        Xmin = min(X(:));
        Xmax = max(X(:));
        Xs = single( (X - Xmin)/(Xmax - Xmin) );

        % de-noise
        Xdn = NLMF(Xs);

        % re-scale
        Xdn = single((Xdn * (Xmax - Xmin)) + Xmin);
        X = Xdn;
    end
    
    X = imresize(X,1/binfactor);
    delayStacks{ii} = X;
end

%% plot averaged curves
clf;
for ii = 1:length(files)
    X = delayStacks{ii};
    [ny, nx, nt] = size(X);
    ds = reshape(X, [ny*nx, nt] );
    ds = ds - mean(ds(:,1)); % subtract mean offset
    x = sum(ds,1);
    x = x/max(x);
    
    plot(t,x);
    hold on;
    
end
legend('bat 1', 'bat 2', 'heart 1', 'heart 2');
xlabel('pump-probe delay, ps');
ylabel('peak-normalized average signal');

%% make a list of delay scans
DS = [];
isheart = [];
isbat = [];
for ii = 1:length(files)
    X = delayStacks{ii};
    
    [ny, nx, nt] = size(X);
    
    ds = reshape(X, [ny*nx, nt] );
    ds = ds - mean(ds(:,1)); % subtract mean offset
    DS = [DS; ds];
    
    % keep track of signal origin
    if( strfind(files{ii},'heart') )
        isheart = [isheart, ones(1,ny*nx)];
    else
        isheart = [isheart, zeros(1,ny*nx)];
    end
       
    if( strfind(files{ii},'bat') )
        isbat = [isbat, ones(1,ny*nx)];
    else
        isbat = [isbat, zeros(1,ny*nx)];
    end
end

%% show the sum-squared signal present in each delay scan
p = sum(DS.^2,2);
plot(p);
xlabel('pixel index')
ylabel('sum-squared signal')

%% histogram of power
hist(p,2048);
xlabel('sum-squared delay scan magnitue');
ylabel('pixel count');

%% independent thresholding of each image, otsu's method
% (does not work well)
clf;
thr = 1;
colormap(gray.^0.5)
for ii = 1:length(files)
    X = delayStacks{ii};
    [ny, nx, nt] = size(X);
    ds = reshape(X, [ny*nx, nt] );
    ds = ds - mean(ds(:,1)); % subtract mean offset
    p = sum(ds.^2,2);
    
    I_p = reshape(p, [ny,nx]);
    %I_thr = reshape(p > thr,[ny,nx]);
    subplot(4,2,ii);
    imagesc(I_p);
    
    thr = graythresh(I_p);
    I_thr = reshape(p > thr,[ny,nx]);
    subplot(4,2,ii+4);
    imagesc(I_thr);
    
end



%% test thresholding
clf;
thr = 0.3;
for ii = 1:length(files)
    X = delayStacks{ii};
    [ny, nx, nt] = size(X);
    
    
    ds = reshape(X, [ny*nx, nt] );
    ds = ds - mean(ds(:,1)); % subtract mean offset
    p = sum(ds.^2,2);
    
    I_thr = reshape(p > thr,[ny,nx]);
    subplot(2,2,ii);
    imagesc(I_thr);
end
annotation('textbox',[0.4,0.9,0.3,0.1],'string',...
    sprintf('threshold %0.2f',thr),'edgecolor','none');




%% threshold
clf;
p = sum(DS.^2,2);
sum((p > thr)) / length(p)
DS_thr = DS(p>thr,:);
isheart_thr = isheart(p>thr);
isbat_thr = isbat(p>thr);

%% run PCA
[U,S,V] = svd(DS_thr.',0);
plot(t,U(:,1:4));
xlabel('probe delay, ps');
ylabel('signal, arb.');
legend('PC1','PC2','PC3','PC4');
title(sprintf('PCA, thresholded at > %f',thr));

%% project onto top 3 PCs and subtract estimated offset
DS_prj = (DS_thr*U(:,1:3))*U(:,1:3).';
DS_thr_offs = DS_thr - repmat(DS_prj(:,1),[1,nt]);
plot(1:size(DS_thr,1),DS_thr(:,1),...
        1:size(DS_thr_offs,1),DS_thr_offs(:,1))
legend('original','offset-corrected');
grid on


%% run PCA on offset-corrected data
[U,S,V] = svd(DS_thr_offs.',0);
plot(t,U(:,1:3));
xlabel('probe delay, ps');
ylabel('signal, arb.');
legend('PC1','PC2','PC3');
title('PCA on thresholded offset-corrected data');

%% project data onto new PCA basis
DS_prj2_coefs = (DS_thr_offs*U(:,1:3));
theta = atan2(DS_prj2(:,1), DS_prj2(:,2));
hist(theta,128)
ylabel('pixel count');
xlabel('atan(PC1 / PC2)')

%% scatter plot
plot(DS_prj2_coefs(:,1),DS_prj2_coefs(:,2),'.')
xlabel('PC1');ylabel('PC2');

%% scatter plot heart vs BAT
clf
plot(DS_prj2_coefs(isheart_thr==1,1),DS_prj2_coefs(isheart_thr==1,2),'o')
hold on
plot(DS_prj2_coefs(isbat_thr==1,1),DS_prj2_coefs(isbat_thr==1,2),'s')
hold off
legend('heart','brown adipose');
xlabel('PC1');ylabel('PC2');
grid on


%% try LDA

% concatenate data for classifier
dat = DS_thr_offs;
dat(:,end+1) = isheart_thr;

%% concatenate for classifier, leave offset

dat = DS_thr;
dat(:,end) = isheart_thr;

%% try applying classifier to images
clf;

for ii = 1:length(files)
    X = delayStacks{ii};
    [ny, nx, nt] = size(X);
    ds = reshape(X, [ny*nx, nt] );
    ds = ds - mean(ds(:,1)); % subtract mean offset
    p = sum(ds.^2,2);
    
    I_p = reshape(p, [ny,nx]);
    %I_thr = reshape(p > thr,[ny,nx]);
    subplot(2,2,ii);
    imagesc(I_p);
    
    yfit = trainedClassifier.predictFcn(ds);
    I_class = reshape(yfit,[ny,nx]);
    
    I_hsv = [];
    I_hsv(:,:,1) = I_class*0.3+0.1;
    I_hsv(:,:,2) = 1;
    I_hsv(:,:,3) = I_p.^0.7;
    imshow(hsv2rgb(I_hsv));
    title(files{ii},'interp','none');
end
