function puprisa_clusterAnalysis( delayScans )
    % pick the strongest n signal points
    %nStrongestSignals = 2000;
    %strengths = sum(abs(delayScans),2);
    %[b,ii] = sort(strengths,'Descend');
    
    
    [U,S,V] = svd( delayScans, 0 );
    size(U)
    size(S)
    size(V)
    
    %s = diag(S);
    
    % normalize each delay scan by area under the curve
    %delayScans = delayScans( ii(1:nStrongestSignals), : );
    %[nscans, ~] = size(delayScans)
    %delayScanSums = repmat(sum(abs(delayScans),2),[1,ndelays]);
    %delayScans = delayScans ./ delayScanSums;
    
    projection = delayScans*V;
    
    norms = sqrt(sum(projection(:,(1:3)).^2,2));
    threshold = 0.1;
    %projection = projection( norms > threshold,: );
    
    % scatter plot of projections with top 3 components
    x = projection(:,1);
    y = projection(:,2);
    z = projection(:,3);
    
    % figure out k means clustering
    % what we really want to do is k means weighted by distance from center
    k = 2;
    [idx,c] = kmeans(projection(:,(1:3)),k,'distance','cosine');
    figure('Name','PC Clustering');
    cm=colormap(lines);
    for ik = 1:k
        clustNorms = sqrt(sum([x(idx==ik),y(idx==ik),z(idx==ik)].^2,2));
       maxClustNorm = max(clustNorms);
        
        % scatter plot
        line(x(idx==ik),y(idx==ik),z(idx==ik),'LineStyle','none',...
            'Marker','.','MarkerSize',1.0,'color',cm(ik,:));
        
        % draw cluster vector
        cvect = [c(ik,1),c(ik,2),c(ik,3)];
        cvect = cvect / sqrt(sum(cvect.^2));
        cvect = cvect * maxClustNorm;
        line([0,cvect(1)],[0,cvect(2)],[0,cvect(3)],'color',cm(ik,:));
    end
    xlabel('PC1');
    ylabel('PC2');
    zlabel('PC3');
    return;
    % broken after cutting by threshold
    figure();
    rr = x;
    gg = y;
    bb = z;
    for ik = 1:k
        rr(idx == ik) = cm(ik,1);
        gg(idx == ik) = cm(ik,2);
        bb(idx == ik) = cm(ik,3);
    end
    rr = reshape(rr,[nrows,ncols]);
    gg = reshape(gg,[nrows,ncols]);
    bb = reshape(bb,[nrows,ncols]);
    imRGB = zeros(nrows,ncols,3);
    imRGB(:,:,1) = rr;
    imRGB(:,:,2) = gg;
    imRGB(:,:,3) = bb;
    
    imshow(imRGB);
end