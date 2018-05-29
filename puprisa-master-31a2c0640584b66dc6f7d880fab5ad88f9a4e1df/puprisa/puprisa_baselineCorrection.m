function [X, baseline] = puprisa_baselineCorrection( t, X, rank )
    % SVD-based delay stack offset correction
    %
    % Determines independent baseline correction for each pixel in a delay
    % stack, using a reduced-rank approximation to de-noise the delay
    % scans.
    %
    % THIS ASSUMES:
    %  negative time delays were sampled (at times when the probe precedes
    %  the pump)
    %
    % typical rank 6 should work for multiexponentials
    %
    % This function works with delay stacks and lists of delay scans
    
    if ndims( X ) == 3
        % the data supplied is a delay stack
    
        [nr,nc,nt] = size(X);
    
        % reshape 3D delay stack into 2D array; each pixel is a row
        D = reshape(X,[nr*nc,nt]);
    elseif ndims( X ) == 2
        % the data supplied is already a list of delay scans
        [~, nt] = size( X );
        D = X;
    else
        error('Data X must be a 2- or 3-dimensional array');
    end

    % singular value decomp
    [U,S,V] = svd(D,0);

    % reduce rank
    s = diag(S);
    s((rank+1):end) = 0;
    S = diag(s);

    % reconstruct delay scans from reduced-rank representation
    D2 = U*S*V';

    [~,it0] = min(t);   % find index of earliest time delay
    baseline = D2(:,it0);
    
    D = D - repmat(baseline,1,nt);
    
    if ndims( X ) == 3
        % re-shape back into a 3D delay stack
        X = reshape(D,[nr,nc,nt]);
    else
        % simply return the list of delay scans
        X = D;
    end