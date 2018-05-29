% puprisa_appendNeighborhood
% append a 3x3 neighborhood to each pixel in the stack
% This is the first step in context-sensitive PCA or other analysis
%
% See Hastie, Tibshirani, Friedman, pg 470-471, 2nd ed

function Xn = puprisa_appendNeighborhood( X )
    [ny0,nx0,nt] = size(X);
    
    % discard edge pixels that have no neighborhood
    ny = ny0 - 2;    
    nx = nx0 - 2;
    
    Xn = single(zeros(ny,nx,nt*9));
    
    for it = 1:nt
        % append neighbors to this image
        Xci = im2col(X(:,:,it),[3 3],'sliding');
        
        % Each column of Xci now is 3x3 neighborhood for each pixel, at
        % time delay it.
        % Each row of Xci is a shifted image, centered at the neighbor
        
        % assign row column of Xci to the appropriate slice of Xn
        Xn(:,:,it       ) = reshape(Xci(1,:),[ny,nx]);
        Xn(:,:,it +   nt) = reshape(Xci(2,:),[ny,nx]);
        Xn(:,:,it + 2*nt) = reshape(Xci(3,:),[ny,nx]);
        Xn(:,:,it + 3*nt) = reshape(Xci(4,:),[ny,nx]);
        Xn(:,:,it + 4*nt) = reshape(Xci(5,:),[ny,nx]);
        Xn(:,:,it + 5*nt) = reshape(Xci(6,:),[ny,nx]);
        Xn(:,:,it + 6*nt) = reshape(Xci(7,:),[ny,nx]);
        Xn(:,:,it + 7*nt) = reshape(Xci(8,:),[ny,nx]);
        Xn(:,:,it + 8*nt) = reshape(Xci(9,:),[ny,nx]);
        
        
    end
end