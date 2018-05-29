function out = meanstd(S, header)

    [~,X] = puprisa_getChannelFromSlices(S, 1, header);
    
    out.mean = mean(X(:));
    out.std = std(X(:));
    out.sum = sum(X(:));