function out = eupheo(S, header)
    [t,X] = puprisa_getChannelFromSlices(S, 1, header);
    [F,names]= puprisa_linearSpectralDecomp(X,t,'batch');
    
    out = [];
    out = setfield(out, names{1}, F(1));
    out = setfield(out, names{2}, F(2));
    
    out.euOverPheo = F(2) / F(1);
    
    out