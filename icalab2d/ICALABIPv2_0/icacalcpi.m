function [G, pi] = icacalcpi( W, Q, A );

% ICACALCPI.M
% Copyright (c) by Krzysztof Siwek

G = abs( W * Q * A ) ;
[Nrows, Ncols] = size( G );
N = Nrows;
[I, J] = max( G' );        
G = pinv(diag(I))*G;
G = G';
if Nrows ~= Ncols,
   pi = -1;
else
   pi = (sum(sum(G)-max(G))/(Nrows*Ncols-length(max(G))));
end
G = G';

