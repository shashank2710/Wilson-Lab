% puprModelFit.m
%
% Fits pump-probe data to a model consisting of a linear combination of an
% instantaneous response, and an arbitrary number of bipolar
% multiexponentials, all convolved with an Gaussian instrument response
% function. This approach is based on the inverse Laplace method outlined
% in Y. Song, “Resolution and uncertainty of Laplace inversion spectrum,” 
% Magnetic Resonance Imaging, vol. 25, no. 4, pp. 445–448, May 2007.

%
% INPUTS:
%   Ymeas:    Measured pump-probe data. For a single measurement, this is
%             just a row vector. For multiple measurements (e.g. all pixels
%             of a delay stack), each measurement is on a separate row.
%
%   T:        Measurement timebase (the probe delays sampled during the
%             measurement that correspond to the columns of Ymeas). Units:
%             picoseconds.
%
%   IRFwidth: Instrument response function width, full width half maximum.
%             This is usually the cross-correlation, as measured by 
%             2-photon absorption in rhodamine 6G. Units: picoseconds.
%
%   minDecayTime: The minimum 1/e decay to consider in the model (typ 0.1
%                 ps) Units: picoseconds.
%
%   maxDecayTime: The maximum 1/e decay to consider in the model. This 
%                 should be greater than the maximum expected decay time. 
%                 (typically 30 ps for melanins) Units: picoseconds.
%
% OUTPUTS:
%
%   Yfit:  the model-constrained curve that best fits the measeurement
%
%   err:   the mean-squared fit error.
%
% MAYBE/SOMEDAY: output an inverse laplace spectrum and model coefficients

function [Yfit, err] = puprModelFit( Ymeas, T, IRFwidth, ...
    minDecayTime, maxDecayTime )

    method = 'least squares';
    %method = 'noise-informed';  % minimizes difference between fit error and estimated noise variance
    
    nT = length(T);

    % generate model timebase
    % pick sampling time that exceeds the measured T to ensure the
    % convolution is valid at the endpoints of T
    tmin = min(T) - 2*IRFwidth;
    tmax = max(T) + 2*IRFwidth;
    nt = 1025;
    t = linspace(tmin, tmax, nt);

    % ensure we have a zero point
    % (otherwise the delta functions at t=0 won't show up)
    [~, mni] = min(abs(t));
    t = t - t(mni);

    % generate the non-convolved model
    n = 32;     % number of different exponential rates in the model
                % (the final model will contain 2*n + 2 curves)

    X = zeros( 2*n+2, nt );

    % generate instantaneous responses (delta functions)
    % the factor of 10.7 is so that these have comparable 
    % magnitude to the exponential decays fter the convolution
    X(1,:) = 1*(t == 0) * 10.7;
    X(2,:) = -1*(t == 0) * 10.7;

    % generate exponential decay responses
    
    % set of decay constants to use
    tau = logspace(log10(minDecayTime),log10(maxDecayTime),n);   
    % for each decay constant, generate a positive- and a negative-
    % signed response
    for ii = 1:n
        X(ii+2,:) = exp(-t/tau(ii)).*(t>=0);
        X(ii+2+n,:) = -exp(-t/tau(ii)).*(t>=0);
    end

    % generate instrument response function matrix for convolution
    % we set up the convolution lags to match the measurement timebase
    h = zeros(nt, nT);
    for ii = 1:nT
        h(:,ii) = exp(-4*log(2)*(t-T(ii)).^2/ IRFwidth^2);
    end

    %plot(t,h);

    % convolve model with instrument response
    Xconvd = X*h;

    % add in offsets to the basis
    Xconvd(end+1,:) = 1;
    Xconvd(end+1,:) = -1;

    % SVD to generate orthonormal basis
    [U,S,V] = svd(Xconvd,0);
    
    % QR decomp to generate orthogonal basis
    [Q,R] = qr(Xconvd,0);
    
    % keep the top 10 components (should be sufficient for our signals)
    %B = V(:,1:2);
    B = V(:,1:(20));
    %B = V;
    %B = Q.';
    assignin('base','B',B);

    % decompose measured data using this basis
    switch method
        case 'least squares'
            C = Ymeas*B;
            assignin('base', 'C', C);
            %Yfit = C*pinv(B);
            Yfit = C*B.';
            Yfit = (B*(pinv(B)*Ymeas.')).';
        case 'noise-informed'
            % do this on a per-pixel basis for now
            [npx,nt] = size(Ymeas);
            
            % estimate RMSE from the first measurement  (negative probe
            % delay)
            RMSE_est = sum( Ymeas(:,1).^2 ) / npx;
            
            % come up with reasonable starting guess from least squares
            C_lsq = Ymeas*B;
            
            c0 = zeros( 1, 20 );
            
            parfor ii = 1:npx
                %c0 = C_lsq(ii, :);
                
                fnModel = @(c) c*pinv(B);
                fnRMSE = @(c) sqrt( sum( ( Ymeas(ii,:) - fnModel(c) ).^2 ) ) / nt;
                
                % a good model producel RMSE error equal to the estimated
                % RMSE
                fnFit = @(c) abs( fnRMSE(c) - RMSE_est );
                
                c = fminsearch( fnFit, c0 );
                
                C(ii,:) = c;
            end
            
            Yfit = C*B.';
    end

    err = sum((Ymeas - Yfit).^2) / nT;
end