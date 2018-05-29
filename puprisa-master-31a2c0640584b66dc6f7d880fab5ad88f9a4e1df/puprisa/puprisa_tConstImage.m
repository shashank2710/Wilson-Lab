function puprisa_tConstImage( X, delays, header, Q0Pool )
    % multiexponential fit to each pixel
    % Q0Pool is a pool of starting points for nonlinear least sq. fitting
    
    [nrow,ncol,nslices] = size(X);
    [nTrials,~] = size(Q0Pool);
    
    %pause;
    
    f = figure();

    % arrays for storing the results
    tConstIm = zeros(nrow,ncol);
    imT1 = zeros(nrow,ncol);
    imT2 = zeros(nrow,ncol);
    imA1 = zeros(nrow,ncol);
    imA2 = zeros(nrow,ncol);
    imOff = zeros(nrow,ncol);
    imChi2 = zeros(nrow,ncol)+inf;
        
    % mask of time points
    tmin = 0.5;
    tmask = (delays >= tmin);
    tmasked = delays(tmask == 1);
    %tmasked = tmasked - min(tmasked); % set first delay to zero
    t = tmasked;
    

    nt = length(tmasked);
    

    % function to which we fit
    % takes as search parameter Q:
    % Q = [amplitude1, amplitude2, timeconst1, timeconst2, offset]
    fn = @(Q,t) Q(1)*exp(-t./Q(3)) + Q(2)*exp(-t./Q(4))+Q(5);
    
    for irow = 1:nrow
        parfor icol = 1:ncol
    
            % get delay scan this pixel
            dscan = squeeze(X(irow,icol,:)).';
            %dscan = squeeze(X(302,422,:)).';
            
            % get only delay scan points after specified 'zero' delay
            dscan = dscan(tmask == 1);
            y = dscan;
            % biexp. fit

            chi2Best = inf;
            QBest = [];
            
            for iTrial = 1:nTrials
                Q0 = Q0Pool(iTrial,:);
                [Q,chi2]= lsqcurvefit(fn, Q0, t, y);

                %ee = fn(AA, x);

                if ( chi2 < chi2Best )
                    chi2Best = chi2;
                    QBest = Q;
                end
            end
            Q = QBest;
            A = [Q(1), Q(2)];
            tau = [Q(3),Q(4)];
            off = Q(5);
            % store the fit only if better than before
            [minTau, mni] = min(tau);
            [maxTau, mxi] = max(tau);
            imT1(irow,icol) = minTau;
            imT2(irow,icol) = maxTau;
            imA1(irow,icol) = A(mni);
            imA2(irow,icol) = A(mxi);
            imOff(irow,icol) = off;
            imChi2(irow,icol) = chi2;   % error

            % get mean tau
            %meanTau = sum(A0.*tau0)./sum(A0);
            meanTau = mean(tau);
            
            tConstIm(irow, icol) = meanTau;
            if 0
                fit = fn(AA,x);
                plot( tmasked, dscan, tmasked, fit );
                drawnow;
            end
            
            
        end
        
        %waitbar( irow/nrow, h );

        if 1
            subplot(221);
            imagesc(imT1);
            subplot(222);
            imagesc(imT2);
            subplot(223);
            imagesc(imA1);
            subplot(224);
            imagesc(imA2);
            drawnow;
        end
        
        % save old values
        imT1_old = imT1;
        imT2_old = imT2;
        imA1_old = imA1;
        imA2_old = imA2;
        imOff_old = imOff;
        imChi2_old = imChi2;
    end
    %close(h);
    
    imagesc(tConstIm);
    
    setappdata(gcf,'imT1', imT1);
    setappdata(gcf,'imT2', imT2);
    setappdata(gcf,'imA1', imA1);
    setappdata(gcf,'imA2', imA2);
end