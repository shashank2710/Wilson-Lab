function mnuMultiexpFit(src, evt)
    % perform a multiexponential fit of the currently-displayed time series
        
    %get paramters from dialog
    fitParm = puprisa_MultiExpFitDialog();
    
    % Parameters for multi-exponential fit
    nExp  = fitParm.NoOfExp; %No. of exponentials
        
    if (nExp == 1)
        LB = [fitParm.InstRespL,fitParm.A1L,fitParm.T1L];
        UB = [fitParm.InstRespU,fitParm.A1U,fitParm.T1U];
        InCon = [fitParm.InstResponse,fitParm.A1,fitParm.T1];
    else
        %bounds for multi-exponential fit
        LB = [fitParm.InstRespL,fitParm.A1L,fitParm.T1L,fitParm.A2L,fitParm.T2L];
        UB = [fitParm.InstRespU,fitParm.A1U,fitParm.T1U,fitParm.A2U,fitParm.T2U];
        %Initial conditions for multi-exponential fit
        InCon = [fitParm.InstResponse,fitParm.A1,fitParm.T1,fitParm.A2,fitParm.T2];
    end
    
    %fit options
    options = optimset('MaxIter',12000,'TolFun',1.0e-9,'MaxFunEvals',8000);
    
       
    l = findobj(gcf,'selected','on'); % get selected object for the line
    x = get(l,'XData');
    y = get(l,'YData');
    
    xfwhm = 0.23;%cross-correlation width in ps
    x0 = 0;%delay offset
    dx_min = min(diff(x));
    min(x) - (10*dx_min),max(x) + (10*dx_min)
     abs(round(max(x) - min(x)/(0.01*dx_min)))
    xNew = linspace(min(x) - (10*dx_min),max(x) + (10*dx_min),...
        abs(round(max(x) - min(x)/(0.03*dx_min))));%create new equally spaced time axis
    xLen = length(xNew);
    [~,ix0] = min(abs(xNew - x0));
    delFunc = 1*(xNew==xNew(ix0));
    stepFunc = 1*(xNew>=x0);
    
    I = zeros(xLen,xLen);
    for ik = 1:xLen
        % calculte impulse response at this lag 
        I(ik,:) = exp(-4*log(2)*((xNew-x0-xNew(ik))/xfwhm).^2);
        I(ik,:) = I(ik,:)/sqrt(sum(I(ik,:).^2)); %unit vector
    end
    
    if (nExp == 1)
        fn = @(A,x) SingleExpFitFunc(A,x,xNew,I,delFunc,stepFunc);
        [AA,RSS,Residual,~,~,~,J] = lsqcurvefit(fn,InCon,x,y,LB,UB,options);
        ci = nlparci(AA,Residual,'jacobian',J);%95% confidence intervals
        AAerror = (ci(:,2)-ci(:,1))/(2*1.96); %standard error
        tau = sprintf('tau = %6.3f %c %6.4f',AA(3),177,AAerror(3));
        amp = sprintf('Amp = %6.3f %c %6.4f',AA(2),177,AAerror(2));
        sprintf('%s \n%s \nInst Response = %6.4f',tau,amp,AA(1))
        ee = fn(AA,x);
    else
        fn = @(A,x) multiExpFitFunc(A,x,xNew,I,delFunc,stepFunc);
        [AA,RSS,Residual,~,~,~,J] = lsqcurvefit(fn,InCon,x,y,LB,UB,options);
        ci = nlparci(AA,Residual,'jacobian',J);%95% confidence intervals
        AAerror = (ci(:,2)-ci(:,1))/(2*1.96); %standard error
        
        tau = sprintf('tau = %6.3f %c %6.4f\t %6.3f %c %6.3f',AA(3),177,AAerror(3),AA(5),177,AAerror(5));
        amp = sprintf('Amp = %6.3f %c %6.4f\t %6.3f %c %6.3f',AA(2),177,AAerror(2),AA(4),177,AAerror(4));
        sprintf('%s \n%s \nInst Response = %6.4f\n',tau,amp,AA(1))
        ee = fn(AA, x);
    end
    line(x,ee);
   
end

function y = SingleExpFitFunc(A,t,tNew,I,delFunc,stepFunc)

yy = A(1)*delFunc...
    + (A(2)*exp(-(tNew.*stepFunc)/A(3))).*stepFunc;

yyConv = (yy*I);

y = interp1(tNew,yyConv,t,'spline');

end

function y = multiExpFitFunc(A,t,tNew,I,delFunc,stepFunc)

yy = A(1)*delFunc...
    + (A(2)*exp(-(tNew.*stepFunc)/A(3))...
    + A(4)*exp(-(tNew.*stepFunc)/A(5))).*stepFunc;

yyConv = (yy*I);

y = interp1(tNew,yyConv,t,'spline');

end
