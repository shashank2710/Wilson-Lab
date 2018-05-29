% testModelFitFunction.m
%
% tests puprModelFit.m on some sample data
%
% Jesse Wilson (2013) jesse.wilson@duke.edu
% (See lab notebook entry, April 23rd, 2013)

%% load some data
fileName = ['\\wlserver\JesseData\2013-04-19\',...
            'EDTAWashed_Sepia_Power9_1.dat'];

[slices, header] = puprisa_readImageStack( fileName );

[T, YY, ~] = puprisa_getChannelFromSlices(slices, 1, header);
Tmin = min(T);
Tmax = max(T);
[nr,nc,nT] = size(YY);
Y = squeeze(sum(sum(YY, 1),2)) / (nr*nc);
Ymeas = Y.';    % make Y a single-row, multiple-column vector
plot(T,Y);

%% fit a model to the sample data
IRFwidth = 0.250;
[Yfit, err] = puprModelFit( Ymeas, T, IRFwidth, ...
    0.1, 30.0 );

clf;
line(T,Y,'LineStyle','none','Marker','x','color','k','MarkerSize',5.0);
line(T,Yfit,'color','b');
xlim([Tmin-0.2,Tmax+0.2]);
xlabel('probe delay, ps');
ylabel('signal (arg)');
legend('measurement','model fit');
title('Test pump-probe model fit');

%% vary the IRF width and see how it affects the fit
widths = linspace(0.1,0.4,64);
errs = zeros(1,length(widths));
Yfits = zeros(ii, nT);
for ii = 1:length(widths)
    [Yfit, err] = puprModelFit( Ymeas, T, widths(ii), ...
    0.1, 30.0 );

    errs(ii) = err;
    Yfits(ii,:) = Yfit;
end

clf;
plot(widths,errs);
xlabel('Instrument response function FWHM (ps)');
ylabel('Model fit mean-sq err');
title('Effect of estimated IRF width on model fit');

%% compare 2 different widths
clf;
[mn,mni]=min(errs);
line(T,Y,'LineStyle','none','Marker','x','color','k','MarkerSize',5.0);
line(T,Yfits(mni,:),'color','b');
line(T,Yfits(33,:),'color','r');
legend('measured data','fit with tfwhm=0.1952ps','fit with tfwhm=0.2524ps');

xlim([Tmin-0.2,Tmax+0.2]);
xlabel('probe delay, ps');
ylabel('signal (arg)');
title('Model fit for different estimated fwhm');