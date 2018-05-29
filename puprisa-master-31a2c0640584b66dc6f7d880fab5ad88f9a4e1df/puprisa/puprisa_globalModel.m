 function varargout = puprisa_globalModel(varargin)
% PUPRISA_GLOBALMODEL M-file for puprisa_globalModel.fig
%      PUPRISA_GLOBALMODEL, by itself, creates a new PUPRISA_GLOBALMODEL or raises the existing
%      singleton*.
%
%      H = PUPRISA_GLOBALMODEL returns the handle to a new PUPRISA_GLOBALMODEL or the handle to
%      the existing singleton*.
%
%      PUPRISA_GLOBALMODEL('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PUPRISA_GLOBALMODEL.M with the given input arguments.
%
%      PUPRISA_GLOBALMODEL('Property','Value',...) creates a new PUPRISA_GLOBALMODEL or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before puprisa_globalModel_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to puprisa_globalModel_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help puprisa_globalModel

% Last Modified by GUIDE v2.5 06-Mar-2012 16:32:43

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @puprisa_globalModel_OpeningFcn, ...
                   'gui_OutputFcn',  @puprisa_globalModel_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before puprisa_globalModel is made visible.
function puprisa_globalModel_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to puprisa_globalModel (see VARARGIN)

% Choose default command line output for puprisa_globalModel
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes puprisa_globalModel wait for user response (see UIRESUME)
% uiwait(handles.figure1);

if length(varargin) == 2
    % parse t, X
    setappdata(hObject, 't', varargin{1});
    setappdata(hObject, 'X', varargin{2});
    
    disp('Data received from PUPRISA');
 
    initAxes(hObject,handles, varargin{2})
    
end

function initAxes(fig, handles, X )
    sumImage = (squeeze(sum(X.^2,3)));
    
    axes(handles.ax_chi2);
    handles.im_chi2 = image('cdata',sumImage,'tag','im_chi2','cdatamapping','scaled');
    set(gca,'climmode','auto');
    axis image;
    
    axes(handles.ax_T1);
    handles.im_T1 = image('cdata',sumImage,'tag','im_T1','cdatamapping','scaled');
    set(gca,'climmode','auto');
    axis image;
    
    axes(handles.ax_T2);
    handles.im_T2 = image('cdata',sumImage,'tag','im_T2','cdatamapping','scaled');
    set(gca,'climmode','auto');
    axis image;
    
    axes(handles.ax_A1);
    handles.im_A1 = image('cdata',sumImage,'tag','im_A1','cdatamapping','scaled');
    set(gca,'climmode','auto');
    axis image;
    
    axes(handles.ax_A2);
    handles.im_A2 = image('cdata',sumImage,'tag','im_A2','cdatamapping','scaled');
    set(gca,'climmode','auto');
    axis image;
    
    axes(handles.ax_B);
    handles.im_B = image('cdata',sumImage,'tag','im_B','cdatamapping','scaled');
    set(gca,'climmode','auto');
    axis image;
    
    axes(handles.ax_C);
    handles.im_C = image('cdata',sumImage,'tag','im_C','cdatamapping','scaled');
    set(gca,'climmode','auto');
    axis image;
    
    axes(handles.ax_U);
    handles.im_U = image('cdata',sumImage,'tag','im_U','cdatamapping','scaled');
    set(gca,'climmode','auto');
    axis image;
    
    colormap(gray);
    
    guidata(gcbf, handles);


% --- Outputs from this function are returned to the command line.
function varargout = puprisa_globalModel_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pbStart.
function pbStart_Callback(hObject, eventdata, handles)
% hObject    handle to pbStart (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles = guihandles(gcbf);

t = getappdata(gcbf,'t').';
X = getappdata(gcbf,'X');

% first generate an initial guess for the impulse response
tfwhm = 0.25;
t0 = 0;
I = exp(-4*log(2)*((t-t0)/tfwhm).^2);
 
dt_min = min(diff(t));
tspan = max(t) - min(t);
ntt  = round(tspan / dt_min);
tt = linspace(min(t),max(t),ntt).';





deltaFn = deltaFunction( tt );
stepFn = (tt>=0);

axes(handles.ax_chi2_iter);
cla;
l_chi2_iter = line('tag','l_chi2_iter','color','k');
chi2_iter = [];
l_chi2_min_iter = line('tag','l_chi2_min_iter','color','b');
chi2_min_iter = [];
l_chi2_max_iter = line('tag','l_chi2_max_iter','color','r');
chi2_max_iter = [];
legend('mean','min','max');

% set up plot of impulse response
l_I = findobj(handles.ax_I,'tag','l_I');
if isempty(l_I)
    %a = gca;
    axes(handles.ax_I);
    l_I = line();
    %axes(a);
end
set(l_I,'xdata',t,'ydata',I,'tag','l_I');

% set up plot of current pixel data point
l_y = findobj(handles.ax_currentLine, 'tag', 'l_y');
if isempty( l_y )
    axes(handles.ax_currentLine)
    l_y = line('color','b','marker','o','linestyle','none','tag','l_y',...
        'xdata',t,'ydata',zeros(1,length(t)));
end

% set up plot of current fitted point
l_yFit = findobj(handles.ax_currentLine, 'tag', 'l_yFit');
if isempty( l_yFit )
    axes(handles.ax_currentLine)
    l_yFit = line('color','k','linestyle','-','tag','l_yFit',...
        'xdata',t,'ydata',zeros(1,length(t)));
end

% set up plot of worst fit data point
l_yWorst = findobj(handles.ax_worstFit, 'tag', 'l_yWorst');
if isempty( l_yWorst )
    axes(handles.ax_worstFit)
    l_yWorst = line('color','b','marker','o','linestyle','none','tag','l_yWorst',...
        'xdata',t,'ydata',zeros(1,length(t)));
end

% set up plot of worst fit
l_yFitWorst = findobj(handles.ax_worstFit, 'tag', 'l_yFitWorst');
if isempty( l_yFitWorst )
    axes(handles.ax_worstFit)
    l_yFitWorst = line('color','k','linestyle','-','tag','l_yFitWorst',...
        'xdata',t,'ydata',zeros(1,length(t)));
end

setappdata(gcf,'doCancel',0);

[nr,nc,nt] = size(X);

% initialize fit model
% parameters: b, a1, a2, t1, t2, c
lb = [-Inf, -Inf, -Inf, 0, 0, -Inf];
ub = [Inf, Inf, Inf, Inf, Inf Inf];

% initialize results images
B = zeros(nr,nc);
A1 = zeros(nr,nc);
A2 = zeros(nr,nc);
T1 = zeros(nr,nc);
T2 = zeros(nr,nc);
C = zeros(nr,nc);
U = zeros(nr,nc);
Chi2 = zeros(nr,nc)+Inf;
negEntr = zeros(nr,nc)+Inf;
S = zeros(nr,nc);

%   p0 = [-1, 0.2, 0.5, 0.4, 3, -.08];
% B0 = B;
% A10 = A1;
% A20 = A2;
% T10 = T1;
% T20 = T2;
% C0 = C;
B0 = zeros(nr,nc)-1;
A10 = zeros(nr,nc)+.2;
A20 = zeros(nr,nc)+.5;
T10 = zeros(nr,nc)+.4;
T20 = zeros(nr,nc)+3;
C0 = zeros(nr,nc)-.08;
U0 = zeros(nr,nc);


thresh = 0;

o = optimset('TolFun', 1e-6, 'TolX', 1e-6, 'MaxFunEvals',250,'display','off',...
    'algorithm','active-set','funvalcheck','off','largeScale','on');
osa = saoptimset('display','off');

chi2_est_old = 0;
chi2_est = .001;

h_im_B = handles.im_B;

iter = 0;
bumpMethod = 1;
anneal_t0 = .1;
anneal_t = anneal_t0;
anneal_cool = 0.999;

    chi2_est =  0.015339;
    
%while( abs(chi2_est_old - chi2_est)/chi2_est_old >.0001 )
while anneal_t > (anneal_t0*1e-5)
    
    iter = iter + 1;
    chi2_est = chi2_est*0.99;
    % impulse response for each delay
    II = zeros(ntt,ntt);
    for ik = 1:ntt
        % calculte impulse response at this lag
        II(ik,:) = exp(-4*log(2)*((tt-t0-tt(ik))/tfwhm).^2);
    end

    
    delayScans = reshape(X,nr*nc, nt);

    % pick which pixels to adjust
    if iter == 1
        jj = 1:(nr*nc);
        thisDelayScans = delayScans;
        thisB = B;
        thisA1 = A1;
        thisA2 = A2;
        thisT1 = T1;
        thisT2 = T2;
        thisC = C;
        thisU = U;
        thisChi2 = Chi2;
        thisNegEntr = negEntr;
        thisB0 = B0;
        thisA10 = A10;
        thisA20 = A20;
        thisT10 = T10;
        thisT20 = T20;
        thisC0 = C0;
        thisU0 = U0;
        
        thisNpx = nr*nc;
    else
        if rand(1,1) > 0.5
           [~,jj] = sort(Chi2(:),'descend');
        else
            % pick randomly, weight by mis-fit
            r = rand(size(Chi2));

            % now weight rand by the fit err
%             if iter > 2
%                 r = r .* abs(Chi2 - chi2_est);
%             else
%                 r = r .* Chi2;
%             end

            % now select the worst fitting (highest chi2)
            [~,jj] = sort(r(:),'descend');
        end

        thisNpx = 8*10;
        jj = jj(1:thisNpx).';
        thisDelayScans = delayScans(jj,:);
        thisB = B(jj);
        thisA1 = A1(jj);
        thisA2 = A2(jj);
        thisT1 = T1(jj);
        thisT2 = T2(jj);
        thisC = C(jj);
        thisU = U(jj);
        thisChi2 = Chi2(jj);
        thisB0 = B0(jj);
        thisA10 = A10(jj);
        thisA20 = A20(jj);
        thisT10 = T10(jj);
        thisT20 = T20(jj);
        thisC0 = C0(jj);
        thisU0 = U0(jj);
        
    end
    
    % make initial guess from previous round, avgd with neighbors
    if iter > 1
        if rand(1,1) < 0.5
             % every other iter, construct guess from neighbors
             %k = ones(9,9);
             %k(2,2) = 0;
             if 0
                 %k = rand(9,9) < 1/9;
                 k = ones(9,9);
                 %k(5,5) = 0;
                 %if( sum(k(:)) == 0 )
                 %    k = ones(9,9);
                 %end
                 k = k / sum(k(:));

                 B0 = conv2(B,k,'same');% + randn(nr,nc)*0.1;
                 A10 = conv2(A1,k,'same');% + randn(nr,nc)*0.1;
                 A20 = conv2(A2,k,'same');% + randn(nr,nc)*0.1;
                 T10 = abs(conv2(T1,k,'same'));% .* (1+randn(nr,nc)*0.1));
                 T20 = abs(conv2(T2,k,'same'));% .* (1+randn(nr,nc)*0.1));
                 C0 = conv2(C,k,'same');% + randn(nr,nc)*0.1;
                 U0 = conv2(U,k,'same');
             else
                 % construct a new guess from the best of each pixel's
                 % neighbors
                 [thisB0,thisA10,thisA20,thisT10,thisT20,thisC0,thisU0] = ...
                     newGuessFromNeighbors(B,A1,A2,T1,T2,C,U, Chi2,jj);
             end
        else
            % every other iter, construct guess by perturbing previous
            % results
            thisB0 = thisB + randn(1,thisNpx)*anneal_t;
            thisA10 = thisA1 + randn(1,thisNpx)*anneal_t;
            thisA20 = thisA2 + randn(1,thisNpx)*anneal_t;
            thisT10 = thisT1 .* abs((1+randn(1,thisNpx)*anneal_t));
            thisT20 = thisT2 .* abs((1+randn(1,thisNpx)*anneal_t));
            thisC0 = thisC + randn(1,thisNpx)*anneal_t;
            thisU0 = thisU + randn(1,thisNpx)*anneal_t;
            
            
            anneal_t = anneal_t*anneal_cool;

        end
    else
        % random initial guess
        B0 = randn(nr,nc)*anneal_t;
        A10 = randn(nr,nc)*anneal_t;
        A20 = randn(nr,nc)*anneal_t;
        T10 = abs((1+randn(nr,nc)*anneal_t));
        T20 = abs((1+randn(nr,nc)*anneal_t));
        C0 = randn(nr,nc)*anneal_t;
        U0 = randn(nr,nc)*anneal_t;
        
        % fixed 'reasonable' initial guess
%         B0 = zeros(nr,nc)-1;
%         A10 = zeros(nr,nc)+.2;
%         A20 = zeros(nr,nc)+.5;
%         T10 = zeros(nr,nc)+.4;
%         T20 = zeros(nr,nc)+3;
%         C0 = zeros(nr,nc)-.08;
    end
    
    disp('starting parfor loop');
    tic;
    
    minCriteriaMethod = get(handles.popup_minCriteria, 'value');
    
    parfor ii = 1:thisNpx %ii = jj % 1:(nr*nc)   
        if 0
        % randomly set exponential amplitudes to zero
        rrr = rand(1,1);
        if rrr < 0.33
            thisA10(ii) = 0;
        elseif rrr< 0.66
            thisA20(ii) = 0;
        end
        end
        
        if 1% (Chi2(ii) > chi2_est) || (iter == 1) % only optimize the worst pixels
            % get delay scan for this pixel
            y = thisDelayScans(ii,:).';

            % set up initial guess for optimization
            p0 = [thisB0(ii),thisA10(ii),thisA20(ii),...
                    thisT10(ii),thisT20(ii),thisC0(ii),thisU0(ii)];

            % find parameters that minimize squared error of fit
            %p = fmincon(@(p)fitFnFast(y,p,t,I,imp),p0,[],[],[],[],lb,ub,[],o);
            switch minCriteriaMethod
                case 1
                    p = fminsearch(@(p)fitFn(y,p,t,tt,II,deltaFn,stepFn),p0,o);
                    %p = fminunc(@(p)fitFn(y,p,t,tt,II,deltaFn,stepFn),p0,o);
                    %p = simulannealbnd(@(p)fitFn(y,p,t,tt,II,deltaFn,stepFn),p0,[],[],osa);
                case 2
                    p = fminsearch(@(p)...
                        abs(chi2_est - fitFn(y,p,t,tt,II,deltaFn,stepFn)),...
                        p0,o);
                case 3
                    %p = fminsearch(@(p)fitFn_maxEntr(y,p,t,tt,II,deltaFn,stepFn),p0,o);
                    %p = fmincon(@(p)fitFn_maxEntr(y,p,t,tt,II,deltaFn,stepFn),p0,[],[],[],[],[],[],@(p)mycon(y,p,t,tt,II,deltaFn,stepFn),o);
                    p = fmincon(@(p)sum(p.^2),p0,[],[],[],[],[],[],@(p)mycon(y,p,t,tt,II,deltaFn,stepFn),o);
            end

            % recalculate square error (redundant; can probably pull this out
            % of fmincon results directly)
            y_fit = modelFn(p,t,tt,II,deltaFn,stepFn);
            chi2 = sum((y - y_fit).^2) / (nt);
            negEntr = fitFn_maxEntr(y,p,t,tt,II,deltaFn,stepFn);

            % sort t1, t2 to keep t1 the shorter time constant
            a1 = p(2);
            a2 = p(3);
            t1 = p(4);
            t2 = p(5);
            if t1 > t2
                % swap values
                tSwap = t1;
                t1 = t2;
                t2 = tSwap;

                aa = a1;
                a1 = a2;
                a2 = aa;
            end

            % save the resulting parameters only if they yield the best fit so
            % far
            keepThis = 0;
            
            
            if iter == 1
                keepThis = 1;
            else
                if minCriteriaMethod == 3
                    if negEntr < thisNegEntr(ii)
                        keepThis = 1;
                    end
                    keepThis = 1;
                else
                    chi2diff = abs((abs(chi2 - chi2_est) - abs(thisChi2(ii) - chi2_est))) / abs(chi2 - chi2_est);

                    if 0 % chi2diff < 0.1
                        % retain only those that minimize amplitudes
                        oldAmp = abs(thisB(ii)) + abs(thisA1(ii)) + abs(thisA2(ii)) + abs(thisC(ii)) + abs(thisU(ii));
                        thisAmp = abs(a1)+abs(a2)+abs(p(1))+abs(p(6))+abs(p(7));

                        if thisAmp < oldAmp
                            keepThis = 1;
                        end
                    else
                        if( minCriteriaMethod == 1 )
                            if chi2 < thisChi2(ii)
                                keepThis = 1;
                            end
                        else
                            if iter > 2
                                if( abs(chi2 - chi2_est) < abs(thisChi2(ii) - chi2_est) )
                                    keepThis = 1;
                                end
                            else
                                if chi2 < thisChi2(ii)
                                    keepThis = 1;
                                end
                            end
                        end
                    end
                end
            end
            % keepThis = 1;
%             if chi2 > .050
%                 chi2
%                 chi2_est_old
%                 thisChi2(ii)
%                 error('foo');
%             end

            
            if keepThis % chi2 < thisChi2(ii)
                thisB(ii) = p(1);
                thisA1(ii) = a1;
                thisA2(ii) = a2;
                thisT1(ii) = t1;
                thisT2(ii) = t2;
                thisC(ii) = p(6);
                thisU(ii) = p(7);
                thisChi2(ii) = chi2;
                thisNegEntr(ii) = negEntr;
            else
                % otherwise, re-use the old values
                % (this is redundant, but parfor doesn't work if we don't do it
                % this way)
                thisB(ii) = thisB(ii);
                thisA1(ii) = thisA1(ii);
                thisA2(ii) = thisA2(ii);
                thisT1(ii) = thisT1(ii);
                thisT2(ii) = thisT2(ii);
                thisC(ii) = thisC(ii);
                thisU(ii) = thisU(ii);
                thisChi2(ii) = thisChi2(ii);
                thisNegEntr(ii) = thisNegEntr(ii);
            end
        end
    end   % parfor
    
    % Take the results from optimizing fit of a few points, and store them
    % in the whole-image results.
    B(jj) = thisB;
    A1(jj) = thisA1;
    A2(jj) = thisA2;
    T1(jj) = thisT1;
    T2(jj) = thisT2;
    C(jj) = thisC;
    U(jj) = thisU;
    Chi2(jj) = thisChi2;
    B0(jj) = thisB0;
    A10(jj) = thisA10;
    A20(jj) = thisA20;
    T10(jj) = thisT10;
    T20(jj) = thisT20;
    C0(jj) = thisC0;
    U0(jj) = thisU0;
    
    t_elapsed = toc;
    t_per_px = t_elapsed / (nr*nc);
    fprintf('Done with parfor loop.\n\tTotal time: %f sec.\n\tAvg time per pixel: %f\n',t_elapsed, t_per_px);
    
    set(handles.im_B,'cdata',B);
    set(handles.im_C,'cdata',C);
    set(handles.im_U,'cdata',U);
    set(handles.im_A1,'cdata',A1);
    set(handles.im_A2,'cdata',A2);
    set(handles.im_T1,'cdata',T1);
    set(handles.im_T2,'cdata',T2);
    set(handles.im_chi2,'cdata',Chi2)

    % now try to optimize impulse response on the best chi2 fit so far
    % (eventually do this on a few more data points)
    %
    % This seems to destabilize the algorithm
    % SNR
    %S = A1.^2+A2.^2+B.^2+ C.^2;

    % find average fit
    [chi2_rowmins,ir]=min(Chi2);
    [chi2_min,ic] = min(chi2_rowmins);
    ir = ir(ic);
    
    % find worst fit
    [chi2_rowmaxs,ir_worst] = max(Chi2);
    [chi2_max, ic_worst] = max(chi2_rowmaxs);
    ir_worst = ir_worst(ic_worst);

    chi2_est_old = chi2_est;
    chi2_est = mean(Chi2(:));
    
    disp(['estimated chi2: ', num2str(chi2_est)]);
    chi2_iter = [chi2_iter, chi2_est];
    chi2_max_iter =[chi2_max_iter, max(Chi2(:))];
    chi2_min_iter =[chi2_min_iter, min(Chi2(:))];
    set(l_chi2_iter,'xdata',1:length(chi2_iter),'ydata',chi2_iter);
    set(l_chi2_min_iter,'xdata',1:length(chi2_min_iter),'ydata',chi2_min_iter);
    set(l_chi2_max_iter,'xdata',1:length(chi2_max_iter),'ydata',chi2_max_iter);

    y = squeeze(X(ir,ic,:)).';
    set(l_y,'ydata',y);
    
    yWorst = squeeze(X(ir_worst,ic_worst,:)).';
    set(l_yWorst,'ydata',yWorst);

    % recall parameters for this best-fitting pixel
    p(1) = B(ir,ic);
    p(2) = A1(ir,ic);
    p(3) = A2(ir,ic);
    p(4) = T1(ir,ic);
    p(5) = T2(ir,ic);
    p(6) = C(ir,ic);
    p(7) = U(ir,ic);

    y_fit = modelFn(p,t,tt,II,deltaFn,stepFn);
    set(l_yFit,'ydata',y_fit);
    
    
    % recall parameters for worst-fitting pixel
    pWorst(1) = B(ir_worst,ic_worst);
    pWorst(2) = A1(ir_worst,ic_worst);
    pWorst(3) = A2(ir_worst,ic_worst);
    pWorst(4) = T1(ir_worst,ic_worst);
    pWorst(5) = T2(ir_worst,ic_worst);
    pWorst(6) = C(ir_worst,ic_worst);
    pWorst(7) = U(ir_worst,ic_worst);

    y_fitWorst = modelFn(pWorst,t,tt,II,deltaFn,stepFn);
    set(l_yFitWorst,'ydata',y_fitWorst);
    
    drawnow;

    I0 = I;

    %tfwhm = lsqcurvefit( @(tfwhm,t)modelFn(p,t,tfwhm,imp), tfwhm, t, y, 0, Inf, o );
    %tfwhm = lsqnonlin(@(tfwhm)fitFn(y,p,t,tfwhm,imp,chi2_est),tfwhm,0,Inf,o);
    tfwhm0 = tfwhm;
    t00 = t0;
    %tfwhm = fminsearch(@(x_tfwhm)fitTotalImage(X,B,A1,A2,T1,T2,C,t,x_tfwhm,imp),tfwhm0,o);
    if rand(1,1) > 0.5
        yy = y;
        pp = p;
    else
        yy = yWorst;
        pp = pWorst;
    end
                        %p = fminsearch(@(p)fitFn(y,p,t,tt,II,deltaFn,stepFn),p0,o);

    %tparms = fminsearch(@(x)fitFn_slow(yy.',pp,t,tt,x(1),x(2),deltaFn,stepFn),[tfwhm0,t00],o);
    %tfwhm = (0.9*tfwhm+0.1*abs(tparms(1)));
    %t0 = (0.9*t0+0.1*tparms(2));
    
    setappdata(gcbf,'tfwhm',tfwhm);
    setappdata(gcbf,'t0',t0);
    
    disp(['estimated tfwhm = ', num2str(tfwhm),'; offset t0 = ', num2str(t0)]);
    I = exp(-log(2)*(t/tfwhm).^2);

    y_fit = modelFn(p,t,tt,II,deltaFn,stepFn);
    set(l_yFit,'ydata',y_fit);

    set(l_I, 'ydata',I);
    drawnow;
    
end % for iter

function [B0,A10,A20,T10,T20,C0,U0] = ...
    newGuessFromNeighbors(B,A1,A2,T1,T2,C,U, Chi2, jj)
% for each pixel, take the best-fitting neighbors, and use their parameters
% for this pixel's next guess


    [nr,nc] = size(Chi2);
    
    methodNum = ceil(rand(1,1)*4);
    
    switch methodNum
        case 1
            maskMethod = 'nearest';
        case 2
            maskMethod = 'randomBest';
        case 3
            maskMethod = 'roulette';    % pick with weight by chi2
        case 4
            maskMethod = 'oneRandom';
    end
    
    maskHeight = 3;
    maskWidth  = 3;
    
    npx = length(jj);
    B0 = zeros(1,npx);
    A10 = zeros(1,npx);
    A20 = zeros(1,npx);
    T10 = zeros(1,npx);
    T20 = zeros(1,npx);
    C0 = zeros(1,npx);
    U0 = zeros(1,npx);
    
    %for ir = 1:nr
    %    parfor ic = 1:nc
    for ii = jj
        if strcmp(maskMethod,'roulette')
            r = rand(size(Chi2));
            
            % now weight rand by the fit err
            r = r .* Chi2;
            
            % now select the lowest
            [~,iMinChi2] = min(r(:));
        elseif strcmp(maskMethod,'oneRandom')
            iMinChi2 = floor(rand(1,1)*(nr*nc))+1;
        else
            ir = mod(ii-1,nr)+1;
            ic = floor((ii-1)/nr)+1;
            switch maskMethod
                case 'nearest'
                    if ir == 1
                        maskRows = ir:(ir+floor(maskHeight/2));
                    elseif ir == nr
                        maskRows = (ir-floor(maskHeight/2)):ir;
                    else
                        maskRows = (ir - floor(maskHeight/2)) : (ir + floor(maskHeight/2));
                    end

                    if ic == 1
                        maskCols = ic:(ic+floor(maskWidth/2));
                    elseif ic == nc
                        maskCols = (ic-floor(maskWidth/2)):ic;
                    else
                        maskCols = (ic - floor(maskWidth/2)) : (ic + floor(maskWidth/2));
                    end

                    % construct a mask
                    mask = zeros(nr,nc);
                    mask(maskRows, maskCols) = 1;
                case 'randomBest'
                    mask = (rand(nr,nc) < 0.01);
            end

            % then mask Chi2
            thisChi2 = Chi2.*mask + Inf.*(1-mask);

            % now find index of best-fitting neighbors
            [~, iMinChi2] = min(thisChi2(:));
        end
        % construct new guess from this result
        B0(jj)  = B(iMinChi2);
        A10(jj) = A1(iMinChi2);
        A20(jj) = A2(iMinChi2);
        T10(jj) = T1(iMinChi2);
        T20(jj) = T2(iMinChi2);
        C0(jj)  = C(iMinChi2);
        U0(jj) = U(iMinChi2);
        
    end
    
    % reshape to return only the selected pixels
    B0 = B0(jj);
    A10 = A10(jj);
    A20 = A20(jj);
    T10 = T10(jj);
    T20 = T20(jj);
    C0 = C0(jj);
    U0 = U0(jj);

function chi2 = fitFn(y,p,t,tt,II,deltaFn,stepFn)
    
    if( (p(4) < 0) || (p(5) < 0) )
        % check for disallowed negative time constants
        % (these would indicate unphysical exponential signal growth)
        % If these are found, return large fit error, and don't bother
        % evaluating the model function.
        chi2 = 1e12;
    else
        % returns fit rmse
        y_fit = modelFn(p,t,tt,II,deltaFn,stepFn);

        chi2 = sum( (y - y_fit).^2 ) / length(t);


        if isinf(chi2)
            chi2 = 1e12;
        end
        if isnan(chi2)
            chi2 = 1e12;
        end
    end
    
function chi2 = fitFn_slow(y,p,t,tt,tfwhm,t0,deltaFn,stepFn)
    % recaclulate II depending on tfwhm, t0
    
    ntt = length(tt);
    
    % impulse response for each delay
    II = zeros(ntt,ntt);
    for ik = 1:ntt
        % calculte impulse response at this lag
        II(ik,:) = exp(-4*log(2)*((tt-t0-tt(ik))/tfwhm).^2);
    end

    chi2 = fitFn(y,p,t,tt,II,deltaFn,stepFn);
    
function negEntr = fitFn_maxEntr(y,p,t,tt,II,deltaFn,stepFn)
    
    if( (p(4) < 0) || (p(5) < 0) )
        % check for disallowed negative time constants
        % (these would indicate unphysical exponential signal growth)
        % If these are found, return large fit error, and don't bother
        % evaluating the model function.
        negEntr = 1e12;
    else
        % returns fit rmse
        y_fit = modelFn(p,t,tt,II,deltaFn,stepFn);

        res = (y - y_fit);
        %res = res - mean(res);
        negEntr = sum(res.*log(res));
    end
    
function [c,ceq] = mycon(y,p,t,tt,II,deltaFn,stepFn)
     if( (p(4) < 0) || (p(5) < 0) )
        % check for disallowed negative time constants
        % (these would indicate unphysical exponential signal growth)
        % If these are found, return large fit error, and don't bother
        % evaluating the model function.
        chi2 = 1e12;
    else
        % returns fit rmse
        y_fit = modelFn(p,t,tt,II,deltaFn,stepFn);

        chi2 = sum( (y - y_fit).^2 ) / length(t);


        if isinf(chi2)
            chi2 = 1e12;
        end
        if isnan(chi2)
            chi2 = 1e12;
        end
     end
    
     c = 0;
     ceq = chi2 - 0.001;

function y = modelFn(p, t, tt, II,deltaFn,stepFn)
% II is the impulse response, each row shifted by a different amount    
% delt is delta function, i.e. [0 0 0 0 1 0 0 0 0]
    
   
    yy = p(1)*deltaFn ...
        + p(2)*exp(-(tt.*stepFn)/p(4)).*stepFn ...
        + p(3)*exp(-(tt.*stepFn)/p(5)).*stepFn ...
        + p(7).*stepFn;
    
    % convolve with impulse function
    yy = II * yy;
    yy = yy + p(6);
    
    % instrument saturation
    %yy(yy>10) = 10; yy(yy< -10.6) = -10.6;
    
    % interpolate
    y = interp1(tt,yy,t);
    
function c = irregularConvolve(t,x2,tfwhm, t0)
% convolve functions with an irregulaly-spaced time axis
    %tfwhm = 0.25; 
    nt = length(t);
    c = zeros(nt,1);
    
   
    II = zeros(ntt,ntt);
    
    for ik = 1:ntt
        % calculte impulse response at this lag
        II(ik,:) = exp(-4*log(2)*((tt-t0-tt(ik))/tfwhm).^2);
    end

    c = II * x2;
    
    
    
function y = modelFnFast(p, t, I, delt)
    % delt is impulse, i.e. [0 0 0 0 1 0 0 0 0]
    
    % this might be made faster by precomputing t>=0

%     y = p(1)*delt ...
%         + p(2)*exp(-t/p(4)).*(t>=0) ...
%         + p(3)*exp(-t/p(5)).*(t>=0) ...
%         + p(7).*(t>=0);
%     
%     %c = conv(yy,I,'same');
%     %c = c(2:(1+length(I)));
%     %y = c + p(6);
%     
%     % instrument saturation
%     y(y>10) = 10; y(y< -10.6) = -10.6;

    y = modelFn(p, t, 0,0, delt);

    
function chi2 = fitTotalImage(X,B,A1,A2,T1,T2,C,U, t, tfwhm, imp)
    % fit entire image
    [nr,nc,~] = size(X);
    for ir = 1:nr
        for ic = 1:nc
            p = [B(ir,ic),A1(ir,ic),A2(ir,ic),T1(ir,ic),T2(ir,ic),C(ir,ic),U(ir,ic)];
            Y(ir,ic,:) = modelFn(p, t, tfwhm, imp);
        end
    end
    
    chi2 = sum( (X(:) - Y(:)).^2 );
    
    if isinf(chi2)
        chi2 = 1e12;
    end
    if isnan(chi2)
        chi2 = 1e12;
    end


% --- Executes on button press in pbPause.
function pbPause_Callback(hObject, eventdata, handles)
% hObject    handle to pbPause (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pbCancel.
function pbCancel_Callback(hObject, eventdata, handles)
% hObject    handle to pbCancel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

setappdata(gcbf,'doCancel',1);


function ed_initFWHM_Callback(hObject, eventdata, handles)
% hObject    handle to ed_initFWHM (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ed_initFWHM as text
%        str2double(get(hObject,'String')) returns contents of ed_initFWHM as a double


% --- Executes during object creation, after setting all properties.
function ed_initFWHM_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ed_initFWHM (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pbResults.
function pbResults_Callback(hObject, eventdata, handles)
% hObject    handle to pbResults (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    handles=guihandles(gcbf);
    B = get(handles.im_B,'cdata');
    A1 = get(handles.im_A1,'cdata');
    A2 = get(handles.im_A2,'cdata');
    T1 = get(handles.im_T1,'cdata');
    T2 = get(handles.im_T2,'cdata');
    C = get(handles.im_C,'cdata');
    U = get(handles.im_U,'cdata');
    props = {'color','k','xtick',[],'ytick',[]};
    
    figure;
    ax_T1 = axes('position',[0,0.5,1/3,0.5],props{:});
    imagesc('cdata',T1,'alphadata',abs(A1),'alphadatamapping','scaled');
    axis image
    
    ax_T2 = axes('position',[1/3,0.5,1/3,0.5],props{:});
    imagesc('cdata',T2,'alphadata',abs(A2),'alphadatamapping','scaled');
    axis image
    
    
    %h=linkprop([ax_T1,ax_T2],{'clim','alim'});
    %setappdata(ax_T1,'linkprop',h);

    ax_B = axes('position',[0,0,1/3,0.5],props{:});
    imagesc('cdata',B);
    axis image
    
    ax_C = axes('position',[1/3,0,1/3,0.5],props{:});
    imagesc('cdata',C)
    axis image
    
    %h2 = linkprop([ax_B,ax_C],'clim');
    
    %setappdata(ax_B,'linkprop',h2);
    
    ax_U = axes('position',[2/3,0,1/3,0.5],props{:});
    imagesc('cdata',U);
    axis image;
    
    h=linkprop([ax_B,ax_C,ax_T1,ax_T2,ax_U],{'xlim','ylim'});
    setappdata(ax_B,'linkprop',h);

% --- Executes on button press in pbResiduals.
function pbResiduals_Callback(hObject, eventdata, handles)
% hObject    handle to pbResiduals (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% get fit parameters
    handles=guihandles(gcbf);
    B = get(handles.im_B,'cdata');
    A1 = get(handles.im_A1,'cdata');
    A2 = get(handles.im_A2,'cdata');
    T1 = get(handles.im_T1,'cdata');
    T2 = get(handles.im_T2,'cdata');
    C = get(handles.im_C,'cdata');
    U = get(handles.im_U,'cdata');
    
    t = getappdata(gcbf,'t').';
    X = getappdata(gcbf,'X');
    
    % reconstruct the results
    Y = 0*X;
    [nr,nc,nt] = size(X);
    
    %tfwhm = 0.250;
    tfwhm = getappdata(gcbf,'tfwhm');
    t0 = getappdata(gcbf,'t0');
    imp = deltaFunction(t);
    
    for ir = 1:nr
        for ic = 1:nc
            p = [B(ir,ic),A1(ir,ic),A2(ir,ic),T1(ir,ic),T2(ir,ic),C(ir,ic),U(ir,ic)];
            Y(ir,ic,:) = modelFn(p, t, tfwhm,t0, imp);
        end
    end
    
    puprisa_viewChannel(Y, t, 'fitted',[],1);

% then subtract from original image
    R = X - Y;
    puprisa_viewChannel(R, t, 'residual',[],1);

% and launch puprisa_viewChannel


function d = deltaFunction( t )
% make a dirac delta function for the impulse response
    [~,it0] = min(abs(t - 0));
    d = zeros(length(t),1);
    d(it0) = 1;


% --- Executes on selection change in popup_minCriteria.
function popup_minCriteria_Callback(hObject, eventdata, handles)
% hObject    handle to popup_minCriteria (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popup_minCriteria contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popup_minCriteria


% --- Executes during object creation, after setting all properties.
function popup_minCriteria_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popup_minCriteria (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
