function varargout = puprisa_exponentialFit(varargin)
% PUPRISA_EXPONENTIALFIT M-file for puprisa_exponentialFit.fig
%      PUPRISA_EXPONENTIALFIT, by itself, creates a new PUPRISA_EXPONENTIALFIT or raises the existing
%      singleton*.
%
%      H = PUPRISA_EXPONENTIALFIT returns the handle to a new PUPRISA_EXPONENTIALFIT or the handle to
%      the existing singleton*.
%
%      PUPRISA_EXPONENTIALFIT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PUPRISA_EXPONENTIALFIT.M with the given input arguments.
%
%      PUPRISA_EXPONENTIALFIT('Property','Value',...) creates a new PUPRISA_EXPONENTIALFIT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before puprisa_exponentialFit_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to puprisa_exponentialFit_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES
%
%
% INTERNAL DATA STRUCTURES (appdata)

% Edit the above text to modify the response to help puprisa_exponentialFit

% Last Modified by GUIDE v2.5 31-Aug-2012 17:09:47

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @puprisa_exponentialFit_OpeningFcn, ...
                   'gui_OutputFcn',  @puprisa_exponentialFit_OutputFcn, ...
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

% --- Executes during object creation, after setting all properties.
function puprisa_exponentialFit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to puprisa_exponentialFit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% --- Executes just before puprisa_exponentialFit is made visible.
function puprisa_exponentialFit_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to puprisa_exponentialFit (see VARARGIN)
%            passing in a delay stack and a time vector will load that
%            stack initially.

% Choose default command line output for puprisa_exponentialFit
handles.output = hObject;

set(hObject,'toolbar','figure');
set(handles.sl_contrastEnhance,'Value',1);

currentfolder = pwd;
[filepath,~, ~] =fileparts( which('puprisa_exponentialFit'));
switch currentfolder 
    case filepath
    otherwise
        cd(filepath)
end

% load an image, if it was passed in the input arguments
if nargin == 5 % data was passed into the function
    setappdata(hObject,'imageStack',varargin{1});
    setappdata(hObject,'delays',varargin{2});
    setappdata(hObject,'numDelays', length(varargin{2}));
    setappdata(hObject,'fileName','puprisaInput.dat');
    % set up sliderbar and simulate slider movement
    setupSlider(handles.sl_delayIndex, length(varargin{2}));
    guidata(hObject, handles);
    sl_delayIndex_Callback(handles.sl_delayIndex, eventdata, handles);
end

% initialize fitParm
if nargin <= 5 % no fit parameter structure was passed into the function
    fitParm.FitModel = 'Double exp';
    set(handles.pm_fitModel,'Value', 2);
    fitParm.GuessMethod = 'Manual';
    set(handles.pm_initialGuess,'Value', 1);
    
    fitParm.GuessVal.Offset = 0;
    set(handles.ed_guessOffset,'String', num2str(fitParm.GuessVal.Offset));
    fitParm.GuessVal.Inst = 1;
    set(handles.ed_guessInst,'String', num2str(fitParm.GuessVal.Inst));
    fitParm.GuessVal.A1 = 1;
    set(handles.ed_guessA1,'String', num2str(fitParm.GuessVal.A1));
    fitParm.GuessVal.A2 = 0.5;
    set(handles.ed_guessA2,'String', num2str(fitParm.GuessVal.A2));
    fitParm.GuessVal.A3 = 0.2;
    set(handles.ed_guessA3,'String', num2str(fitParm.GuessVal.A3));
    fitParm.GuessVal.T1 = 0.5;
    set(handles.ed_guessT1,'String', num2str(fitParm.GuessVal.T1));
    fitParm.GuessVal.T2 = 2;
    set(handles.ed_guessT2,'String', num2str(fitParm.GuessVal.T2));
    fitParm.GuessVal.T3 = 10;
    set(handles.ed_guessT3,'String', num2str(fitParm.GuessVal.T3));
    
    fitParm.ThreshVariation = 15;
    set(handles.ed_threshVar,'String', num2str(fitParm.ThreshVariation));

    fitParm.InputVal.CrossCorrWidth = 0.25;
    set(handles.ed_crossCorrWidth,'String', num2str(fitParm.InputVal.CrossCorrWidth));
    fitParm.InputVal.DelayOffset = 0;
    set(handles.ed_delayOffset,'String', num2str(fitParm.InputVal.DelayOffset));
    fitParm.InputVal.FitStart = 0.4;
    set(handles.ed_fitStart,'String', num2str(fitParm.InputVal.FitStart));

    fitParm.BinSize = 20;
    set(handles.ed_binSize,'String', num2str(fitParm.BinSize));
    
    setappdata(hObject, 'fitParm', fitParm);
end

% initialize fitParm
if nargin == 6 % fit parameter and data structure was passed in -> fit
    error('Batch processing not yet implemented');
end

% Update handles structure
guidata(hObject, handles);


% UIWAIT makes puprisa_exponentialFit wait for user response (see UIRESUME)
% uiwait(handles.puprisa_exponentialFit);


% --- Outputs from this function are returned to the command line.
function varargout = puprisa_exponentialFit_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pb_loadImageStack.
function pb_loadImageStack_Callback(hObject, eventdata, handles)

if ispref('puprisa','imageWorkingDirectory')
    pathName = getpref('puprisa','imageWorkingDirectory');
else
    pathName = pwd();
    addpref('puprisa','imageWorkingDirectory',pathName);
end

% prompt user for imageStack file to be loaded
[fileName, pathName] = uigetfile( [pathName,filesep,'*.dat'], ...
                                  'Pick an image stack .dat file',...
                                  'multiselect','off');
% fileName = 'LapisDelayStack_10500_5.dat';
% pathName = 'D:\User\Martin\Documents\Data\Art\Pigments\Lapis\10500\';
                              
if isequal(fileName, 0); return; end

setpref('puprisa','imageWorkingDirectory',pathName);

setappdata(gcbf,'fileName',fileName);

[S,header] = puprisa_readImageStack( [pathName,filesep,fileName] );
[t,X] = puprisa_getChannelFromSlices(S,1,header);

setappdata(gcbf,'imageStack',X);
setappdata(gcbf,'delays',t);
setappdata(gcbf,'numDelays', length(t));
if isappdata(gcbf,'fitResults'); rmappdata(gcbf,'fitResults'); end

% set up sliderbar and simulate slider movement
setupSlider(handles.sl_delayIndex, length(t));
sl_delayIndex_Callback(handles.sl_delayIndex, eventdata, handles);


% sets up the slider properties
function setupSlider(handle, numDelays)

set(handle, 'Value', 1);
set(handle, 'Min', 1);
set(handle, 'Max', numDelays);
set(handle, 'SliderStep', 1/(numDelays-1)*[1 1]);



% --- Executes on slider movement.
function sl_delayIndex_Callback(hObject, eventdata, handles)
% hObject    handle to sl_delayIndex (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% instead of gcbo, use the parent of the passed handle, since I also
% call this function directly (not as a callback)
figH = get(hObject,'parent');

ind = round(get(hObject,'Value'));
setappdata(figH,'currentSlice', ind);
data = getappdata(figH,'imageStack');

% create or update image
im = findobj(figH,'Tag','rawImage');
if isempty(im) % doesn't exist yet -> create
    axes(handles.ax_imageData);
    im = imagesc(data(:,:,ind));
    set(im,'Tag','rawImage');
    axis image;
    axis off;
    colormap(cmap(256,[0,1,1],[0,0,1],[0,0,0],[1,0,0],[1,1,0]));

    % create ROI
    roi = imrect(gca, [10 10 20 20]);
    setappdata(figH,'roiHandle', roi)
    fcn = makeConstrainToRectFcn('imrect',get(handles.ax_imageData,'XLim'),get(handles.ax_imageData,'YLim'));
    setPositionConstraintFcn(roi,fcn); 
     
else % image exists -> update data only
    set(im,'CData', data(:,:,ind));
end

% --- Executes on selection change in pm_fitDisplay.
function pm_fitDisplay_Callback(hObject, eventdata, handles)
% hObject    handle to pm_fitDisplay (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

figH = get(hObject,'parent');
if ~isappdata(figH,'fitResults')
    error('Need to run grid fit first');
end

ind = get(hObject,'Value');
currDisp = getappdata(figH,'currentDisplay');

% if currentDisplay not defined yet or if it has changed, set 
if (isempty(currDisp) || (currDisp ~= ind))
    minValDisplay = -Inf;
    maxValDisplay = Inf;
    setappdata(figH,'minValDisplay', minValDisplay);
    setappdata(figH,'maxValDisplay', maxValDisplay);
    set(handles.ed_minValDisp,'String','-Inf')
    set(handles.ed_maxValDisp,'String','Inf')
else        
    minValDisplay = getappdata(figH,'minValDisplay');
    maxValDisplay = getappdata(figH,'maxValDisplay');
end    
    
setappdata(figH,'currentDisplay', ind);

fitParm = getappdata(figH,'fitParm');
fitResults = getappdata(figH,'fitResults');
[numBinsX, numBinsY, numParam] = size(fitResults);
I = getappdata(figH,'intensityImage');
fitMask = getappdata(gcbf, 'fitMask');
M = expand2DSet (fitMask, fitParm.BinSize);


contrEnh = get(handles.sl_contrastEnhance,'Value');

% create or update image
iCode = get(handles.cb_intensityCoded,'Value');
if iCode
    I = I(1:numBinsX*fitParm.BinSize,1:numBinsY*fitParm.BinSize).*M;
else
    I = M;
end
I_normd = I - min(I(:));
I_normd = I_normd/max(I_normd(:));
I_normd = I_normd.^(contrEnh);
%I_normd = round(I_normd*255);
%cm = contrast(I_normd,256);
%IRGB = ind2rgb(I_normd,cm);

C = expand2DSet (fitResults(:,:,ind), fitParm.BinSize);
if minValDisplay == -Inf
    C_normd = C - min(C(:));
else
    C_normd = C - minValDisplay;
end
if maxValDisplay == Inf
    C_normd = C_normd/max(C_normd(:));
else
    C_normd = C_normd/maxValDisplay;
end
C_normd = round(C_normd*255);

% clip under and overflow pixels
lowBound = (C_normd<0);
highBound = (C_normd>255);
C_normd(lowBound) = 0;
C_normd(highBound) = 255;

im = findobj(figH,'Tag','fitResultImage');
if isempty(im) % doesn't exist yet -> create
    axes(handles.ax_fitResults);
    im = image('CData',C_normd,'alphadata',I_normd,'alphadatamapping','scaled');
%    im = image('CData',C_normd);
    set(im,'Tag','fitResultImage');
    set(handles.ax_fitResults,'ydir','reverse');
    %axis square;
    axis image;
    set(gca,'xtick',[],'ytick',[]);
    colormap(cmap(256,[0,1,1],[0,1,0],[1,0,0],[1,1,0]));
    
    % create pointer
    pointer = impoint(gca, [10 10]);
    setappdata(figH,'pointerHandle', pointer)
    fcn = makeConstrainToRectFcn('impoint',get(gca,'XLim'),get(gca,'YLim'));
    setPositionConstraintFcn(pointer,fcn); 
    
else % image exists -> update data only
    set(im,'CData', C_normd);
    set(im,'alphadata',I_normd);
end



% --- Executes on button press in pb_fitROI.
function pb_fitROI_Callback(hObject, eventdata, handles)
% hObject    handle to pb_fitROI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% first, clear the plots
cla(handles.ax_transientData);

X = getappdata(gcbf, 'imageStack');
roiHandle = getappdata(gcbf,'roiHandle');
p = getPosition(roiHandle);

% find X within ROI box
iLower = ceil(p(1));
iUpper = floor(p(1)+p(3));
jLower = ceil(p(2));
jUpper = floor(p(2)+p(4));
Xroi = X(  jLower:jUpper,iLower:iUpper, : );
XroiSum = squeeze(sum(sum(Xroi,1),2));
lineout = (XroiSum /((jUpper-jLower)*(iUpper-iLower)))'; % normalize

% now plot averaged data
axes(handles.ax_transientData);
delays = getappdata(gcbf,'delays');
plot(delays, lineout);
hold on;

% pass data and fit parameters to the fit function
fitParm = getappdata(gcbf, 'fitParm');
[fitResults, fitRange, fitCurve] = fitExpCurve(delays, lineout, fitParm);
plot(fitRange, fitCurve, 'r');
axis tight;

line = formatFitResults(fitResults.Parameters, fitResults.Values(:,1),...
    fitResults.Values(:,2));
set(handles.tx_fitResults,'String',line);


function line = formatFitResults(names, values, errors)
numElem = length(names);
pmSign = cellstr(char(ones(numElem,1)*177));
eqSign = cellstr(char(ones(numElem,1)*61));
vals = cellstr(num2str(values));
errors = cellstr(num2str(errors));
line = strcat(names', eqSign, vals, pmSign, errors);


% --- Executes on button press in pb_checkFit.
function pb_checkFit_Callback(hObject, eventdata, handles)
% hObject    handle to pb_checkFit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% first, clear the plots
cla(handles.ax_transientData);

fitMask = getappdata(gcbf, 'fitMask');
binnedData = getappdata(gcbf,'binnedData');
binnedSize = size(binnedData);
fitCurves = getappdata(gcbf,'fitCurves');
fitParm = getappdata(gcbf,'fitParm');
fitRange = getappdata(gcbf,'fitRange');

%set(handles.ed_binSize,'String', num2str(fitParm.BinSize));


pointerHandle = getappdata(gcbf,'pointerHandle');
p = getPosition(pointerHandle);

% find the right location corresponding to pointer
% pointer returns range [0.5..max+0.5]
maxPosX = get(handles.ax_fitResults,'XLim');
maxPosY = get(handles.ax_fitResults,'YLim');
% coordinates in original image
posX = 1+round((p(1)-maxPosX(1))/maxPosX(2)*(fitParm.BinSize*binnedSize(1)-1));
posY = 1+round((p(2)-maxPosY(1))/maxPosY(2)*(fitParm.BinSize*binnedSize(2)-1));

% coordinate in binned data
binX = 1+floor((posX-1)/fitParm.BinSize);
binY = 1+floor((posY-1)/fitParm.BinSize);
%orig = squeeze(binnedData(binX, binY, :));
%fit = squeeze(fitCurves(binX, binY, :));
orig = squeeze(binnedData(binY, binX, :));
fit = squeeze(fitCurves(binY, binX, :));


% now plot fitted data
axes(handles.ax_transientData);
delays = getappdata(gcbf,'delays');
plot(delays, orig');
hold on;
plot(fitRange, fit', 'r');
% plot(fitRange, fitCurve, 'r');
axis tight;

names = getappdata(gcbf,'fitValueKey');
results = getappdata(gcbf,'fitResults');
errors = getappdata(gcbf,'fitErrors');

line = formatFitResults(names, squeeze(results(binY, binX,:)),...
    squeeze(errors(binY, binX,:)));
set(handles.tx_fitResults,'String',line);




% --- Executes on button press in pb_fitGrid.
function pb_fitGrid_Callback(hObject, eventdata, handles)
% hObject    handle to pb_fitGrid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% clear plots
cla(handles.ax_transientData);
cla(handles.ax_fitResults);

fitParm = getappdata(gcbf, 'fitParm');
X = getappdata(gcbf, 'imageStack');
delays = getappdata(gcbf,'delays');

[fitResults, fitErrors, mask, fitRange, binnedData, fitCurves, valueKey] = ...
    fitExpGrid(X, delays, fitParm, 1);

setappdata(gcbf,'fitResults', fitResults);
setappdata(gcbf,'fitErrors', fitErrors);
setappdata(gcbf,'fitMask', mask);
setappdata(gcbf,'fitRange', fitRange);
setappdata(gcbf,'binnedData', binnedData);
setappdata(gcbf,'fitCurves', fitCurves);
setappdata(gcbf,'fitValueKey', valueKey);

% set up popupmenu
set(handles.pm_fitDisplay,'String',valueKey);
set(handles.pm_fitDisplay,'Value',1);

% get abs data
I = sqrt(sum(X.^2,3));
setappdata(gcbf,'intensityImage',I);

pm_fitDisplay_Callback(handles.pm_fitDisplay, eventdata, handles);


function [fitResults, fitErrors, mask, fitRange, binnedData, fitCurves, valueKey] = ...
    fitExpGrid(data, delays, fitParm, bar)
% bin data
binnedData = bin3DSet (data, fitParm.BinSize);
newSize = size(binnedData);

% do some kind of thresholding on the data
range = max(max(sqrt(sum(binnedData.^2,3))));
mask = (sqrt(sum(binnedData.^2,3)) >= fitParm.ThreshVariation/100*range);

% do a test run to find the number of fit parameters and range
curve = 0*squeeze(binnedData(1, 1, :));
[results, fitRange, ~] = fitExpCurve(delays, curve', fitParm);
numFitParm = length(results.Parameters);
fitResults = zeros(newSize(1), newSize(2), numFitParm);
fitErrors = zeros(newSize(1), newSize(2), numFitParm);
fitCurves = zeros(newSize(1), newSize(2), length(fitRange));

%fit the data
if bar wb = waitbar(0,'Fitting data...'); end
for i1 = 1:newSize(1)
    for i2 = 1:newSize(2)
        if bar waitbar( ((i1-1)*newSize(2)+i2)/(newSize(1)*newSize(2)) ); end
        if mask(i1,i2)
            curve = squeeze(binnedData(i1, i2, :));
            [results, ~, trace] = fitExpCurve(delays, curve', fitParm);
            fitResults(i1, i2, :) = results.Values(:,1);
            fitErrors(i1, i2, :) = results.Values(:,2);
            fitCurves(i1, i2, :) = trace;
        end
    end % i2
end % i1
if bar close(wb); end

valueKey = results.Parameters;


% --- bin matrix into larger bin (averages over bin x bin pixels)
function Out = bin3DSet (In, binFactor)
inDim = size(In);
binnedDim = [floor(inDim(1)/binFactor),floor(inDim(2)/binFactor),inDim(3)];
subMatrix = In(1:binnedDim(1)*binFactor,1:binnedDim(2)*binFactor,:);
M = sum(reshape(subMatrix, binFactor, [], inDim(3)),1);
M = reshape(M,binnedDim(1),[], inDim(3)); 
M = permute(M,[2 1 3]);
M = sum(reshape(M,binFactor,[],inDim(3)),1);
M = reshape(M,binnedDim(1),binnedDim(2),inDim(3));
Out = permute(M,[2 1 3])/(binFactor^2);


% --- expand matrix to larger dimension
function Out = expand2DSet (In, expandFactor)
inDim = size(In);
M = In';
M = reshape(M,1,[]);
M = repmat(M,expandFactor,1);
M = reshape(M,expandFactor*inDim(2),[]);
M = M';
M = reshape(M,1,[]);
M = repmat(M,expandFactor,1);
Out = reshape(M,expandFactor*inDim(1),[]);


% --- single curve fit
function [fitResults, fitRange, fitCurve] = fitExpCurve(xData, yData, fitParm)

% select fit range
mask = (xData >= fitParm.InputVal.FitStart);
x = xData(mask);
y = yData(mask);

% create a short name, otherwise the function def is really long (and unreadable)
sig = fitParm.InputVal.CrossCorrWidth;

if (strcmp(fitParm.FitModel,'Single exp'))
    fitResults.Parameters = {'A1', 'T1'};
    LB = [-Inf, 0];
    UB = [Inf, Inf];
    InitCond = [fitParm.GuessVal.A1, fitParm.GuessVal.T1];
    fn = @(A,x) A(1)*exp(-x./A(2));
elseif (strcmp(fitParm.FitModel,'Offset + single exp'))
    fitResults.Parameters = {'Offset', 'A1', 'T1'};
    LB = [-Inf, -Inf, 0];
    UB = [Inf, Inf, Inf];
    InitCond = [fitParm.GuessVal.Offset, fitParm.GuessVal.A1, fitParm.GuessVal.T1];
    fn = @(A,x) A(1)+A(2)*exp(-x./A(3));
elseif (strcmp(fitParm.FitModel,'Double exp'))
    fitResults.Parameters = {'A1', 'T1', 'A2', 'T2'};
    LB = [-Inf, 0, -Inf, 0];
    UB = [Inf, Inf, Inf, Inf];
    InitCond = [fitParm.GuessVal.A1, fitParm.GuessVal.T1, fitParm.GuessVal.A2, fitParm.GuessVal.T2];
    fn = @(A,x) A(1)*exp(-x./A(2))+A(3)*exp(-x./A(4));
elseif (strcmp(fitParm.FitModel,'Offset + double exp'))
    fitResults.Parameters = {'Offset', 'A1', 'T1', 'A2', 'T2'};
    LB = [-Inf, -Inf, 0, -Inf, 0];
    UB = [Inf, Inf, Inf, Inf, Inf];
    InitCond = [fitParm.GuessVal.Offset, fitParm.GuessVal.A1, fitParm.GuessVal.T1,...
        fitParm.GuessVal.A2, fitParm.GuessVal.T2];
    fn = @(A,x) A(1)+A(2)*exp(-x./A(3))+A(4)*exp(-x./A(5));
elseif (strcmp(fitParm.FitModel,'Inst + single exp'))
    fitResults.Parameters = {'Inst','A1','T1'};
    LB = [-Inf,-Inf,0];
    UB = [Inf,Inf,Inf];
    InitCond = [fitParm.GuessVal.Inst, fitParm.GuessVal.A1, fitParm.GuessVal.T1];
    fn = @(A,x) A(1)/(sqrt(pi)*sig)*exp(-x.^2/sig^2)...
        +A(2)*exp(-(4*x*A(3)+sig^2)/(4*A(3)^2)).*(1+erf(x/sig-sig/2/A(3)));
elseif (strcmp(fitParm.FitModel,'Inst + double exp'))
    fitResults.Parameters = {'Inst','A1','T1','A2','T2'};
    LB = [-Inf,-Inf,0,-Inf,0];
    UB = [Inf,Inf,Inf,Inf,Inf];
    InitCond = [fitParm.GuessVal.Inst, fitParm.GuessVal.A1, fitParm.GuessVal.T1,...
        fitParm.GuessVal.A2, fitParm.GuessVal.T2];
    fn = @(A,x) A(1)/(sqrt(pi)*sig)*exp(-x.^2/sig^2)...
        +A(2)*exp(-(4*x*A(3)+sig^2)/(4*A(3)^2)).*(1+erf(x/sig-sig/2/A(3)))...
        +A(4)*exp(-(4*x*A(5)+sig^2)/(4*A(5)^2)).*(1+erf(x/sig-sig/2/A(5)));
elseif (strcmp(fitParm.FitModel,'Offset + inst + single exp'))
    fitResults.Parameters = {'Offset','Inst','A1','T1'};
    LB = [-Inf,-Inf,-Inf,0];
    UB = [Inf,Inf,Inf,Inf];
    InitCond = [fitParm.GuessVal.Offset, fitParm.GuessVal.Inst,...
        fitParm.GuessVal.A1, fitParm.GuessVal.T1];
    fn = @(A,x) A(1) + A(2)/(sqrt(pi)*sig)*exp(-x.^2/sig^2)...
        +A(3)*exp(-(4*x*A(4)+sig^2)/(4*A(4)^2)).*(1+erf(x/sig-sig/2/A(4)));
elseif (strcmp(fitParm.FitModel,'Offset + inst + double exp'))
    fitResults.Parameters = {'Offset','Inst','A1','T1','A2','T2'};
    LB = [-Inf,-Inf,-Inf,0,-Inf,0];
    UB = [Inf,Inf,Inf,Inf,Inf,Inf];
    InitCond = [fitParm.GuessVal.Offset, fitParm.GuessVal.Inst,...
        fitParm.GuessVal.A1, fitParm.GuessVal.T1,...
        fitParm.GuessVal.A2, fitParm.GuessVal.T2];
    fn = @(A,x) A(1) + A(2)/(sqrt(pi)*sig)*exp(-x.^2/sig^2)...
        +A(3)*exp(-(4*x*A(4)+sig^2)/(4*A(4)^2)).*(1+erf(x/sig-sig/2/A(4)))...
        +A(5)*exp(-(4*x*A(6)+sig^2)/(4*A(6)^2)).*(1+erf(x/sig-sig/2/A(6)));
else
    error('model not implemented');
end

%fit options
%options = optimset('MaxIter',12000,'TolFun',1.0e-9,'MaxFunEvals',8000);
options = optimset('MaxIter',1200,'TolFun',1.0e-7,'MaxFunEvals',800,...
            'Display','Off');

% dummy fit to bypass fitting and just return the fit parameter types
% and fit range
if sum(yData.^2)==0 
    fitRange = x;
    fitResults.Values = [zeros(length(fitResults.Parameters),1)...
        zeros(length(fitResults.Parameters),1)];
    fitCurve = yData;
    return
end

% fitting
% 
startPoint = InitCond;

[AA,~,Residual,~,~,~,J] = lsqcurvefit(fn,startPoint,x,y,LB,UB,options);

% [AA,~,Residual,~,~,~,J] = lsqcurvefit(fn,InitCond,x,y,LB,UB,options);
ci = nlparci(AA,Residual,'jacobian',J);%95% confidence intervals
AAerror = (ci(:,2)-ci(:,1))/(2*1.96); %standard error



fitRange = x;
fitResults.Values = [AA' AAerror];
fitCurve = fn(AA,x);


% --- Executes on button press in pb_save.
function pb_save_Callback(hObject, eventdata, handles)
% hObject    handle to pb_save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

binnedData = getappdata(gcbf,'binnedData');
fitParm = getappdata(gcbf,'fitParm');
fitMask = getappdata(gcbf, 'fitMask');
fitResults = getappdata(gcbf,'fitResults');
errors = getappdata(gcbf,'fitErrors');
fitCurves = getappdata(gcbf,'fitCurves');
fitRange = getappdata(gcbf,'fitRange');
names = getappdata(gcbf,'fitValueKey');

pathName = getpref('puprisa','imageWorkingDirectory');

fileName = getappdata(gcbf,'fileName');
saveName = regexprep(fileName, '.dat', '_fit.mat')

[fileName, pathName] = uiputfile([pathName,filesep,saveName],'Save file name');

save([pathName,filesep,fileName], 'binnedData', 'fitParm', 'fitMask', 'fitResults',...
     'errors', 'fitCurves', 'fitRange', 'names');



% --- Executes on selection change in pm_fitModel.
function pm_fitModel_Callback(hObject, eventdata, handles)
% hObject    handle to pm_fitModel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

fitParm = getappdata(gcbf, 'fitParm');
contents = cellstr(get(hObject,'String'));
fitParm.FitModel = contents{get(hObject,'Value')};
setappdata(gcbf, 'fitParm', fitParm);


% --- Executes on selection change in pm_initialGuess.
function pm_initialGuess_Callback(hObject, eventdata, handles)
% hObject    handle to pm_initialGuess (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

fitParm = getappdata(gcbf, 'fitParm');
contents = cellstr(get(hObject,'String'));
fitParm.GuessMethod = contents{get(hObject,'Value')};
setappdata(gcbf, 'fitParm', fitParm);


function setGuessParms_Callback(hObject, eventdata, handles)
% hObject    handle to ed_guessInst (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

fitParm = getappdata(gcbf, 'fitParm');
fitParm.GuessVal.Offset = str2double(get(handles.ed_guessOffset,'String'));
fitParm.GuessVal.Inst = str2double(get(handles.ed_guessInst,'String'));
fitParm.GuessVal.A1 = str2double(get(handles.ed_guessA1,'String'));
fitParm.GuessVal.A2 = str2double(get(handles.ed_guessA2,'String'));
fitParm.GuessVal.A3 = str2double(get(handles.ed_guessA3,'String'));
fitParm.GuessVal.T1 = str2double(get(handles.ed_guessT1,'String'));
fitParm.GuessVal.T2 = str2double(get(handles.ed_guessT2,'String'));
fitParm.GuessVal.T3 = str2double(get(handles.ed_guessT3,'String'));
setappdata(gcbf, 'fitParm', fitParm);


function ed_crossCorrWidth_Callback(hObject, eventdata, handles)
% hObject    handle to ed_crossCorrWidth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

fitParm = getappdata(gcbf, 'fitParm');
fitParm.InputVal.CrossCorrWidth = str2double(get(hObject,'String'));
setappdata(gcbf, 'fitParm', fitParm);


function ed_delayOffset_Callback(hObject, eventdata, handles)
% hObject    handle to ed_delayOffset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

fitParm = getappdata(gcbf, 'fitParm');
fitParm.InputVal.DelayOffset = str2double(get(hObject,'String'));
setappdata(gcbf, 'fitParm', fitParm);


function ed_fitStart_Callback(hObject, eventdata, handles)
% hObject    handle to ed_fitStart (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ed_fitStart as text
%        str2double(get(hObject,'String')) returns contents of ed_fitStart as a double

fitParm = getappdata(gcbf, 'fitParm');
fitParm.InputVal.FitStart = str2double(get(hObject,'String'));
setappdata(gcbf, 'fitParm', fitParm);


function ed_binSize_Callback(hObject, eventdata, handles)
% hObject    handle to ed_binSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

fitParm = getappdata(gcbf, 'fitParm');
fitParm.BinSize = str2double(get(hObject,'String'));
setappdata(gcbf, 'fitParm', fitParm);


function ed_threshVar_Callback(hObject, eventdata, handles)
% hObject    handle to ed_threshVar (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

fitParm = getappdata(gcbf, 'fitParm');
fitParm.ThreshVariation = str2double(get(hObject,'String'));
setappdata(gcbf, 'fitParm', fitParm);


function ed_changeLimFitDisp_Callback(hObject, eventdata, handles)
% hObject    handle to ed_minValDisp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

min = str2double(get(handles.ed_minValDisp,'String'));
setappdata(gcbf,'minValDisplay', min);

max = str2double(get(handles.ed_maxValDisp,'String'));
setappdata(gcbf,'maxValDisplay', max);

pm_fitDisplay_Callback(handles.pm_fitDisplay, eventdata, handles);



% --- other stuff, like various create functions

% --- Executes on button press in cb_intensityCoded.
function cb_intensityCoded_Callback(hObject, eventdata, handles)
% hObject    handle to cb_intensityCoded (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

pm_fitDisplay_Callback(handles.pm_fitDisplay, eventdata, handles);

% --- Executes on slider movement.
function sl_contrastEnhance_Callback(hObject, eventdata, handles)
% hObject    handle to sl_contrastEnhance (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

pm_fitDisplay_Callback(handles.pm_fitDisplay, eventdata, handles);

% --- Executes during object creation, after setting all properties.
function sl_delayIndex_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sl_delayIndex (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes during object creation, after setting all properties.
function ed_crossCorrWidth_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ed_crossCorrWidth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes during object creation, after setting all properties.
function ed_delayOffset_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ed_delayOffset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes during object creation, after setting all properties.
function pm_initialGuess_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pm_initialGuess (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function ed_guessA1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ed_guessA1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function ed_guessT1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ed_guessT1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function ed_guessA2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ed_guessA2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function ed_guessT2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ed_guessT2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function ed_guessA3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ed_guessA3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function ed_guessT3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ed_guessT3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function ed_guessInst_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ed_guessInst (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function ed_guessOffset_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ed_guessOffset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function pm_fitModel_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pm_fitModel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function ed_fitStart_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ed_fitStart (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function ed_binSize_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ed_binSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes during object creation, after setting all properties.
function ed_threshVar_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ed_threshVar (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes during object creation, after setting all properties.
function pm_fitDisplay_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pm_fitDisplay (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes during object creation, after setting all properties.
function ed_minValDisp_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ed_minValDisp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function ed_maxValDisp_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ed_maxValDisp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function sl_contrastEnhance_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sl_contrastEnhance (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


