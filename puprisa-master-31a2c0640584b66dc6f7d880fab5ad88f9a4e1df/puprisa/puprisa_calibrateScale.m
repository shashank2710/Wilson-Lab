function varargout = puprisa_calibrateScale(varargin)
% PUPRISA_CALIBRATESCALE M-file for puprisa_calibrateScale.fig
%      PUPRISA_CALIBRATESCALE, by itself, creates a new PUPRISA_CALIBRATESCALE or raises the existing
%      singleton*.
%
%      H = PUPRISA_CALIBRATESCALE returns the handle to a new PUPRISA_CALIBRATESCALE or the handle to
%      the existing singleton*.
%
%      PUPRISA_CALIBRATESCALE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PUPRISA_CALIBRATESCALE.M with the given input arguments.
%
%      PUPRISA_CALIBRATESCALE('Property','Value',...) creates a new PUPRISA_CALIBRATESCALE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before puprisa_calibrateScale_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to puprisa_calibrateScale_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
%   may optionally pass in a list of calibration profile names,
%   so that when user clicks 'done', we can double-check that the requested
%   profile name won't overrwrite an existing profile.
%
%   returns x and y scale in volts per micron
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help puprisa_calibrateScale

% Last Modified by GUIDE v2.5 13-Mar-2012 15:43:38

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @puprisa_calibrateScale_OpeningFcn, ...
                   'gui_OutputFcn',  @puprisa_calibrateScale_OutputFcn, ...
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


% --- Executes just before puprisa_calibrateScale is made visible.
function puprisa_calibrateScale_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to puprisa_calibrateScale (see VARARGIN)

% Choose default command line output for puprisa_calibrateScale
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% set colormap
colormap(gray);

if nargin == 1
    % passed in a cell array of saved calibration profiles
    setappdata(gcbf,'profileNames',varargin{1});
end

% UIWAIT makes puprisa_calibrateScale wait for user response (see UIRESUME)
 uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function [x,y] = puprisa_calibrateScale_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
if ( handles.output == 1 )
    x = getappdata(hObject,'scaleX_V_per_micron');
    y = getappdata(hObject,'scaleY_V_per_micron');
else
    x = [];
    y = [];
end

delete(hObject);

% --- Executes on button press in pbLoadImage.
%       Loads the first image
function pbLoadImage_Callback(hObject, eventdata, handles)
% hObject    handle to pbLoadImage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    loadImage( 1 );


% --- Executes on button press in pbLoadShiftedImage.
function pbLoadShiftedImage_Callback(hObject, eventdata, handles)
% hObject    handle to pbLoadShiftedImage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    loadImage( 2 );


function edShiftX_micron_Callback(hObject, eventdata, handles)
% hObject    handle to edShiftX_micron (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edShiftX_micron as text
%        str2double(get(hObject,'String')) returns contents of edShiftX_micron as a double


% --- Executes during object creation, after setting all properties.
function edShiftX_micron_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edShiftX_micron (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pbDone.
function pbDone_Callback(hObject, eventdata, handles)
% hObject    handle to pbDone (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    handles.output = 1;
    guidata(gcbf, handles);
    uiresume(gcbf);


% --- Executes on button press in pbCancel.
function pbCancel_Callback(hObject, eventdata, handles)
% hObject    handle to pbCancel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    handles.output = 0;
    guidata(gcbf, handles);
    uiresume(gcbf);


function edShiftY_micron_Callback(hObject, eventdata, handles)
% hObject    handle to edShiftY_micron (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edShiftY_micron as text
%        str2double(get(hObject,'String')) returns contents of edShiftY_micron as a double


% --- Executes during object creation, after setting all properties.
function edShiftY_micron_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edShiftY_micron (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edScanRangeX_Callback(hObject, eventdata, handles)
% hObject    handle to edScanRangeX (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edScanRangeX as text
%        str2double(get(hObject,'String')) returns contents of edScanRangeX as a double


% --- Executes during object creation, after setting all properties.
function edScanRangeX_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edScanRangeX (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edScanRangeY_Callback(hObject, eventdata, handles)
% hObject    handle to edScanRangeY (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edScanRangeY as text
%        str2double(get(hObject,'String')) returns contents of edScanRangeY as a double


% --- Executes during object creation, after setting all properties.
function edScanRangeY_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edScanRangeY (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pbCalcShift.
%   Calculates shift between two images
function pbCalcShift_Callback(hObject, eventdata, handles)
% hObject    handle to pbCalcShift (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    
    % get image data
    im1 = get( findobj(handles.axImage1,'Type','Image'), 'CData' );
    im2 = get( findobj(handles.axImage2,'Type','Image'), 'CData' );
    
    [nx, ny] = size(im1);
    
    if( get(handles.cbHighPass, 'value') == 1 )
        % high pass filter
        im1 = ifilt_highPass( im1 );
        im2 = ifilt_highPass( im2 );
    end
    
    % get shift amount
    shiftX_micron = str2num( get( handles.edShiftX_micron, 'String' ) );
    shiftY_micron = str2num( get( handles.edShiftY_micron, 'String' ) );
   
    % check for valid numbers
    if( isempty(shiftX_micron) || isempty(shiftY_micron ) )
        errordlg('Invalid entry in measured x/y shift');
    else


        % calculate shift
        [output, shifted] = dftregistration(...
                fft2(im1), fft2(im2),100);

        diffphase = output(2);
        row_shift = output(3);
        col_shift = output(4);

        scaleX_px_per_micron = abs(col_shift) / shiftX_micron
        scaleY_px_per_micron = abs(row_shift) / shiftY_micron
        
        scanRangeX = str2num(get(handles.edScanRangeX,'String'));
        scanRangeY = str2num(get(handles.edScanRangeY,'String'));
        
        scaleX_V_per_micron = scaleX_px_per_micron * scanRangeX / nx;
        scaleY_V_per_micron = scaleY_px_per_micron * scanRangeY / ny;
        
        setappdata(gcbf,'scaleX_V_per_micron',scaleX_V_per_micron);
        setappdata(gcbf,'scaleY_V_per_micron',scaleY_V_per_micron);
        
        
        handles.output(1) = scaleX_V_per_micron;
        handles.output(2) = scaleY_V_per_micron;
        
        guidata(gcbf, handles);
        
        % show the results
        set(handles.edScaleX, 'String', num2str(scaleX_V_per_micron));
        set(handles.edScaleY, 'String', num2str(scaleY_V_per_micron));
        
        % shift the original images and overlay
        % restore image data
        im1 = get( findobj(handles.axImage1,'Type','Image'), 'CData' );
        im2 = get( findobj(handles.axImage2,'Type','Image'), 'CData' );
        displayShiftedOverlay( im1, im2, diffphase, row_shift, col_shift);
    end



function edScaleX_Callback(hObject, eventdata, handles)
% hObject    handle to edScaleX (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edScaleX as text
%        str2double(get(hObject,'String')) returns contents of edScaleX as a double


% --- Executes during object creation, after setting all properties.
function edScaleX_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edScaleX (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edScaleY_Callback(hObject, eventdata, handles)
% hObject    handle to edScaleY (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edScaleY as text
%        str2double(get(hObject,'String')) returns contents of edScaleY as a double


% --- Executes during object creation, after setting all properties.
function edScaleY_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edScaleY (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%==========================================================================
%
% Utility functions

% --- Loads image into specified axis
% whichimage = 1 for original, whichImage = 2 for shifted iamge
function loadImage( whichImage )
    % first prompt user for filename
    
    % recall working directory
    if ispref('puprisa','calibrateScaleWorkingDirectory')
        wd = getpref('puprisa','calibrateScaleWorkingDirectory');
    else
        % or use MATLAB's working directory if none has been set yet
        wd = pwd();
    end
    
    
    % or use MATLAB's working directory if none has been set yet
    if isempty(wd)
        wd = pwd();
    end
    
    [fileName, pathName] = uigetfile([wd,'/*.dat'],...
        'Select an ImageScanner Stack' );
    
    if ~isequal(fileName,0)
        % remember the working directory for later
        setpref('puprisa','calibrateScaleWorkingDirectory',pathName);

        
        
        % load the image stack
        [slices, header] = puprisa_readImageStack( [pathName,'/',fileName] );
        
        % store loaded data
        slicesArray = getappdata(gcbf,'slicesArray');
        if isempty(slicesArray)
            clear slicesArray;
        end
        headerArray = getappdata(gcbf,'headerArray');
        if isempty( headerArray )
            clear headerArray;
        end
        
        slicesArray(whichImage) = slices;
        headerArray(whichImage) = header;
        
        setappdata(gcbf,'slicesArray',slicesArray);
        setappdata(gcbf,'headerArray',headerArray);
        setappdata(gcbf,'scanRangeX', header.scanrangex);
        setappdata(gcbf,'scanRangeY', header.scanrangey);
        
        [nx, ny] = size(slices(1).imageData{1});
        setappdata(gcbf,'nPxX', nx);
        setappdata(gcbf,'nPxY', ny);
        
        
        updateImageDisplay();

        % get voltage information, etc
        setScanRange( header.scanrangex,header.scanrangey );
        
    end

% --- Set X and Y scan range, in Volts
function setScanRange( scanRangeX, scanRangeY )
    h = guihandles(gcbf);

    set(h.edScanRangeX,'String', num2str(scanRangeX));
    set(h.edScanRangeY,'String', num2str(scanRangeX));
    
    
function updateImageDisplay()
    h = guihandles(gcbf);
    chan = str2num(get(h.edChannel,'String'));
    
    if isempty(chan)
        errordlg('Invalid channel');
    else
        % get loaded data
        slicesArray = getappdata(gcbf,'slicesArray');
        headerArray = getappdata(gcbf,'headerArray');
        
        if length(slicesArray) > 0
        % show first image
        slices = slicesArray(1);
        
        imageData = slices(1).imageData{chan};
        
        % and display it
        axes(h.axImage1);
        imagesc( imageData,'parent',h.axImage1 );
        end
        
        if length(slicesArray) > 1
        % show second image
        slices = slicesArray(2);
        imageData = slices(1).imageData{chan};
        
        % and display it
        axes(h.axImage2);
        imagesc( imageData );
        end

    end

function edChannel_Callback(hObject, eventdata, handles)
% hObject    handle to edChannel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edChannel as text
%        str2double(get(hObject,'String')) returns contents of edChannel as a double


% --- Executes during object creation, after setting all properties.
function edChannel_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edChannel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    handles.output = 0;
    guidata(gcbf,handles);
    uiresume(gcbf);


% --- Executes on button press in pbHiPass.
function pbHiPass_Callback(hObject, eventdata, handles)
% hObject    handle to pbHiPass (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in cbHighPass.
function cbHighPass_Callback(hObject, eventdata, handles)
% hObject    handle to cbHighPass (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of cbHighPass



function edProfileName_Callback(hObject, eventdata, handles)
% hObject    handle to edProfileName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edProfileName as text
%        str2double(get(hObject,'String')) returns contents of edProfileName as a double


% --- Executes during object creation, after setting all properties.
function edProfileName_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edProfileName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- display shifted images in an overlay
function displayShiftedOverlay( im1, im2, diffphase, row_shift, col_shift)
    h = guihandles( gcbf );
    im2Shifted = dftapplyshift( im2, diffphase, row_shift, col_shift );
    
    cscale1 = get( h.axImage1, 'CLim' );
    cscale2 = get( h.axImage2, 'CLim' );
    
    X = stretchValues( im1, cscale1, [0, 1], 1 );
    Y = stretchValues( im2Shifted, cscale2, [0, 1], 1 );
    
    [szx, szy] = size(im1);
    
    % assume im2 is same size
    
    imRGB = zeros(szx, szy, 3);
    
    imRGB(:,:,1) = X;
    imRGB(:,:,2) = Y;
    
    axes(h.axAligned);
    imshow(imRGB, 'Parent', h.axAligned);
    set(gca,'ydir','normal');
    
    % draw a scale bar
    scanRangeX = getappdata(gcbf,'scanRangeX');
    scanRangeY = getappdata(gcbf,'scanRangeY');
        
    
    nx = getappdata(gcbf,'nPxX');
    ny = getappdata(gcbf,'nPxY');
    
    scaleX_V_per_micron = getappdata(gcbf,'scaleX_V_per_micron');
    scaleY_V_per_micron = getappdata(gcbf,'scaleY_V_per_micron');
    
    scaleX_px_per_micron = scaleX_V_per_micron * nx / scanRangeX;

    hsb = puprisa_scaleBar( 100, scaleX_px_per_micron);
    drawnow;



function edName_Callback(hObject, eventdata, handles)
% hObject    handle to edName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edName as text
%        str2double(get(hObject,'String')) returns contents of edName as a double


% --- Executes during object creation, after setting all properties.
function edName_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pbSave.
function pbSave_Callback(hObject, eventdata, handles)
% hObject    handle to pbSave (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

name = get(handles.edName,'string');

% create new calibration structure, if necessary
if ~ispref('puprisa','scaleCalibrationSet')
    addpref('puprisa','scaleCalibrationSet',[]);
end

s = getpref('puprisa','scaleCalibrationSet');
ns = length(s);

this.name = name;
this.scaleX_V_per_micron = getappdata(gcbf,'scaleX_V_per_micron');
this.scaleY_V_per_micron = getappdata(gcbf,'scaleY_V_per_micron');

if ns == 0
    s = this;
else
    s(ns+1) = this;
end

setpref('puprisa','scaleCalibrationSet',s);