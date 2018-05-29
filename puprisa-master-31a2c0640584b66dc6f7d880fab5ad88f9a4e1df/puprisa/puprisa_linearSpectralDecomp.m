function varargout = puprisa_linearSpectralDecomp(varargin)
% PUPRISA_LINEARSPECTRALDECOMP M-file for puprisa_linearSpectralDecomp.fig
%      PUPRISA_LINEARSPECTRALDECOMP, by itself, creates a new PUPRISA_LINEARSPECTRALDECOMP or raises the existing
%      singleton*.
%
%      H = PUPRISA_LINEARSPECTRALDECOMP returns the handle to a new PUPRISA_LINEARSPECTRALDECOMP or the handle to
%      the existing singleton*.
%
%      PUPRISA_LINEARSPECTRALDECOMP('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PUPRISA_LINEARSPECTRALDECOMP.M with the given input arguments.
%
%      PUPRISA_LINEARSPECTRALDECOMP('Property','Value',...) creates a new PUPRISA_LINEARSPECTRALDECOMP or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before puprisa_linearSpectralDecomp_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to puprisa_linearSpectralDecomp_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES
%
%
% INTERNAL DATA STRUCTURES (appdata)
%
% baselineChannels  Indices to entries in the standards list that indicate
%                   baseline/offset correction channels.
%
% baselineMethod    Setting for the method used for baseline (offset)
%                   correction method. This is a string, and can be one of
%                   the following:
%
%                     'BASELINE_REDUCEDRANK'
%                       
%                       Calculates baseline by making a reduced-rank
%                       approximation of the image stack, and uses the
%                       first image slice as a baseline for each pixel (we 
%                       assume this is the most-negative time delay, i.e. 
%                       the probe precedes the pump in this slice)
%
%                     'BASELINE_INDEPENDENTCHANNEL'
%
%                       Rather than correcting for baseline in the
%                       preprocessing steps, add a baseline channel to the
%                       reference standards.
%
%                     'BASELINE_FRAMEAVERAGE'
%
%                       Use the average of the slice with most-negative
%                       time delay as the baseline of all pixels.
%
%                     'BASELINE_NONE'
%
%                       Don't do any baseline corrections.
%
%                   baselineMethod is remembered in preferences as
%                   'lindecomp_baselineMethod'.
%
% delays            A vector of probe delays, each element corresponding to 
%                   a slice of imageStack: The slice at imageStack(:,:,k)
%                   was recorded at probe delay delays(k).
%
% delayIndex        A vector of indices, used to keep track of which probe
%                   delays correspond to individual, appended image stacks.
%                   e.g. for two delay stacks with 5 delays sampled each,
%                   delayIndex = [1 1 1 1 1 2 2 2 2 2]; for a single delay
%                   stack with 6 delays sampled, delayIndex = [1 1 1 1 1 1]
%
% enablePhasor      Enable phasor segmentation. Setting this to 1 causes
%                   phasors to be calculated after running linear unmixing.
%
% fileNames         A cell array of filenames corresponding to the image
%                   stacks that have been loaded and appended together for 
%                   processing.
%
% imageStack        The delay stack, a 3D array of doubles, size is
%                   [nr,nc,nt], where nr is no. rows in the image, nc is
%                   no. cols in the image, and nt is no. of probe delays
%                   sampled. For images acquired with multiple wavelength
%                   combinations, the size of the appended image stack is
%                   [nr,nc, nt1+nt2+...], where nt1, nt2, etc., are the no.
%                   of probe delays for each measurement at a different
%                   wavelength combination.
%
% imageStack_original   A backup copy of imageStack, so we can undo changes
%                       (such as reduced-rank approximation, offset
%                       correction, etc).
%
% negValRejectMethod    Method used for rejecting negative values after
%                       unmixing. Chemical concentrations, obviously,
%                       cannot be negative. The question is whether to
%                       simply reject only the negative information on a 
%                       per-channel basis (e.g. a pixel with 0.5 eumelanin
%                       content and -0.2 pheomelanin will be transformed to
%                       simply having 0.5 eumelanin and 0 pheomelanin) or
%                       on a whole-pixel basis (e.g. the previous example
%                       would be transformed to zero content for all
%                       channels). Rejecting on a whole-pixel basis tends
%                       to reject pixels that have extra components that
%                       are poorly accounted for by the model under
%                       consideration. Options are 
%                       'NEGREJECT_INDEPCHANNEL', 'NEGREJECT_WHOLEPX', or
%                       'NEGREJECT_NONE'.
%
% nStacksLoaded     the number of delay stacks loaded. Normally this is 1,
%                   unless multiple stacks have been appended (as in when
%                   analyzing data acquired with multiple pump/probe
%                   wavelength combinations)
%
% projThreshEnable  0 or 1, to disable or enable projection thresholding
%                   after unmixing
%
% projThreshVal     threshold value for rejecting pixels that have very low
%                   identifiable composition
%
% reduceRankEnable  0 or 1, to disable or enable rank reduction in the
%                   preprocessing steps (used to de-noise the image stack)
%
% reduceRankK       integer >= 1, the number of principal components to 
%                   retain in rank reduction. Typical value of 6 should be
%                   fine for most multiexponential signals.
%
% residualThresEnable   0 or 1, to disalbe/enable discarding pixels having
%                       pump-probe responses that are not well-described by
%                       the model
%
% residualThreshVal     RMS error threshold for residual thresholding
%
% standards         An array of pump-probe standards. Each entry is a
%                   struct with the following fields:
%
%   standards.name      Name used to label the standard
%
%   standards.color     RGB color to display the standard
%
%   standards.scale     scaling factor
%
%   standards.fileName  .mat file in the standards directory that contains
%                       the pump-probe response curve for this standard
%
%   standards.t         probe delays
%
%   standards.x         pump-probe response
%
%   standards.wlpairs   (optional, for multi-wavelength) a 2D array of the
%                       pump and probe wavelengths used. Each row contains
%                       a pump/probe wavelength pair, in nm. 1st column
%                       lists pump wavlenegths; 2nd column list probe
%                       wavelengths.
%
%   standards.tindex    (optional, for multi-wavelength) similar to 
%                       delayIndex, but applied to the individual standards
%                       rather than the delay stack; see description of 
%                       delayIndex above for more info.
%
%   standards.lineobj   handle to line object that shows the standard
%
% standards_original copy of standards, to allow changes to be undone
%
% standardsLoaded   0 or 1, depending on whether valid standards have been
%                   loaded. 
%
% standardsNormalized 0 or 1, depending on whether standards have been
%                     normalized (as unit vectors)
%
% unmixedCoeffs     After decomposition is successfully run, this variable
%                   stores the mixing coefficients for each pixel. It is a
%                   3D matrix of size [nr, nc, nch], where nr, nc are the
%                   original image's number of rows and columns, and nch is
%                   the number of channels (number of standards that were 
%                   used to generate the orthonormal basis).
%
% This program written by Jesse Wilson (jw295@duke.edu) 2012 for data 
% analysis in the lab of Warren S. Warren. Contributions by Tana E. 
% Villafana, Mary Jane Simpson, and Francisco E. Robles.

% Edit the above text to modify the response to help puprisa_linearSpectralDecomp

% Last Modified by GUIDE v2.5 29-Sep-2013 11:02:34

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @puprisa_linearSpectralDecomp_OpeningFcn, ...
                   'gui_OutputFcn',  @puprisa_linearSpectralDecomp_OutputFcn, ...
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


% --- Executes just before puprisa_linearSpectralDecomp is made visible.
function puprisa_linearSpectralDecomp_OpeningFcn(hObject, eventdata, handles, varargin)
    % This function has no output args, see OutputFcn.
    % hObject    handle to figure
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    % varargin   command line arguments to puprisa_linearSpectralDecomp (see VARARGIN)
    %            passing in a delay stack and a time vector will load that
    %            stack initially.

    % Choose default command line output for puprisa_linearSpectralDecomp
    handles.output = hObject;

    set(hObject,'toolbar','figure');
    setappdata(hObject,'standardsLoaded',0);

    currentfolder = pwd;
    [filepath,~, ~] =fileparts( which('puprisa_linearSpectralDecomp'));

    switch currentfolder 
        case filepath
        otherwise
            cd(filepath)
    end

    % Update handles structure
    guidata(hObject, handles);

    batchMode = 0;
    if nargin > 4
        if isstr( varargin{ nargin - 3 } )
            if strcmp( varargin{ nargin -3 }, 'batch' )
                batchMode = 1;
            end
        end
    end
    setappdata(hObject,'batchMode',batchMode);

    % load an image, if it was passed in the input arguments
    if nargin >= 5
        setappdata(hObject,'imageStack',varargin{1});
        setappdata(hObject,'imageStack_original',varargin{1});
        setappdata(hObject,'delays',varargin{2});
    end
    
    % restore preferences
    restorePreferences( hObject, handles );
    

    % UIWAIT makes puprisa_linearSpectralDecomp wait for user response (see UIRESUME)
    % uiwait(handles.figure1);

    handles.output = 2;
    guidata(hObject,handles);
% end  % puprisa_linearSpectralDecomp_OpeningFcn()



function restorePreferences( hObject, handles )
    
    % recall baseline correction settings
    restoreEnumPreference( hObject, 'lindecomp_baselineMethod', ...
        'baselineMethod', handles.pmnu_baselineCorrection, ...
        {'BASELINE_REDUCEDRANK', ...
         'BASELINE_INDEPENDENTCHANNEL', ...
         'BASELINE_FRAMEAVERAGE', ...
         'BASELINE_NONE'}, 1 ); 

    % recall prefs on normalizing standards
    restorePreference( hObject, 'normalizeStandards', ...
        'standardsNormalized', handles.cb_normalizeStandards, 1 );

    % load most recently used spectral standards
    if ispref('puprisa','standardsSetsDir') && ispref('puprisa','standardsFile')
        pathName = getpref('puprisa','standardsSetsDir');
        standardsFileName = getpref('puprisa','standardsFile');
        if (~isempty( pathName )) && (~isempty( standardsFileName ));
            loadStandards( hObject, [pathName,filesep,standardsFileName]);
        end
    end
    
    % restore settings for reduced-rank approximation
    restorePreference( hObject, 'lindecomp_reduceRankEnable', ...
        'reduceRankEnable', handles.cb_reduceRank, 0 );
    restorePreference( hObject, 'lindecomp_reduceRankK', ...
        'reduceRankK', handles.ed_rank, 6 );
    
    % restore settings for projection thresholding
    restorePreference( hObject, 'lindecomp_projThreshEnable', ...
        'projThreshEnable', handles.ProjThresh, 0 );
    restorePreference( hObject, 'lindecomp_projThreshVal', ...
        'projThreshVal', handles.ThresholdVal, 0.1 );

    % restore settings for phasor analysis
    restorePreference( hObject, 'lindecomp_enablePhasor', ...
        'enablePhasor', handles.cb_enablePhasor, 0 );
    
    % restore settings for residual thresholding
    restorePreference( hObject, 'lindecomp_residualThreshEnable', ...
        'residualThreshEnable', handles.cb_threshold, 0 );
    restorePreference( hObject, 'lindecomp_residualThreshVal', ...
        'residualThreshVal', handles.ed_threshold, 0.01 );
    
    % restore settings for negative value rejection
    restoreEnumPreference( hObject, 'lindecomp_negValRejectMethod', ...
        'negValRejectMethod', handles.NegValIgnore, ...
        {'NEGREJECT_INDEPCHANNEL', 'NEGREJECT_WHOLEPX', ...
         'NEGREJECT_NONE'}, 2 ); 
    
% end % restorePreferences()

% --- called by restorePreferences()
%     Restore an individual preference, given its name, the name it will
%     have in appdata, a handle to the ui control that sets this
%     preference, and a default value
function restorePreference( fig, prefName, appdataName, uihandle, default )
    
    if ispref( 'puprisa', prefName )
        % if this has been set before, restore its previous value
        val = getpref('puprisa', prefName );
    else
        % otherwise, set it to the default
        val = default;
    end
    
    % set the UI element's value accordingly
    if ~isempty( uihandle )
        if ischar( val )
            set( uihandle, 'String', val );
        else
            switch( get(uihandle,'Style') )
                case 'checkbox'
                    set( uihandle, 'Value', val);
                case 'edit'
                    set( uihandle, 'String', num2str(val));
                otherwise
                    error(['UIControl type ', get(uihandle,'Style'), ...
                        ' not yet supported by restorePreference()']);
            end
        end
    end
    
    % set appdata value accordingly
    setappdata( fig, appdataName, val );
    
% end % restorePreference()

% --- called by restorePreferences, to restore a setting that can be one of
%     a specified number of options (enumerated). Similar to
%     restorePreference(), except it has an extra argument, optionsList,
%     which is a cell array of strings. The argument default is a number
%     that specifies which entry in options list is the default. uihandle
%     must be a handle to a popupmenu control
function restoreEnumPreference( fig, prefName, appdataName, uihandle, ...
    optionsList, default )

    if ispref( 'puprisa', prefName )
        val = getpref( 'puprisa', prefName );
    else
        val = optionsList{default};
    end
        
    % find index of the current value in the list of options
    k = find( strcmp( optionsList, val ) );
        
    % if we cannot find the current value in the list of options, then
    % raise an error
    if isempty(k)
        errstr = sprintf(['Invalid preference: ',...
            prefName, '.',...
            'This sometimes happens after upgrading PUPRISA ',...
            'to a new version, if the old version''s ',...
            'settings are incompatible. Try the following ',...
            'command to fix: \n\n',...
            '  rmpref(''puprisa'',''lindecomp_baselineMethod'')']);
        error( errstr );
    end
    
    % save this setting in appdata
    setappdata( fig, appdataName, val );

    % set the popupmenu's current selection to the one specified
    if ~isempty( uihandle )
        
        % double-check that a popupmenu was passed in
        if( strcmp( get( uihandle, 'Style' ), 'popupmenu' ) == 0 )
            error(['restoreEnumPreference() must be given a handle to ',...
                'a popupmenu control when passing in uihandle.']);
        end
        
        % set the value of the popup menu
        set( uihandle, 'Value', k );
        
        % store the options list in this control's appdata
        setappdata( uihandle, 'enumOptionsList', optionsList );
    end
% end % restoreEnumPreference()

% --- called by ui controls, when changing a setting, store it in appdata
%     as well as the puprisa preferences group
function setWithPref( fig, prefName, appdataName, val )
    setappdata( fig, appdataName, val );
    if ispref( 'puprisa', prefName )
        setpref( 'puprisa', prefName, val );
    else
        addpref( 'puprisa', prefName, val );
    end
% end % setWithPref()

% --- called by popupmenu controls, when changing a setting: store it in 
%     appdata as well as the puprisa preferences group
function setWithPrefEnum( fig, prefName, appdataName, uicontrol )
    % figure out which option was selected
    k = get( uicontrol, 'Value' );
    optionsList = getappdata( uicontrol, 'enumOptionsList' );
    
    val = optionsList{k};
    
    setappdata( fig, appdataName, val );
    if ispref( 'puprisa', prefName )
        setpref( 'puprisa', prefName, val );
    else
        addpref( 'puprisa', prefName, val );
    end
    
% end % setWithPrefEnum()

% --- Outputs from this function are returned to the command line.
function varargout = puprisa_linearSpectralDecomp_OutputFcn(hObject, eventdata, handles) 
    % varargout  cell array for returning output args (see VARARGOUT);
    % hObject    handle to figure
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)

    % Get default command line output from handles structure

    if getappdata(hObject,'batchMode' )
        % go ahead and run decomposition, if we're set up to do so
        if getappdata(hObject,'standardsLoaded') == 1
            set(handles.pb_decompose,'Enable','on');
            pb_decompose_Callback(handles.pb_decompose, eventdata, handles);
        end

        varargout{1} = getappdata(gcf,'fractionalComposition');

        ll = findobj(handles.ax_components,'type','line');
        varargout{2} = get(ll,'DisplayName');

        close( hObject );
    end
% end % puprisa_linearSpectralDecomp_OutputFcn()

% --- Executes on button press in pb_loadImageStack.
function pb_loadImageStack_Callback(hObject, eventdata, handles)
    
    % prompt user for an image stack and load the data
    [t,X,fileName] = LoadImageStack();

    setappdata( gcbf, 'fileNames', {fileName} );
    setappdata( gcbf, 'imageStack', X );
    setappdata( gcbf, 'imageStack_original', X );
    setappdata( gcbf, 'delays', t );
    setappdata( gcbf, 'nStacksLoaded', 1 );
    setappdata( gcbf, 'delayIndex', ones(1, length(t)) );
    
    % enable append button
    set( handles.pb_appendStack, 'Enable', 'on' );

    % go ahead and run decomposition, if we're set up to do so
    if getappdata(gcbf,'standardsLoaded') == 1
        set(handles.pb_decompose,'Enable','on');
        
        standards = getappdata( gcbf, 'standards' );
        if isfield( standards(1), 'wlpairs' )
            [nWavelengthPairs,~] = size( standards(1).wlpairs );
        else
            nWavelengthPairs = 1;
        end

        if nWavelengthPairs == 1
            pb_decompose_Callback(hObject, eventdata, handles);
        end
    end
% end % pb_loadImageStack_Callback()
    
% --- Executes on button press in pb_appendStack.
function pb_appendStack_Callback(hObject, eventdata, handles)

    % prompt user for an image stack and load the data
    [t,X,fileName] = LoadImageStack();
    
    % double-check these stacks have the same number of rows and columns
    X0 = getappdata(gcbf, 'imageStack');
    [nr0, nc0, ~] = size(X0);
    [nr, nc, ~] = size(X);
    if ( (nr0 ~= nr) || (nc0 ~= nc) )
        msgbox(['Error: cannot append a stack with a different ', ...
            'number of rows and columns!','Error Appending Stack', ...
            'error']);
        return;
    end

    % update number of loaded stacks
    nStacksLoaded = getappdata(gcbf, 'nStacksLoaded');
    nStacksLoaded = nStacksLoaded + 1;
    setappdata(gcbf, 'nStacksLoaded', nStacksLoaded );

    % append filename to the list
    fileNames = getappdata(gcbf, 'fileNames');
    fileNames{end+1} = fileName;
    setappdata(gcbf, 'fileNames', fileNames);

    % append the image stack itself
    Xa = cat(3, X0, X);
    setappdata(gcbf, 'imageStack', Xa);
    
    % store a copy as the 'original'
    setappdata(gcbf, 'imageStack_original', Xa);
    
    % append time delays
    t0 = getappdata(gcbf, 'delays');
    ta = [t0, t];
    setappdata(gcbf, 'delays', ta);
    
    % append the delay index
    tindex0 = getappdata( gcbf, 'delayIndex');
    tindex = [tindex0, zeros(1, length(t))+nStacksLoaded];
    setappdata( gcbf, 'delayIndex', tindex );
    
    % go ahead and run decomposition, if we're set up to do so
%     if getappdata(gcbf,'standardsLoaded') == 1
%         set(handles.pb_decompose,'Enable','on');
%         pb_decompose_Callback(hObject, eventdata, handles);
%     end
% end % pb_appendStack_Callback()
    
function [t,X, fileName] = LoadImageStack()
% prompt user for image stack data, and load it
% this function is called by the pushbuttons "Load Image Stack" and 
% "Append Stack"
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
    if isequal(fileName, 0); return; end

    setpref('puprisa','imageWorkingDirectory',pathName);


    [S,header] = puprisa_readImageStack( [pathName,filesep,fileName] );
    [t,X] = puprisa_getChannelFromSlices(S,1,header);
% end % LoadImageStack()

% --- Executes on button press in pb_viewOrthBasis.
function pb_viewOrthBasis_Callback(hObject, eventdata, handles)
% perform graham-schmidt orthogonalization, then show results
    [Q,D] = orthonormalize( gcbf, handles );

    figure;
    plot(Q);

    disp('component to orthonormal basis mixing matrix:')
    disp(D);

function [Q,D] = orthonormalize( fig, handles )
% generate an orthonormal basis from the standards we have loaded, after
% interpolating the standards to the same timebase used to acquire the
% image stack

    % retrieve pump-probe standards
    standards = getappdata( gcf, 'standards' );
    nStandards = length( standards );

    % retrieve the timebase used for the image stack
    t = getappdata( gcf, 'delays' );
    nt = length( t );
    delayIndex = getappdata( gcf, 'delayIndex' );
    
    % recall how many pump-probe stacks are loaded
    % (this should be the same as the number of pump/probe wavelength 
    % combinations that were used in total)
    nStacksLoaded = getappdata( gcf, 'nStacksLoaded' );

    % initialize array used to hold interpolated standards
    A = zeros( nt, nStandards );

    % construct a matrix to feed to QR decomp; each component in a column
    for ic = 1:nStandards
        
        % check to see if the standards were produced with the same number
        % of wavelength combinations as the loaded image stacks
        if isfield( standards(ic), 'wlpairs' )
            [nWavelengthPairs,~] = size( standards(ic).wlpairs );
        else
            nWavelengthPairs = 1;
        end
        
        if nWavelengthPairs ~= nStacksLoaded
            errStr = sprintf( ['Error interpolating standards to ',...
                'delay stack timebase. Cannot decompose a set of %i ',...
                'appended image stacks with pump-probe standards that ',...
                'were measured with %i wavelength pairs!'], ...
                nStacksLoaded, nWavelengthPairs );
            
            msgbox(errStr,'Error','error');
            error(errStr);
        end
        
        tt = standards(ic).t;
        x = standards(ic).x;
        
        if nWavelengthPairs == 1
            % for single-wavelength standards, simply interpolate
            A(:,ic) = interp1(tt,x,t,'spline');
        else
            % for multi-wavelength standards do something else
            for ip = 1:nWavelengthPairs
                % index to time delays for this standard and wavelength
                % pair
                stdInd = find( standards(ic).tindex == ip );
                
                % index to time delays for delay stack acquired with this
                % wavlength pair
                stkInd = find( delayIndex == ip );
                
                A(stkInd, ic) = interp1( ...
                    tt(stdInd), x(stdInd), t(stkInd), ...
                    'spline' );
            end
        end
        
    end

    setappdata(fig,'Y',A); 

    % QR orthogonalize
    [Q,~] = qr(A,0);

    % normalize
    %Q = Q ./ repmat(sqrt(sum(Q.^2)) / nt, [nt,1]);

    setappdata(fig,'Q',Q);

    % work out projection back to the original components
    D = (Q.'*A);
    setappdata(fig,'D',D);
% end % orthonormalize()

% --- Executes on button press in pb_decompose.
function pb_decompose_Callback(hObject, eventdata, handles)
    
    % get figure from hObject's parent rather than gcbf,
    % allowing this function to operate even if called by some means other
    % than MATLAB event handling
    fig = get(hObject,'Parent');
%     --> This Makes Phasor segmentation crash. 
%         Changed some getappdata handles to gcbf to run callbacks - FER

    % get settings
    doReducedRankAppx = getappdata( gcbf, 'reduceRankEnable' );
    baselineMethod = getappdata( gcbf, 'baselineMethod' );
    

    % retrieve pump-probe standards for each decomposition channel
    standards = getappdata( gcbf, 'standards' );
    ncomp = length(standards);

    % generate an orthonormal basis from the standards
    [Q,D] = orthonormalize( gcbf, handles );

    % arrange imagestack into an npx by nt array
    X = getappdata(gcbf,'imageStack');
    
    % store an unlatered copy of the image stack
    X0 = X;
    [nr,nc,nt] = size(X);

    % if preprocessing requires reduced-rank approximation
    if doReducedRankAppx ...
            || strcmp( baselineMethod, 'BASELINE_REDUCEDRANK' )
        
        r = getappdata( fig, 'reduceRankK');
        X_delays = reshape(X,[nr*nc, nt]);

        [nr,nc,nt] = size(X);
        [U,S,V] = svd(X_delays,0);
        s = diag(S);
        s((r+1):end) = 0;
        S = diag(s);
        X_delays = U*S*V';
        Xrr = reshape(X_delays,[nr,nc,nt]);
    end
    
    if doReducedRankAppx
        % replace image stack with reduced-rank approximation
        X = Xrr;
    end

    % find indices at beginnings of delay stacks
    nStacksLoaded = getappdata( fig, 'nStacksLoaded' );
    tind = getappdata( fig, 'delayIndex' );
    
    % Aubtract baseline/offset according to user's preference of method
    % (After this block of code, both X and X0 will have the baseline
    % subtracted. We do this for both stacks, so that at the end, the
    % RMS error calculations take into account this baseline subtraction.)
    switch baselineMethod
        case {'BASELINE_INDEPENDENTCHANNEL', 'BASELINE_REDUCEDRANK'}
            % estimate offset for each appended delay stack independently

            % if multiple stacks have been appended, we estimate the
            % baseline of each one independently.
            for k = 1:nStacksLoaded
                % find index to first slice in this stack
                iFirst = find(tind == k,1,'first');
                
                % get number of slices in this stack
                nt_this = sum( tind == k );
                
                if strcmp( baselineMethod, 'BASELINE_REDUCEDRANK' )
                    % estimate offset from the reduced-rank form
                    off = Xrr(:,:,iFirst);
                else
                    % otherwise, estimate offset from original stack
                    off = X(:,:,iFirst);
                end

                % subtract the offset
                X(:, :, tind == k)= X(:, :, tind == k) ...
                                        - repmat(off,[1,1,nt_this]);

                % also subtract offset from original data, so we can make a
                % fair comparison later for measuring RMS fit error
                X0(:, :, tind == k)= X0(:, :, tind == k) ...
                                        - repmat(off,[1,1,nt_this]);
            end
            
        case 'BASELINE_FRAMEAVERAGE'
            % if multiple stacks have been appended, we estimate the
            % baseline of each one independently.
            for k = 1:nStacksLoaded
                % find index to first slice in this stack
                iFirst = find(tind == k,1,'first');

                % use the average of the first frame as the offset
                off = mean(mean(X(:,:,iFirst)));
                
                % subtract the offset
                X(:, :, tind == k)= X(:, :, tind == k) - off;
                
                % also subtract offset from original data, so we can make a
                % fair comparison later for measuring RMS fit error
                X0(:, :, tind == k)= X0(:, :, tind == k) - off;

            end
    end % switch baselineMethod

    % convert image stack into list of delay scans,
    X_delays = reshape(X,[nr*nc, nt]);
    
    % initialize xDotz; each row of this will contain each delay scan's
    % projection onto each of the orthonormal basis vectors
    xDotz = zeros((nr*nc), ncomp);

    % project delay scans onto each component
    % (This can probably be replaced by a single matrix multiplication
    % instead of a for loop)
    for icomp = 1:ncomp
        % get this basis vector
        z = Q(:,icomp);

        % calculate dot product with this basis vector
        xDotz(:,icomp) = X_delays*z;

        % show the image representing the dot product
        axes(handles.ax_decomposed);
        imagesc(reshape(xDotz(:,icomp),[nr,nc]));
        colormap(gray);
        drawnow;
    end

    % now unmix each pixel. C is an npx by ncomp matrix
    % where each column is the pixels' projection onto the component
   
    if get( handles.Force_Pos_Proj, 'Value' )
        % After projections, force unmixing to be positive.
%        C = zeros(size(xDotz,1),size(xDotz,2));
%         h = waitbar(0,'Please wait...');
%         steps = size(xDotz,1);
%         for nonNegLoop = 1: steps;
%             C(nonNegLoop,:) = lsqnonneg(D,xDotz(nonNegLoop,:).').';
%             waitbar(nonNegLoop / steps)
%         end
%         close(h)
% THIS TAKES TOO LONG
        C = (D\xDotz.').';
%         
        
    else
        C = (D\xDotz.').';
    end
    
    if getappdata( fig, 'projThreshEnable' )
        ThreshProjInt = getappdata( fig, 'projThreshVal' );
        
        % the following code rejects pixels lacking a specified amount of
        % any component
        ThresholdMask = zeros(nr,nc);
        for ii=1:ncomp;
            ThresholdMask = ThresholdMask | (C2D(:,:,ii) > ThreshProjInt);
        end
        C2D = C2D.* repmat(ThresholdMask,[1,1,ncomp]); 

        % the following code is used for phasor masking
        if get( handles.ApplyProjThresh, 'Value' )
            X = X.* repmat(ThresholdMask,[1,1,size(X,3)]); 
        end

    end

    % post-processing:
    baselineChannels = getappdata( gcbf, 'baselineChannels' );
    activeChannels = setdiff( 1:ncomp, baselineChannels );
    
    % reject pixels with negative values
    negValRejectMethod = getappdata( gcbf, 'negValRejectMethod' );
    
    switch negValRejectMethod
        
        case 'NEGREJECT_INDEPCHANNEL'
            % set all negative values in non-baseline channels to zero
            
            C_ac = C( :, activeChannels);
            C_ac( C_ac < 0 ) = 0;
            C(:,activeChannels) = C_ac;

        case 'NEGREJECT_WHOLEPX'
            % go through each pixel, if any channels for that pixel are
            % negative, set all channels for that pixel to zero
            
            posvalMask = ones(nr*nc,1);

            for ii = activeChannels;
                posvalMask = posvalMask.*(C(:,ii) > 0.);
            end
            C = C.*repmat(posvalMask,[1,ncomp]);
            
            %C2D( C2D < 0 ) = 0; % reject negative mixing coefficients
            
            % THE FOLLOWING IS BROKEN; REQUIRES C2D
            %if get( handles.ApplyProjThresh, 'Value' )
            %    X = X.* repmat(posvalMask,[1,1,size(X,3)]);
            %end
            
        case 'NEGREJECT_NONE'
            % don't do any negative value rejection
            
            % Note that negative values are ignores in RGB image, thus
            % absolute value is taken.  No way to disntiguish between
            % a positive and negative projection with same magnitude--FER
            
        otherwise
            error( ['Invalid setting for negValRejectMethod: ', ...
                negValRejectMethod] );
    end
    
    

    % generate a new image stack, one slice per component
    C2D = reshape(C,[nr,nc,ncomp]);
    setappdata( gcbf, 'unmixedCoeffs', C2D );

    C2D_old = C2D;

     
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % calculate model fit RMS error (residuals)
    
    % figure out residuals
    % re-mix from composition
    Xremix = zeros((nr*nc),nt);
    Y = getappdata(gcbf,'Y');
    for icomp = 1:ncomp
        y = Y(:,icomp);
        Xremix = Xremix + ...
            repmat(C(:,icomp),[1,nt]).* repmat(y.',[(nr*nc),1]);
    end

    % remix from orthonormal basis
    Xremix_PC = zeros((nr*nc),nt);
    for icomp = 1:ncomp
        z = Q(:,icomp);
        Xremix_PC = Xremix_PC + ...
            repmat(xDotz(:,icomp),[1,nt]).*repmat(z.',[(nr*nc),1]);
    end

    %Xremix = Xremix_PC;
    X_delays = reshape(X0,[nr*nc, nt]);
    
    % square-error
    se = (X_delays - Xremix).^2;
    
    % image of root mean square error for each pixel
    rmse = sqrt(sum(se,2))/nt;
    
    % average RMS error of entire image
    rmse_total = sqrt( sum (se(:)) ) / (nt*nr*nc);
    
    fprintf('RMS Error: %f\n', rmse_total);
    
    % ??
    ser = (sqrt(sum(X_delays.^2,2))./rmse)/nt;
    axes(handles.ax_residuals);
    imagesc(reshape(rmse,[nr,nc]));
    set(gca,'ydir','normal');
    colorbar;

    IRGB = zeros(nr,nc,3);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
    % create a color composite image!
    % use only non-baseline channels
    axes(handles.ax_decomposed);

    for icomp = activeChannels
        % retrieve color from s
        color = standards(icomp).color;
        
        colorMatrix = repmat(shiftdim(color,-1),[nr,nc,1]);

        IRGB = IRGB + repmat(abs(C2D(:,:,icomp)),[1,1,3]).*colorMatrix;
        % Have to take absolute value. Other wise, if the 'do nothing'
        % option for negative values is selected, then it will be equal to
        % the 'whole pixel' option  
        % If values are negative, then options above will enable one to
        % ignore specific channel that is negative or whole pixel.  This
        % has to be done explicitly. 
        %  - FER
    end

    discardInk = 0;
    if discardInk
        mask = (IRGB(:,:,3) > IRGB(:,:,1)) | (IRGB(:,:,3) > IRGB(:,:,2));
        IRGB = IRGB.*repmat((1-mask),[1,1,3]);
        IRGB(:,:,3) = 0;
    end

    if 0
        IRGB = IRGB / max(IRGB(:));
    else
        m = mean(IRGB(:));
        s = std(IRGB(:));
        IRGB = IRGB / (m+10*s);
    end
    % Then display resulting orthonormal basis image
    imshow(IRGB);
    set(gca,'ydir','normal');

    % re-label the image axes after calls to imshow and imagesc
    set(handles.ax_decomposed,'tag','ax_decomposed');


    setappdata(fig,'Xoffsetcorrected',X);
    setappdata(fig,'Xremix',Xremix);
    setappdata(fig,'Xresiduals',X_delays - Xremix);
    setappdata(fig,'SER',ser);
    setappdata(fig,'RMSE',rmse);
    setappdata(fig,'rundecomp',0);
    setappdata(fig,'ProjectedData',C2D); 
    setappdata(fig,'ProjectedData_NoThresholding',C2D_old); 

    if getappdata( fig, 'enablePhasor' )
        ComputePhasors( fig )
    end

    % quantify composition
    % C is an npx by ncomp matrix
    F = quantifyComposition( C );
    setappdata(gcf,'fractionalComposition',F);
% end % pb_decompose_Callback()


% --- Executes on button press in pb_loadStandards.
function pb_loadStandards_Callback(hObject, eventdata, handles)
% prompt user for standards .xml file, and then load it
    if ispref('puprisa','standardsSetsDir')
        standardsSetsDir = getpref('puprisa','standardsSetsDir');
    else
        % locate puprisa install directory
        puprisaDir = which('puprisa');
        ifs = strfind(puprisaDir,filesep);
        puprisaDir = puprisaDir(1:ifs(end));

        standardsSetsDir = [puprisaDir,filesep,'..',filesep,'standardsSets'];
        addpref('puprisa','standardsSetsDir',standardsSetsDir);
    end

    % prompt user for standards set
    [fileName, pathName] ...
        = uigetfile('*.xml', ...
                    'Select pump-probe spectral standards set',...
                    [standardsSetsDir,filesep,'*.xml']);
    
    % abort if user clicked 'cancel'
    if isequal(fileName,0); return; end;

    % double check that a .xml file was specified
    if ~strcmp(fileName((end-3):end),'.xml')
        errordlg('Invalid file type. Please select a set of standards specified by a .xml file from the puprisa standardsSets directory.');
        return;
    end

    % save to preferences
    setpref('puprisa','standardsSetsDir',pathName);
    if ~ispref('puprisa','standardsFile')
        addpref('puprisa','standardsFile',[])
    end 
    setpref('puprisa','standardsFile',fileName);

    % load standards from the specified file
    loadStandards( gcbf, [pathName,filesep,fileName] );
% end % pb_loadStandards_Callback()

function loadStandards( fig, setFileName )
% read a set of standards specified in an xml file

    % clear currently-loaded standards
    setappdata( fig, 'standards_original', [] );
    setappdata( fig, 'standards', [] );
    setappdata(fig,'standardsLoaded',0);

    % open the XML document
    xDoc = xmlread(setFileName);
    
    standardsSetsDir = getpref('puprisa','standardsSetsDir');

    standardsDir = [standardsSetsDir,filesep,'..',filesep,'standards'];
    
    % get top node
    pumpProbeStandardsSetNode = xDoc.getDocumentElement;
    
    % get all entry nodes
    topNode = pumpProbeStandardsSetNode.getChildNodes;
    
    % get nodes for each 'standard' entry
    standardNodes = topNode.getElementsByTagName('pumpProbeStandard');
    
    
    for inode = 0:(standardNodes.getLength - 1)
        node = standardNodes.item(inode);
        
        % name for labeling this standard
        s{inode+1}.name = ...        
            char(node.getElementsByTagName('name').item(0).getTextContent);
        
        % color
        s{inode+1}.color = [ ...
            str2num(node.getElementsByTagName('color_red').item(0).getTextContent), ...
            str2num(node.getElementsByTagName('color_green').item(0).getTextContent), ...
            str2num(node.getElementsByTagName('color_blue').item(0).getTextContent), ...
            ];
        
        % scaling factor
        s{inode+1}.scale = ...        
            str2num(node.getElementsByTagName('scale').item(0).getTextContent);
        
        % name of the .mat file
        fileName = ...        
            char(node.getElementsByTagName('fileName').item(0).getTextContent);
        
        % load the standard from the .mat file
        dat = load([standardsDir,filesep,fileName]);
        s{inode+1}.t = dat.t;
        s{inode+1}.x = dat.x;
        
        % is this a multi-wavelength standard?
        if isfield(dat,'wlpairs') && isfield(dat,'tindex')
            s{inode+1}.wlpairs = dat.wlpairs;
            s{inode+1}.tindex = dat.tindex;
        end
    end
    
    handles = guihandles(fig);
    axes( handles.ax_components );
    cla;

    for is = length(s):-1:1
        s{is}.lineobj = ...
            line(s{is}.t, s{is}.x*s{is}.scale, 'color', s{is}.color, ...
                'DisplayName', s{is}.name );
    end
    
    % save these standards in appdata
    setappdata(fig,'standardsLoaded',1);
    setappdata(fig,'standards',cell2mat(s));
    
    % add extra standards if necessary, for baseline channels
    if strcmp( getappdata( fig, 'baselineMethod' ), ... 
                'BASELINE_INDEPENDENTCHANNEL' )
        
        addBaselineChannels( fig );
    end
    
    % update legend in the channels/standards plot
    updateChanLegend( fig );
    
    % now enable decompose button if an image has been loaded
    X = getappdata(fig,'imageStack');
    if ~isempty(X)
        set(handles.pb_decompose,'Enable','on');
    end

    % update the list of channel names
    S = cell2mat(s);
    channelNames = {S.name};
    updateChannelNames( fig, channelNames  );
    
    % normalize standards, if requested
    if getappdata( fig, 'standardsNormalized' )
        normalizeStandards( fig );
    end

% --- Called by any functions that make changes to the pump/probe standards 
%     or baseline channels,
function updateChanLegend( fig )
    handles = guihandles(fig);
    axes( handles.ax_components );
    legend('off');
    l=legend('show');
    tx=findobj(l,'type','text');
    set(tx,'color','w');
    
% end % updateChanLegend()

function updateChannelNames( fig, channelNames )
    % update the listbox of channel names
    % and enable appropriate user interface elements
    handles = guihandles( fig );
    set( handles.lb_channelNames, 'String', channelNames );
    set( handles.lb_channelNames, 'Enable', 'on' );
    
    
    
    
    
% --- Executes during object creation, after setting all properties.
function figure1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on button press in pbSaveImage.
function pbSaveImage_Callback(hObject, eventdata, handles)

    % generate a default filename, using the first loaded image stack as a
    % starting point
    fileNames = getappdata(gcbf,'fileNames');
    fileName = fileNames{1};
    
    fileName = [fileName(1:(end-4)),'.png'];

    % recall last accessed directory
    if ispref('puprisa','imageExportPath')
        pathName = getpref('puprisa','imageExportPath');
    else
        pathName = pwd();
        addpref('puprisa','imageExportPath',pathName);
    end

    % prompt user for imageStack file to be loaded
    [fileName, pathName] = uiputfile( '*.png','Export PNG as:',...
                                      [pathName,filesep,fileName]);
    if isequal(fileName, 0); return; end

    setpref('puprisa','imageExportPath',pathName);

    axes(handles.ax_decomposed);

    I = getimagergb();

    imwrite(I,[pathName,filesep,fileName],'PNG');



% --- Executes on selection change in lb_channelNames.
function lb_channelNames_Callback(hObject, eventdata, handles)
% hObject    handle to lb_channelNames (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns lb_channelNames contents as cell array
%        contents{get(hObject,'Value')} returns selected item from lb_channelNames


% --- Executes during object creation, after setting all properties.
function lb_channelNames_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lb_channelNames (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox1.
function checkbox1_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox1


% --- Executes on button press in checkbox2.
function checkbox2_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox2


% --- Executes on button press in pushbutton6.
function pushbutton6_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in checkbox3.
function checkbox3_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox3



function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double


% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pbViewResiduals.
function pbViewResiduals_Callback(hObject, eventdata, handles)
% hObject    handle to pbViewResiduals (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
X = getappdata(gcf,'imageStack');
Xresiduals = getappdata(gcf,'Xresiduals');
[nr,nc,nt] = size(X);
puprisa_viewChannel(reshape(Xresiduals,[nr,nc,nt]),getappdata(gcbf,'delays'),'residuals',[]);


% --- Executes on button press in pbViewFit.
function pbViewFit_Callback(hObject, eventdata, handles)
% hObject    handle to pbViewFit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
X = getappdata(gcf,'imageStack');
Xremix = getappdata(gcf,'Xremix');
[nr,nc,nt] = size(X);
puprisa_viewChannel(reshape(Xremix,[nr,nc,nt]),getappdata(gcbf,'delays'),'fit',[]);


% --- Executes on button press in cb_reduceRank.
function cb_reduceRank_Callback(hObject, eventdata, handles)
% hObject    handle to cb_reduceRank (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    reduceRankEnable = get(hObject,'Value');
    setWithPref( gcbf, 'lindecomp_reduceRankEnable', ...
        'reduceRankEnable', reduceRankEnable );
% end % cb_reduceRank_Callback()


function ed_rank_Callback(hObject, eventdata, handles)
% hObject    handle to ed_rank (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    % recall previous setting
    r0 = getappdata( gcbf, 'reduceRankK' );

    % convert to string
    r = str2double( get( hObject, 'String' ) );
    
    % valid entries are positive integers
    if isnan( r ) || (round(r) ~= r) || (r < 1)
        % if user input is invalid, just restore the textbox with the
        % previous setting
        set( hObject, 'String', r0 );
    else
        % if entry was valid, change it
        setWithPref( gcbf, 'lindecomp_reduceRankK', ...
            'reduceRankK', r );
    end
% end % ed_rank_Callback() 

% --- Executes during object creation, after setting all properties.
function ed_rank_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ed_rank (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox5.
function checkbox5_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox5


% --- Executes on button press in pb_projectStandards.
function pb_projectStandards_Callback(hObject, eventdata, handles)
% hObject    handle to pb_projectStandards (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    Y = getappdata(gcbf,'Y');
    ll = findobj(handles.ax_components,'type','line');

    X = getappdata(gcf,'imageStack');
    [nr,nc,nt] = size(X);
    X_delays = reshape(X,[nr*nc,nt]);

    [U,S,V] = svd(X_delays,0);

    % full rank
    M = V'*Y;
    ncomp = length(ll);

    rank = ncomp+1;

    Y = V(:,1:rank)*M(1:rank,:);
    setappdata(gcbf,'Y',Y);
    t = getappdata(gcbf,'delays');
    for icomp = 1:ncomp
        set(ll(icomp),'xdata',t,'ydata',Y(:,icomp));
    end
    drawnow;
% end % pb_projectStandards_Callback()

% --- Executes during object creation, after setting all properties.
function ed_threshold_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ed_threshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in NegValIgnore.
function NegValIgnore_Callback(hObject, eventdata, handles)
% hObject    handle to NegValIgnore (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


    setWithPrefEnum( gcbf, 'lindecomp_negValRejectMethod', ...
            'negValRejectMethod', hObject );

    setappdata(gcbf,'rundecomp',1);
% end % NegValIgnore_Callback()



% --- Executes during object creation, after setting all properties.
function NegValIgnore_CreateFcn(hObject, eventdata, handles)
% hObject    handle to NegValIgnore (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in ShowHist.
function ShowHist_Callback(hObject, eventdata, handles)
% hObject    handle to ShowHist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of ShowHist


    if getappdata(gcbf,'rundecomp')
        pb_decompose_Callback(hObject, eventdata, handles);
    end

    C2D = getappdata(gcbf,'ProjectedData');

    DiscardAllNegValPix = getappdata(gcbf,'IgnoreNegativeValue');
    if isempty(DiscardAllNegValPix)
        DiscardAllNegValPix = 'Single Neg element';
    end
    
    
    if size(C2D,3)==2;% Assume its eu- and pheo-melanin, respectively.
    
      totalImage = (sum(C2D,3));
      fracImageMasked = C2D(:,:,1)./totalImage;


      figure('Name','Quantitative Eu/Pheo Image and Histogram');

        subplot(2,3,[1,2,4,5]);

          title('Fractional Eumelanin')
            im = imagesc(fracImageMasked*100);
            set(gca,'ydir','normal','xtick',[],'ytick',[]);

            set(gca,'color','k');
            colormap(cmap(256,[0 1 0], [1 1 0], [1 0 0]));
            caxis([10,100]);
            axis image;
            h=colorbar;
            ylabel(h, 'Eumelanin percentage');
         % set up transparency to indicate melanin concentration
            set(im,'alphadata',totalImage);
            set(im,'alphadatamapping','scaled');
            if min(totalImage(:)) < 0
                alim([0,max(totalImage(:))/5]);
            else
                alim([min(totalImage(:)),max(totalImage(:))/5]);
            end


        subplot(2,3,3);
           % figure;

            [n,bin] = histc(fracImageMasked(:), linspace(-.1,1.1));

            % now accumulate weights
            w = 0*n;
            for ii = 1:length(n);%2:length(n)-1;%
                w(ii) = sum(totalImage( bin == ii) );
            end

            h = bar(linspace(-.1,1.1),w);
            set(h,'tag','histogramBarGraph');
            xlabel('Eumelanin fraction');
            ylabel('Sum pixel intensity');
            title(['Intensity-weighted Histogram, ignoring Neg val: ', DiscardAllNegValPix])
            xlim([-.1,1.1])
            axis square
            
            CoM = sum(w'.*linspace(-.1,1.1))./sum(w);
           subplot(2,3,[1,2,4,5]);
            title(['Mean Eu frac: ', num2str(CoM*100),'%']);
            
         subplot(2,3,6);   
            imagesc(1-(totalImage==0))
            axis xy;axis square
    else
        beep
        display('Only for Eu and Pheo Mixtures')
        display('use euPheo_cuvette-pu720pr810 Standard Set')

    end
% end % ShowHist_Callback()

% --- Executes on button press in ProjThresh.
function ProjThresh_Callback(hObject, eventdata, handles)
% hObject    handle to ProjThresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    x = get(hObject,'Value');
    
    setWithPref( gcbf, 'lindecomp_projThreshEnable', ...
        'projThreshEnable', x );
    
% end % ProjThresh_Callback()

function ThresholdVal_Callback(hObject, eventdata, handles)
% hObject    handle to ThresholdVal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    val0 = getappdata( gcbf, 'projThreshVal' );
    
    val = str2double( get( hObject,'String' ) ); 
    
    if ~isnan( val )
        % if a valid number has been entered, remember the setting
        setWithPref( gcbf, 'lindecomp_projThreshVal', 'projThreshVal', ...
            val );
    else
        % otherwise, overwrite the textbox with the last good setting
        val = val0;
    end
    
    set( hObject, 'string', num2str( val ) );
% end % ThresholdVal_Callback()
    

% --- Executes during object creation, after setting all properties.
function ThresholdVal_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ThresholdVal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in Freq.
function Freq_Callback(hObject, eventdata, handles)
% hObject    handle to Freq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns Freq contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Freq

switch( get(hObject,'Value') )
    case 1
       cm = 0.25;
    case 2
       cm = 0;
end

setappdata(gcbf,'PhasorFrequency',cm);

% --- Executes during object creation, after setting all properties.
function Freq_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Freq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in DisplayMethod.
function DisplayMethod_Callback(hObject, eventdata, handles)
% hObject    handle to DisplayMethod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    switch( get(hObject,'Value') )
        case 1
            cm = 'Int weighted';
        case 2
            cm = 'Threshold';
    end

    setappdata(gcbf,'PhasorDisplayMethod',cm);
% end % DisplayMethod_Callback()

% --- Executes during object creation, after setting all properties.
function DisplayMethod_CreateFcn(hObject, eventdata, handles)
% hObject    handle to DisplayMethod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in MaskType.
function MaskType_Callback(hObject, eventdata, handles)
% hObject    handle to MaskType (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    switch( get(hObject,'Value') )
        case 1
            cm = 'Skin eu/pheo';
        case 2
           cm = 'Ocular eu/pheo';
        case 3
            cm = 'Ocular eu/pheo/Hb';
        case 4
            cm = 'Ocular Hb';
        case 5
            cm = 'Manual';

    end

    setappdata(gcbf,'PhasorMaskType',cm);
% end % MaskType_Callback()


% --- Executes during object creation, after setting all properties.
function MaskType_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MaskType (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in ViewPhasors.
function ViewPhasors_Callback(hObject, eventdata, handles)
% hObject    handle to ViewPhasors (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    N = getappdata(gcbf,'NPhasorCounts');
    DispMethod = getappdata(gcbf,'PhasorDisplayMethod');

    x = linspace(-1,1,256);

    figure('Name', 'Phasor Plot');
    title(['Histogram using method: ',DispMethod])
    imagesc('xdata',x,'ydata',x,'cdata',N.'); 
    load 'contourColorMap'
    colormap(contourColorMap);
    xlim([-1,1]);ylim([-1,1]);
    xlabel('g(\omega)'); ylabel('s(\omega)');
    xx = linspace(0,1,512);
    semicircle = sqrt(xx.*(1-xx));
    line(xx,semicircle,'LineStyle','--','Color','k');
    line(-xx,-semicircle,'LineStyle','--','Color','k');

    for drawlines = -.8:.2:.8;
        line(drawlines,-1:.05:1,'Color','k','Linewidth',1);
        line(-1:.05:1,drawlines,'Color','k','Linewidth',1);
    end
% end % ViewPhasors_Callback()


% --- Executes on button press in ViewPhasorMask.
function ViewPhasorMask_Callback(hObject, eventdata, handles)
% hObject    handle to ViewPhasorMask (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


    f = getappdata(gcbf,'PhasorFrequency');

    MaskType = getappdata(gcbf,'PhasorMaskType');
 
    if isempty(MaskType)
        MaskType = 'Skin eu/pheo';
    end
        
    if f ==0.25;
                
         switch MaskType

             case 'Skin eu/pheo';
                    load PhasorMasks\SkinEuPheoPhasorROI.mat
             case 'Ocular eu/pheo';
                    load PhasorMasks\OcularEuPheoPhasorROI.mat
             case 'Ocular eu/pheo/Hb';
                    load PhasorMasks\OcularHbEuPheoPhasorROI.mat
             case 'Ocular Hb';
                    load PhasorMasks\OcularHbPhasorROI.mat
             case 'Manual';
                 beep
                 display('No standard Phasor Mask available')
                 display('Run Segment to define phasor mask')
                 display('Change Mask type ')
                 ViewPhasors_Callback(hObject, eventdata, handles)
                 [PhasorMask,xi,yi] = roipoly;

         end
         
         if size(PhasorMask,1)>1;
           ViewPhasors_Callback(hObject, eventdata, handles)
           line(xi,yi,'LineWidth',3,'color','r');
           
         end 
         
    elseif f==0;
        beep
                 display('No standard Phasor Mask available')
                 display('Run Segment to define phasor mask')
                 display('or change frequency')
                 display('to define new f, Decompose or Segment')
                 PhasorMask = 0; xi=0; yi=0;
    else
        beep
        ViewPhasors_Callback(hObject, eventdata, handles)
        [PhasorMask,xi,yi] = roipoly;    
        line(xi,yi,'LineWidth',3,'color','r');
    end
 
    setappdata(gcbf,'PhasorMask',PhasorMask)
    setappdata(gcbf,'xi',xi)
    setappdata(gcbf,'yi',yi)
% end % ViewPhasorMask_Callback()

% --- Executes on button press in segment.
function segment_Callback(hObject, eventdata, handles,ApplyMask)
% hObject    handle to segment (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    % run decompose to ensure that all values are up to date
    % this also updates 'computephasors'
 
   if nargin<4; ApplyMask  = 0 ;end
 
    pb_decompose_Callback(hObject, eventdata, handles);
    figureID = gcf;

    g = getappdata(gcbf,'g');
    s = getappdata(gcbf,'s');

    ViewPhasorMask_Callback(hObject, eventdata, handles);

    PhasorMask = getappdata(gcbf,'PhasorMask');
%     xi = getappdata(gcbf,'xi');
%     yi = getappdata(gcbf,'yi');


    [py,px] = size(PhasorMask);
    npx = sqrt(size(g,1));% g is nr*nc,  and nr=nc
    G = reshape(g,[npx, npx]);
    S = reshape(s,[npx, npx]);

    ImMask = zeros(npx,npx);
    for pxx = 1:npx;
        for pyy = 1:npx;
            gp = round((G(pyy,pxx) * px/2) + px/2) ; %maps from phasor val to pixel val of phasor plot
            sp = round((S(pyy,pxx) * py/2) + py/2) ;
            if gp==0; gp=1;elseif sp==0; sp=1;end
            ImMask(pyy,pxx) = PhasorMask( sp, gp);

        end
    end

 if ApplyMask;

      X = getappdata(gcbf,'imageStack');
      X = X.*repmat(ImMask,[1,1,size(X,3)]);
      setappdata(gcbf,'imageStack',X);  
      
     figure(figureID);
     
     pb_decompose_Callback(hObject, eventdata, handles);
     

 else

        if get( handles.ApplyProjThresh, 'Value' )
            C2D = getappdata(gcbf,'ProjectedData'); 
        else
            C2D = getappdata(gcbf,'ProjectedData_NoThresholding'); 
        end


        if size(C2D,3)==2;% Assume its pheo- and eu-melanin, respectively.

            C2D(:,:,1)= C2D(:,:,1).*ImMask;
            C2D(:,:,2)= C2D(:,:,2).*ImMask;

            totalImage = (sum(C2D,3));
            fracImageMasked = C2D(:,:,1)./totalImage;% assume eumlanin is on first channel

          ShowHist_Callback(hObject, eventdata, handles)

            %%%%
            figure('Name','Phasor Segmented Quantitative Eu/Pheo Image and Histogram');

            subplot(2,3,[1,2,4,5]);

                im2 = imagesc(fracImageMasked*100);
                set(gca,'ydir','normal','xtick',[],'ytick',[]);

                set(gca,'color','k');
                colormap(cmap(256,[0 1 0], [1 1 0], [1 0 0]));
                caxis([10,100]);
                axis image;
                h=colorbar;
                ylabel(h, 'Eumelanin percentage');
             % set up transparency to indicate melanin concentration
                set(im2,'alphadata',totalImage);
                set(im2,'alphadatamapping','scaled');
                if min(totalImage(:)) < 0
                    alim([0,max(totalImage(:))/5]);
                else
                    alim([min(totalImage(:)),max(totalImage(:))/5]);
                end


            subplot(2,3,3);
            % figure;

                [n,bin] = histc(fracImageMasked(:), linspace(-.1,1.1));

                % now accumulate weights
                w = 0*n;
                for ii = 1:length(n);%2:length(n)-1;%
                    w(ii) = sum(totalImage( bin == ii) );
                end

                h = bar(linspace(-.1,1.1),w);
                xlabel('Eumelanin fraction');
                ylabel('Sum pixel intensity');
                title('Intensity-weighted Histogram')
                xlim([-.1,1.1])
                axis square

                CoM = sum(w'.*linspace(-.1,1.1))./sum(w);
              subplot(2,3,[1,2,4,5]);
                title(['Mean Eu frac: ', num2str(CoM*100),'%']);

            subplot(2,3,6);

                if get( handles.ApplyProjThresh, 'Value' )
                     title('Phasor and Proj Thersh Mask')
                else
                    title('Phasor Mask')
                end
                imagesc(ImMask.*(1-(totalImage==0)));
                axis xy;axis square

        %         figure;image(totalImage)
        end
 end
% end % segment_Callback()

function ComputePhasors(fig)

    f = getappdata(fig,'PhasorFrequency');

    if f == 0;
        beep
        f = input('input frequency (in THz):  ');
        setappdata(fig,'PhasorFrequency',f);
    elseif isempty(f);
        f = 0.25;
        setappdata(fig,'PhasorFrequency',f);
    end

    display(['freq: ', num2str(f)]);


    X = getappdata(fig,'Xoffsetcorrected');

    t = getappdata(fig,'delays');

    [nr,nc,tcomp] = size(X);

    omega = 2*pi*f;

    X = reshape(X,[nr*nc,tcomp]);

    Xint = (sum(abs(X),2));

    g = sum(X.*repmat(cos(omega*t),[nr*nc,1]),2)./Xint;
    s = sum(X.*repmat(sin(omega*t),[nr*nc,1]),2)./Xint;

    g(Xint == 0) = 0;
    s(Xint == 0) = 0;

       DispMethod = getappdata(fig,'PhasorDisplayMethod');
        if isempty(DispMethod)
            DispMethod = 'Int weighted';
        end

     switch DispMethod
        case 'Int weighted'
            I = sum(sqrt(X.^2),2);
            x = linspace(-1,1,256);

            gr = interp1(x,1:numel(x),g,'nearest')';
            sr = interp1(x,1:numel(x),s,'nearest')';

            N = accumarray([gr' sr'], I.^2, [numel(x) numel(x)]);

        case 'Threshold'
            NN = 1;
            X2sum = sum(abs(X),2);
            X2mean = mean(X2sum(:));
            X2stdev = std(X2sum(:));
            IntImMask = (X2sum > (X2mean + NN*X2stdev));

            x=linspace(-1,1,256);
            ctrs{1} = x;
            ctrs{2} = x;
            [N,C] = hist3([g(IntImMask>0),s(IntImMask>0)],ctrs);

     end

    setappdata(fig,'g',g);
    setappdata(fig,'s',s);
    setappdata(fig,'NPhasorCounts',N);
% end % ComputePhasors()
 
 

% --- Executes on button press in ApplyProjThresh.
function ApplyProjThresh_Callback(hObject, eventdata, handles)
% hObject    handle to ApplyProjThresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of ApplyProjThresh


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over ApplyProjThresh.
function ApplyProjThresh_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to ApplyProjThresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in ProjHist.
function ProjHist_Callback(hObject, eventdata, handles)
% hObject    handle to ProjHist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    C2D = getappdata(gcbf,'ProjectedData'); 
    C2Deu = C2D(:,:,2);
    C2Dpheo = C2D(:,:,1);
    C2D_old = getappdata(gcbf,'ProjectedData_NoThresholding'); 

    C2D_oldeu = C2D_old(:,:,2);
    C2D_oldpheo = C2D_old(:,:,1);

    ThreshProjInt = str2num( get( handles.ThresholdVal, 'String' ) );

    DiscardAllNegValPix = getappdata(gcbf,'IgnoreNegativeValue');

    if isempty(DiscardAllNegValPix)
        DiscardAllNegValPix = 'Single Neg element';
    end
        

    figure('Name','Hist of Projection to components');
        subplot(2,2,1);
            [n , xout] = hist(C2D_oldeu(:),2*1024);
            bar(xout,n,'hist');
            title('original Eu hist projection ')
            xlim([mean(C2D_oldeu(:))-5*std(C2D_oldeu(:)),mean(C2D_oldeu(:))+5*std(C2D_oldeu(:))])
            ylim([0,max(n)+2*std(n)])
           text(mean(C2D_oldeu(:))-4.9*std(C2D_oldeu(:)), max(n)+1.25*std(n),...
               ['StDev of all Proj Int data: ',num2str(std(C2D_old(:)))],...
                'BackgroundColor',[.9 .9 .9]);     
            
        subplot(2,2,2);
            [n , xout] = hist(C2Deu(:),2*1024);
            bar(xout,n,'hist');

        if get( handles.ProjThresh, 'Value' )
            title('Eu hist with threshold projection')
            xlim([ThreshProjInt,mean(C2Deu(:))+5*std(C2Deu(:))])
        else
            title('Eu hist with Neg Val ommitted ')
            
            switch DiscardAllNegValPix
                case 'Whole Pixel'
                    xlim([xout(2),mean(C2Deu(:))+5*std(C2Deu(:))])
                case 'Single Neg element'
                    xlim([xout(2),mean(C2Deu(:))+5*std(C2Deu(:))])
                case 'None'
                    xlim([mean(C2Deu(:))-5*std(C2Deu(:)),mean(C2Deu(:))+5*std(C2Deu(:))])
             end

        end

        subplot(2,2,3);
            hist(C2D_oldpheo(:),2*1024);title('original Pheo projection hist')
            xlim([mean(C2D_oldpheo(:))-5*std(C2D_oldpheo(:)),mean(C2D_oldpheo(:))+5*std(C2D_oldpheo(:))])
        subplot(2,2,4);
            [n , xout] = hist(C2Dpheo(:),2*1024);
            bar(xout,n,'hist');

        if get( handles.ProjThresh, 'Value' )
            title('Pheo hist with threshold projection')
            xlim([ThreshProjInt,mean(C2Dpheo(:))+5*std(C2Dpheo(:))])
        else
            title('Pheo hist with Neg Val ommitted ')
            switch DiscardAllNegValPix
                case 'Whole Pixel'
                    xlim([xout(2),mean(C2Dpheo(:))+5*std(C2Dpheo(:))])
                case 'Single Neg element'
                    xlim([xout(2),mean(C2Dpheo(:))+5*std(C2Dpheo(:))])
                case 'None'
                    xlim([mean(C2Dpheo(:))-5*std(C2Dpheo(:)),mean(C2Dpheo(:))+5*std(C2Dpheo(:))])
             end

        end
% end % ProjHist_Callback()

function F = quantifyComposition( C )
    % calculates overall fraction of each component in the image stack
    % C is an npx by ncomp matrix calculated by linear unmixing
    
    [nPx, nComp] = size( C );
    
    % discard any remaining negative pixels
    C( C < 0 ) = 0;
    
    F = sum( C );
    F = F / sum( F );
    
% end % quantifyComposition()

% --- Executes on button press in pb_testImage.
function pb_testImage_Callback(hObject, eventdata, handles)
    % make a test image from the standards we have loaded
    
    ll = findobj(handles.ax_components,'type','line');
    fig = get(hObject,'Parent');

    ncomp = length(ll);
    
    t = get(ll(1),'XData');
    nt = length(t);
    
    imageStack = zeros(ncomp, ncomp, nt);
    
    for ir = 1:ncomp
        rowComponent = get(ll(ir),'YData');
        for ic = 1:ncomp
            colComponent = get(ll(ic), 'YData');
            
            imageStack(ir,ic,:) = 0.5*rowComponent + 0.5*colComponent;
        end
    end
        
    setappdata(gcbf,'imageStack',imageStack);
    setappdata(gcbf,'imageStack_original',imageStack);
    setappdata(gcbf,'delays',t);

    % go ahead and run decomposition, if we're set up to do so
    if getappdata(gcbf,'standardsLoaded') == 1
        set(handles.pb_decompose,'Enable','on');
        pb_decompose_Callback(hObject, eventdata, handles);
    end
% end % pb_testImage_Callback ()


% --- Executes on button press in pb_exportChannels.
function pb_exportChannels_Callback(hObject, eventdata, handles)
% hObject    handle to pb_exportChannels (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    C = getappdata(gcbf,'unmixedCoeffs');
    if ~isempty( C )
        chans = getappdata(gcbf,'standards');

        for iChan = 1:numel(chans)
            n = chans(iChan).name;
            assignin('base',n,C(:,:,iChan));
        end
    end
% end % pb_exportChannels_Callback()



function ed_channelGain_Callback(hObject, eventdata, handles)
% hObject    handle to ed_channelGain (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ed_channelGain as text
%        str2double(get(hObject,'String')) returns contents of ed_channelGain as a double


% --- Executes during object creation, after setting all properties.
function ed_channelGain_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ed_channelGain (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pb_viewChannels.
%     Opens a new figure window, and plots the decomposition channels in
%     separate subplots. Sometimes it's easier to diagnose unmixing
%     problems with this view, rather than the composite RGB image
function pb_viewChannels_Callback(hObject, eventdata, handles)
% hObject    handle to pb_viewChannels (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    
    standards = getappdata( gcbf, 'standards' );
    nStandards = length( standards );
    baselineChans = getappdata( gcbf, 'baselineChannels' );
    
    % retreive mixing coefficients found by linear decomposition
    C = getappdata( gcbf, 'unmixedCoeffs' );
    
    nChnPerFig = 4;    % how many channels to display in one fig
    nCol = 2;
    nRow = 2;
    
    % number of figure windows needed to display all channels
    nFig = ceil( nStandards / nChnPerFig );
    
    % generate each new figure, and display the channels
    startChnThisFig = 1; % first channel to display in the current figure
    for iFig = 1:nFig
        
        % the last channel to be displayed in this figure:
        endChnThisFig = min( startChnThisFig + nChnPerFig - 1, nStandards );
        
        % make a new figure, with a title that tells which channels
        %  are displayed
        if( (endChnThisFig - startChnThisFig) > 1 )
            figName = sprintf('channels %i -- %i', startChnThisFig, ...
                                endChnThisFig );
        else
            figName = sprintf('channels %i', startChnThisFig );
        end
            
            
        figure( 'name', figName );
        colormap(gray);
        
        % iterate through all channels in this figure
        for iChan = 1:nChnPerFig
            
            thisChan = startChnThisFig + iChan - 1;
            isBaselineChan = ismember( thisChan, baselineChans );
            
            if thisChan > nStandards
                % stop if there are no more channels left to display
                break;
            end
            
            subplot(nRow, nCol, iChan);
            
            % retrieve the image for this channel
            C_this = C(:,:,thisChan);
            
            % plot it
            imagesc(C_this);
            
            % leave baseline channels auto-scaled
            if ~isBaselineChan
                % all other channels should be scaled from zero to
                % 4 standard devs above mean
                m = mean(C_this(:));
                if m > 0
                    % only attempt this auto-scaling if mean is > 0
                    % (which might not be the case if negative-value
                    % rejection was disabled)
                    caxis([0, m + 4*std(C_this(:))]);
                end
            end
            
            % label things
            axis image
            title( standards(thisChan).name );
            
        end
        
        startChnThisFig = startChnThisFig + nChnPerFig; 
    end
% end % pb_viewChannels_Callback()
    


% --- Executes on selection change in pmnu_baselineCorrection.
function pmnu_baselineCorrection_Callback(hObject, eventdata, handles)
% hObject    handle to pmnu_baselineCorrection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%     contents = cellstr(get(hObject,'String'));
%     switch( contents{get(hObject,'Value')} )
%         case 'Reduced-rank estimate'
%             baselineMethod = 'BASELINE_REDUCEDRANK';
%         case 'Independent baseline channels'
%             baselineMethod = 'BASELINE_INDEPENDENTCHANNEL';
%         case 'Frame average'
%             baselineMethod = 'BASELINE_FRAMEAVERAGE';
%         case 'None'
%             baselineMethod = 'BASELINE_NONE';
%         otherwise
%             error(['Invalid baseline correction method; someone may ',...
%                 'have edited the string of pmnu_baselineCorrection ',...
%                 'without updating the source code.']);
%     end
%   

    baselineMethod_old = getappdata( gcbf, 'baselineMethod' );

    setWithPrefEnum( gcbf, 'lindecomp_baselineMethod', ...
        'baselineMethod', hObject );
    
    baselineMethod_new = getappdata( gcbf, 'baselineMethod' );

    changeBaselineMethod( baselineMethod_old, baselineMethod_new );
% end % pmnu_baselineCorrection_Callback()


% --- Executes during object creation, after setting all properties.
function pmnu_baselineCorrection_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pmnu_baselineCorrection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

    % Hint: popupmenu controls usually have a white background on Windows.
    %       See ISPC and COMPUTER.
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end
    
% end % pmnu_baselineCorrection_CreateFcn()

% --- Executes after user changes baseline method. Stores the new
% preference, and ensures data structures are in place for the selected
% method.
function changeBaselineMethod( baselineMethod_old, baselineMethod_new )    
    
    if strcmp( baselineMethod_new, baselineMethod_old ) == 0
        % create/destroy baseline channels as needed
        if strcmp( baselineMethod_new, 'BASELINE_INDEPENDENTCHANNEL' )
            addBaselineChannels( gcbf );
        end

        if strcmp( baselineMethod_old, 'BASELINE_INDEPENDENTCHANNEL' )
            removeBaselineChannels();
        end
    end
% end % updateBaselineMethod()

% --- Executes when user switches baseline method, or when new standards
% are loaded. Adds channels to account for baseline offsets, as needed.
function addBaselineChannels( fig )
    % add baseline channels
    
    
    % if standards are normalized, undo it first, then restore when we're
    % done, so that the baseline channels will be normalized too
    standardsNormalized = getappdata( fig, 'standardsNormalized' );
    if standardsNormalized
        unnormalizeStandards( fig );
    end
    
    % get info
    % (this must be done /after/ unnnormalizing)
    s = getappdata( fig, 'standards' );
    nStandards = length(s);
    
    handles = guihandles( fig );
    axes( handles.ax_components );

    baselineChannels = [];
    if isfield( s(1), 'wlpairs' )
        % if multiple wavelength pairs, add one channel per wavelength pair
        wlpairs = s(1).wlpairs;
        [nWlPairs, ~] = size( wlpairs );
        
        for k = 1:nWlPairs
            % generate info required for an offset channel
            ind = nStandards+k;
        
            s(ind).name = sprintf( 'baseline%i', k );
            s(ind).color = [1 1 1];
            s(ind).scale = 1;
            s(ind).t = s(1).t;
            s(ind).x = double(s(1).tindex == k);
            s(ind).wlpairs = s(1).wlpairs;
            s(ind).tindex = s(1).tindex;
            
            % plot it along with the other standards
            s(ind).lineobj = ...
                line(s(ind).t, s(ind).x*s(ind).scale, ...
                    'color', s(ind).color, 'DisplayName', s(ind).name );

            % append index of this baseline channel to the list
            baselineChannels = [baselineChannels, ind];
        end
    else
        % otherwise just one master offset
        ind = nStandards+1;
        
        s(ind).name = 'baseline';
        s(ind).color = [1 1 1];
        s(ind).scale = 1;
        s(ind).t = s(1).t;
        s(ind).x = s(1).t * 0 + 1;

        % plot it along with the other standards
        s(ind).lineobj = ...
            line(s(ind).t, s(ind).x*s(ind).scale, ...
                'color', s(ind).color, 'DisplayName', s(ind).name );

        % remember index of this baseline channel
        baselineChannels = ind;
    end
    
    % update entries in appdata
    setappdata( fig, 'baselineChannels', baselineChannels );
    setappdata( fig, 'standards', s );
    
    if standardsNormalized
        normalizeStandards( fig );
    end
        
    
    % update legend in the channels/standards plot
    updateChanLegend( fig );
    
% end % addBaselineChannels()

% --- Executes when user switches baseline method, or when new standards
% are loaded. Removes channels to account for baseline offsets when these 
% are no longer needed.
function removeBaselineChannels()
    % baseline channels are marked
    baselineChannels = getappdata( gcbf, 'baselineChannels' );
    
    if ~isempty( baselineChannels )
        % retrieve standards in use 
        s = getappdata( gcbf, 'standards' ); 
        
        % retrieve backup copy of standards (in case the originals were
        % normalized, we'll need to remove the baseline channels from these 
        % too)
        s0 = getappdata( gcbf,'standards_original' );
        
        for k = 1:length( baselineChannels )
            % remove the line from the standards plot
            delete( s(baselineChannels(k)).lineobj );
        end
        
        % remove all the baseline standards entries
        s = s( setdiff( 1:length(s), baselineChannels ) );
        
        if ~isempty( s0 )
            s0 = s0( setdiff( 1:length(s0), baselineChannels ) );
        end
        
        setappdata( gcbf, 'standards', s );
        setappdata( gcbf, 'standards_original', s0 );
        setappdata( gcbf, 'baselineChannels', [] );
        
        % update legend in the channels/standards plot
        updateChanLegend( gcbf );
    end
% end % addBaselineChannels()


% --- Executes on button press in cb_normalizeStandards.
function cb_normalizeStandards_Callback(hObject, eventdata, handles)
% hObject    handle to cb_normalizeStandards (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    standardsNormalized = get( hObject, 'Value' );
    
    % save this preference
    setWithPref( gcbf, 'normalizeStandards', 'standardsNormalized', ...
        standardsNormalized );
    
    if standardsNormalized
        normalizeStandards( gcbf )
    else
        unnormalizeStandards( gcbf )
    end
% end % cb_normalizeStandards_Callback()

% --- Normalizes the pump-probe reference standards. Assumes standards are
% already loaded and displayed, etc.
function normalizeStandards( fig )
    
    % check s0 to see if we're already normalized
    s0 = getappdata( fig, 'standards_original' );
    
    if isempty(s0)
        s = getappdata( fig, 'standards' );

        % save a copy so we can later undo the normalization
        setappdata( fig, 'standards_original', s );

        % iterate through, normalize each to be unit vectors
        for ind = 1:length(s)
            x = s(ind).x(:);
            x = x / sqrt(x'*x);
            s(ind).x = x;

            set( s(ind).lineobj, 'YData', x );
        end

        setappdata( fig, 'standards', s );
    end
% end % normalizeStandards()

function unnormalizeStandards( fig )
    % retrieve the un-normalized standards
    s = getappdata( fig, 'standards_original' );
    
    if ~isempty( s )

        % iterate through, updating the plot
        for ind = 1:length(s)
            set( s(ind).lineobj, 'YData', s(ind).x );
        end

        setappdata( fig, 'standards', s );
        setappdata( fig, 'standards_original', [] );
        
    end
% end % unnormalizeStandards()


% --- Executes on button press in pushbutton20.
function pushbutton20_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton20 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in cb_enablePhasor.
function cb_enablePhasor_Callback(hObject, eventdata, handles)
% hObject    handle to cb_enablePhasor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    x = get(hObject,'Value');
    setWithPref( gcbf, 'lindecomp_enablePhasor', 'enablePhasor', x );
% end % cb_enablePhasor_Callback()


% --- Executes on button press in cb_threshold.
function cb_threshold_Callback(hObject, eventdata, handles)
% hObject    handle to cb_threshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    val = get(hObject,'Value');
    
    setWithPref( gcbf, 'lindecomp_residualThreshEnable', ...
        'residualThreshEnable', val );
% end % cb_threshold_Callback()

function ed_threshold_Callback(hObject, eventdata, handles)
% hObject    handle to ed_threshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    x0 = getappdata( gcbf, 'residualThreshVal' );
    x = str2double(get(hObject,'String'));
    
    if ~isnan(x) && (x >= 0)
        setWithPref( gcbf, 'lindecomp_residualThreshVal', ...
        'residualThreshVal', x0 );
    else
        x = x0;
    end
    
    set( hObject, 'String', num2str( x ) );
% end % ed_threshold_Callback()


% --- Executes on button press in Force_Pos_Proj.
function Force_Pos_Proj_Callback(hObject, eventdata, handles)
% hObject    handle to Force_Pos_Proj (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Force_Pos_Proj


% --- Executes on button press in ApplyMask.
function ApplyMask_Callback(hObject, eventdata, handles)
% hObject    handle to ApplyMask (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
segment_Callback(hObject, eventdata, handles,1)
