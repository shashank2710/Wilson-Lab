function varargout = puppy(varargin)
% PUPPY MATLAB code for puppy.fig
%      PUPPY, by itself, creates a new PUPPY or raises the existing
%      singleton*.
%
%      H = PUPPY returns the handle to a new PUPPY or the handle to
%      the existing singleton*.
%
%      PUPPY('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PUPPY.M with the given input arguments.
%
%      PUPPY('Property','Value',...) creates a new PUPPY or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before puppy_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to puppy_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help puppy

% Last Modified by GUIDE v2.5 08-Oct-2012 16:34:03

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @puppy_OpeningFcn, ...
                   'gui_OutputFcn',  @puppy_OutputFcn, ...
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


% --- Executes just before puppy is made visible.
function puppy_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to puppy (see VARARGIN)

% Choose default command line output for puppy
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes puppy wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = puppy_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --------------------------------------------------------------------
function mnu_file_Callback(hObject, eventdata, handles)
% hObject    handle to mnu_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function mnu_help_Callback(hObject, eventdata, handles)
% hObject    handle to mnu_help (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function mnu_about_Callback(hObject, eventdata, handles)
% hObject    handle to mnu_about (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    scrsz = get(0,'ScreenSize');
    
    
    % load fun image
    I = imread(['..',filesep,'imagesIconsAndLogos',filesep,'puppy.png']);
    [imht,imwd,~] = size(I);
    
    txtht = 80;
    
    fwid = 340;
    fht = imht+txtht;
    
    f=figure('WindowStyle','modal','Units','pixels',...
        'Position',[scrsz(3)/2 - fwid/2,scrsz(4)/2-fht/2, fwid, fht],...
        'Name','About puppy','NumberTitle','off');
    imshow(I);
    set(gca,'units','pixels','Position',[fwid/2-imwd/2,txtht,imwd,imht]);
    
    
    % LEFT OFF HERE
    uicontrol('Style','Edit','String',['Puppy: fetching quantitative metrics from piles of ',...
            'pump-probe image stacks. Created 2012 by Jesse W. Wilson ',...
            '(syrex314@gmail.com), Warren Lab, Duke University.'],...
            'Units','Pixels','Position',[0,0,fwid,txtht],'Enable','Inactive',...
            'Max',2,'Min',0,'HorizontalAlignment','left');
        
        
        


% --- Executes on button press in pb_loadImage.
function pb_loadImage_Callback(hObject, eventdata, handles)
% hObject    handle to pb_loadImage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    % recall last accessed directory
    if ispref('puppy','imageStackPath')
        pathName = getpref('puppy','imageStackPath');
    else
        pathName = pwd();
        addpref('puppy','imageStackPath',pathName);
    end

    % prompt user for imageStack file to be loaded
    [fileName, pathName] = uigetfile( [pathName,filesep,'*.dat'], ...
                                      'Pick an image stack .dat file',...
                                      'multiselect','on');
    if isequal(fileName, 0); return; end

    % strip trailing filesep from pathName
    pathName = pathName(1:(end-1));
    
    setpref('puppy','imageStackPath',pathName);
    
    addFiles( gcbf, fileName, pathName, handles );

function addFiles( fig, fileName, pathName, handles ) 

    
    tableData = get(handles.tbl_main,'Data');
    [nr,~] = size(tableData);
    
    if iscell(fileName)
        for ii = 1:length(fileName)


            % add to list
            tableData{nr+ii,1} = '';
            tableData{nr+ii,2} = pathName;
            tableData{nr+ii,3} = fileName{ii};
            tableData{nr+ii,4} = '';

        end

    else
        % add to list
        tableData{nr+1,1} = '';
        tableData{nr+1,2} = pathName;
        tableData{nr+1,3} = fileName;
        tableData{nr+1,4} = '';
        
        set(fig,'pointer','watch');
        drawnow;
        
        loadImage(fig, [pathName,filesep,fileName]);
        
        set(fig,'pointer','arrow');
    end


        cw = get( handles.tbl_main,'ColumnWidth');
        set(handles.tbl_main,'Data',tableData);
        set(handles.tbl_main,'ColumnWidth',cw);
        

function loadImage( fig, fileName )
    handles = guihandles(fig);

    [S,header] = puprisa_readImageStack( fileName );
    [t,X] = puprisa_getChannelFromSlices(S,1,'');

    setappdata(fig,'fullImageStack',S);
    setappdata(fig,'header',header)
    setappdata(fig,'imageStack',X);
    setappdata(fig,'imageStack_original',X);
    setappdata(fig,'delays',t);

    %preprocessImageStack( fig );
    im = get(handles.ax_preview,'UserData');
    C = sum(X,3);
    [C,I] = max(abs(X),[],3);    % max intensity projection
    lim = max(abs(X(:)));

    [nr,nc,nd] = size(X);
    for ir = 1:nr
        for ic = 1:nc
            C(ir,ic) = X(ir,ic,I(ir,ic));
        end
    end
    
    s = std(C(:));
    
    x = 1:nc;
    y = 1:nr;
    
    lim = 5*s;
    
    
    set(im,'XData',x,'YData',y,'CData',C);
    set(handles.ax_preview,'clim',[-1,1]*lim);
    set(handles.ax_preview,'xlim',[1,nc]);
    set(handles.ax_preview,'ylim',[1,nr]);
    
    drawnow;

% --- Executes on button press in pb_loadDir.
function pb_loadDir_Callback(hObject, eventdata, handles)
% hObject    handle to pb_loadDir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% recall last accessed directory
    if ispref('puppy','imageStackPath')
        pathName = getpref('puppy','imageStackPath');
    else
        pathName = pwd();
        addpref('puppy','imageStackPath',pathName);
    end

    % bring up uigetfile to pick a directory
    pathName = uigetdir(pathName);

    % abort if user clicked 'cancel'
    if isequal( pathName, 0 ); return; end;

    % remember this path
    setpref('puppy','imageStackPath',pathName);

    % then drill down through directories within,
    fileNames = recurs_findImageStacks( pathName );

    addFiles( gcbf, fileNames, '', handles );

function fileNames = recurs_findImageStacks(pathName)
    fileNames = {};
    % locate all .dat files (assume only delay stacks for now)
    dd = dir([pathName,filesep,'*.dat']);
    for ii = 1 : length(dd)
        disp(['===> ', dd(ii).name]);
        fileNames{length(fileNames)+1} = [pathName,filesep,dd(ii).name];
    end

    % recurse through subdirectories
    dd = dir(pathName);
    dd_dirs = dd(find([dd.isdir]));

    for id = 1 : length(dd_dirs)
        s = dd_dirs(id).name;
        if ~strcmp(s,'.') && ~strcmp(s,'..')
            disp(dd_dirs(id).name);
            fileNamesThisDir = recurs_findImageStacks([pathName,filesep,dd_dirs(id).name]);

            fileNames = {fileNames{:}, fileNamesThisDir{:}};
        end
    end

% --- Executes on button press in pb_remove.
function pb_remove_Callback(hObject, eventdata, handles)
% hObject    handle to pb_remove (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    
    % locate selected line of table data, and remove it
    td = get(handles.tbl_main,'Data');
    [nr,~] = size(td);
    ii = getappdata(gcbf,'iSelected');


    % remove element ii from s
    if ii == nr
        td = td(1:(ii-1),:);    % remove the last element
        ii = ii-1;          % and decrement the index to the new 'last' element
    elseif ii == 1      
        td = td(2:end,:);       % remove the first element
                            % leave the index at 1
    else
        td = [td(1:(ii-1),:); td((ii+1):end,:)];   % cut out the middle element ii
                                            % leave the index at ii
    end

    % store changes in the listbox
    
    cw = get( handles.tbl_main,'ColumnWidth');
    set(handles.tbl_main,'Data',td);
    set(handles.tbl_main,'ColumnWidth',cw);

    


% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pb_go.
function pb_go_Callback(hObject, eventdata, handles)
% hObject    handle to pb_go (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    td = get(handles.tbl_main,'Data');
    [nr,~] = size(td);
    
    runBatchScript( 1:nr, handles );

function runBatchScript( whichRows, handles )


    td = get(handles.tbl_main,'Data');
    [nr,~] = size(td);
    
    updateProgress( 0, handles );
    
    
    % get custom matlab code
     puppyScriptDir = getPuppyScriptDir();
     
     dataDir = getappdata(gcbf,'rootDataDir');
     outDir = getappdata(gcbf,'outDir')
     
     puppyScriptFn = getappdata(gcf,'selectedScript');
     puppyScriptFn = puppyScriptFn(1:(end-2));  % get rid of the extra .m in the filename
     
     wd0 = pwd();
     cd(puppyScriptDir);
    
    for ii = 1:length(whichRows)
        ir = whichRows(ii);
        
        % construct filename
        if ~isempty( td{ir,2} )
            fullFileName = [dataDir,filesep,td{ir,2}, filesep, td{ir,3}];
        else
            fullFileName = td{ir,3};
        end
       
        disp(['Beginning batch processing on entry ', num2str(ir)]);
        disp(fullFileName);
        
        
        loadImage( gcbf, fullFileName );
        
        S = getappdata( gcbf, 'fullImageStack' );
        header = getappdata( gcbf, 'header' );
        
        % run custom code
        cd(puppyScriptDir);
        out = feval(puppyScriptFn, S, header);
        

        
        % if the first time around, find out names of all outputs
        % and label columns
        if ir == 1
            outputNames = fieldnames( out );
            tableHeader = get( handles.tbl_main,'columnName' );
            newCols = outputNames( ~ismember( outputNames, tableHeader ) );
            
            %tableHeader = unique( [tableHeader; outputNames], 'first' );
            tableHeader = [tableHeader; newCols];
            
            set( handles.tbl_main, 'ColumnName', tableHeader );
        end
        
        % store everything in an array
        batchOut(ir) = out;
        
        % store into the table
        for iField = 1 : length( outputNames )
            outName = outputNames{iField};
            iCol = find(cellfun(@any,strfind(tableHeader,outName)));
            
            
            outVal = getfield(out,outName);
            
            % store in table if it's a single number
            if( length(outVal) == 1 )
                td{ir, iCol} ...
                    = getfield(out,outName);
            end
            
            % if it looks like an image, write to output directory
            if( isa( outVal , 'uint8' ) )
                % TODO: check if it is a 2D array
                outFileName = [outDir,filesep,td{ir,3}(1:(end-4)) ,'_',outName,'.png'];
                
                imwrite(outVal, outFileName, 'PNG');
                
            end
        end
        
        set(handles.tbl_main,'Data',td);
        setappdata(gcbf,'batchOut',batchOut);

        updateProgress( ir/nr, handles );
    end
    
    set(handles.tbl_main,'Data',td);
    
    % Group results and box plot
    classes = unique(td(:,4));
    
    % SORT RESULTS INTO APPROPRIATE BINS FOR BOXPLOT
    axes(handles.ax_scatterplot);
    V = zeros(nr, length(classes));
    V = V + NaN;
    for iClass = 1:length(classes)
        % find all entries for this class
        for ir = 1:nr
            if strcmp(td{ir, 4}, classes{iClass})
                V(ir,iClass) = td{ir,5};
            end
        end
        
    end
    
    axpos = get(handles.ax_scatterplot,'Position');
    boxplot(handles.ax_scatterplot,V,'labels',classes);
    set(handles.ax_scatterplot,'Position',axpos);
    
    for ii = 1:length(classes)
        v = V(:,ii);
        vm = 1-isnan(v);
        vv = v(vm==1);
        line(zeros(1,length(vv))+ii+0.2+rand(1,length(vv))*0.03,vv,'marker','.','linestyle','none');
    end
        
    % restore working directory
    cd( wd0 );

% --- Executes on button press in pushbutton6.
function pushbutton6_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton7.
function pushbutton7_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function ed_matlabCode_Callback(hObject, eventdata, handles)
% hObject    handle to ed_matlabCode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ed_matlabCode as text
%        str2double(get(hObject,'String')) returns contents of ed_matlabCode as a double


% --- Executes during object creation, after setting all properties.
function ed_matlabCode_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ed_matlabCode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pb_test.
function pb_test_Callback(hObject, eventdata, handles)
% hObject    handle to pb_test (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    
    runBatchScript( 1, handles );


% --- Executes on button press in pushbutton10.
function pushbutton10_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes during object creation, after setting all properties.
function tbl_main_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tbl_main (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
td = get(hObject,'Data');
td = {};
set(hObject,'Data',td);


% --- Executes when selected cell(s) is changed in tbl_main.
function tbl_main_CellSelectionCallback(hObject, eventdata, handles)
% hObject    handle to tbl_main (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Indices: row and column indices of the cell(s) currently selecteds
% handles    structure with handles and user data (see GUIDATA)
if length(eventdata.Indices) > 0
    setappdata(gcbf,'iSelected',eventdata.Indices(1));
else
    setappdata(gcbf,'iSelected',0);
end


% --- Executes during object creation, after setting all properties.
function ax_preview_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ax_preview (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate ax_preview
    im = imagesc(0);
    puprisa_colorMap('pumpProbe');
    set(hObject,'UserData',im);
    set(hObject,'Tag','ax_preview');


% --- Executes during object creation, after setting all properties.
function ax_progress_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ax_progress (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate ax_progress

    % this progress bar operates via a patch object

    p = patch([0,0,1,1,0],[0,1,1,0,0],'b');
    set(hObject,'UserData',p);
    
function updateProgress( fraction, handles )
    % update progress bar given a new fraction and the handles object
    % (progress is on scale of 0 to 1)
    p = get(handles.ax_progress, 'UserData');
    
    x = [0,0,fraction,fraction,0];
    set(p,'XData', x);
    drawnow;
    
    
    


% --------------------------------------------------------------------
function Untitled_4_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function mnu_file_newTable_Callback(hObject, eventdata, handles)
% hObject    handle to mnu_file_newTable (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function mnu_file_openTable_Callback(hObject, eventdata, handles)
% hObject    handle to mnu_file_openTable (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    % prompt user for a file to read
    if ~ispref( 'puppy', 'dataTablePath' )
        addpref('puppy', 'dataTablePath', [pwd(),filesep] );
    end
    pathName = getpref('puppy','dataTablePath');
    
    [fileName, pathName] = uigetfile( '*.csv', 'Open data table',[pathName,filesep,'*.csv'] );
    
    if isequal( fileName, 0 )
        return;
    end
    
    fileName = [pathName,fileName];

    
    setpref('puppy', 'dataTablePath', pathName );
        
    file = fopen(fileName,'r');
    
    % make sure it opened OK
    if( isempty( file ) )
        errorDlg('Unable to open file.');
        error(['Unable to open file: ', file]);
    end
    
    % read first row as the header
    tline = fgetl( file );
    headerColumns = textscan(tline,'%s','Delimiter',',');
    headerColumns = headerColumns{1};
    nCol = length( headerColumns );
    
    set( handles.tbl_main, 'ColumnName', headerColumns );
    
    % read in each row as data
    iRow = 1;
    tline = fgetl( file )

    while ~isequal( tline, -1 )
        
        columns = textscan(tline,'%s','Delimiter',',')
        
        for iCol = 1 : nCol
            % determine if it's text or numeric
            s = columns{1}{iCol}
            k = str2double( s )
            if isnan( k )
                % store as a string
                data{iRow, iCol} = s;
            else
                % store as number
                data{iRow, iCol} = k;
            end
        end
        

        iRow = iRow + 1;

        tline = fgetl( file )
    end
    
    set( handles.tbl_main, 'Data', data );

    
    fclose(file);
    
    setActiveFileName( fileName );
    


% --------------------------------------------------------------------
function mnu_file_saveTable_Callback(hObject, eventdata, handles)
% hObject    handle to mnu_file_saveTable (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    
    % since we will save in a comma-separated value (CSV) format,
    % commas are not allowed in any of the data entries; search for commas
    % and alert the user
    data = get( handles.tbl_main, 'Data' );
    
    if( ~ all(all(cellfun('isempty',strfind(data(cellfun(@ischar,data)),','))) ))
        errordlg( 'Data table contains commas; these are not allowed when saving to .csv format.' );
        return;
    end
    
    header = get( handles.tbl_main, 'ColumnName' );
    
    if( ~ all(cellfun('isempty',strfind(header,','))) )
        errordlg( 'Table header contains commas; these are not allowed when saving to .csv format.' );
        return;
    end

    fileName = getappdata(gcbf,'activeFileName' );
    saveTable( fileName, handles )


% --------------------------------------------------------------------
function mnu_file_saveTableAs_Callback(hObject, eventdata, handles)
% hObject    handle to mnu_file_saveTableAs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    % prompt user for a file to write
    if ~ispref( 'puppy', 'dataTablePath' )
        addpref('puppy', 'dataTablePath', [pwd(),filesep] );
    end
    pathName = getpref('puppy','dataTablePath');
    
    [fileName, pathName] = uiputfile( '*.csv', 'Save data table as:',[pathName,filesep,'Untitled.csv'] );
    fileName = [pathName,fileName];
    
    if isequal( fileName, 0 )
        return;
    end
    
    setpref('puppy', 'dataTablePath', pathName );
    
    saveTable( fileName, handles);
    
    
function saveTable( fileName, handles )

    % read table data and header
    data = get( handles.tbl_main, 'Data' );
    header = get( handles.tbl_main, 'ColumnName' );
    
    
    
    
    % figure out column data types
    [nRow, nCol] = size( data );
    
    
    % open a file
    file = fopen( fileName, 'w' );
    
    % write header
    for iCol = 1 : nCol
        fprintf( file, '%s', header{iCol} );
        if iCol < nCol
            fprintf( file, ',' );
        end
    end
    fprintf(file,'\n');
    
    % iterate through and write each row
    for iRow = 1 : nRow
        for iCol = 1 : nCol
            % get data
            d = data{iRow, iCol};
            
            % determine data type and write to file
            if( isnumeric( d ) )
                fprintf( file, '%f', d );
            else
                fprintf( file, '%s', d ); 
            end
            
            
            if iCol < nCol
                fprintf( file, ',' );
            end
        end
        fprintf(file,'\n');
        
    end
    
    % All done. Close the file.
    fclose( file );
    
    setActiveFileName( fileName );


% --------------------------------------------------------------------
function mnu_file_close_Callback(hObject, eventdata, handles)
% hObject    handle to mnu_file_close (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

function setActiveFileName( fileName )
% remember the file name we're working with, so users can save without
% having to specify the filename again
    setappdata(gcbf,'activeFileName', fileName );
    
    % enable the 'save' menu uption
    handles = guihandles(gcbf);
    set(handles.mnu_file_saveTable,'Enable','on');



function ed_rootDataDir_Callback(hObject, eventdata, handles)
% hObject    handle to ed_rootDataDir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ed_rootDataDir as text
%        str2double(get(hObject,'String')) returns contents of ed_rootDataDir as a double


% --- Executes during object creation, after setting all properties.
function ed_rootDataDir_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ed_rootDataDir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end

    setappdata(get(hObject,'parent'), 'rootDataDir', get(hObject,'string'));


% --- Executes on selection change in pmnu_batchScript.
function pmnu_batchScript_Callback(hObject, eventdata, handles)
% hObject    handle to pmnu_batchScript (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns pmnu_batchScript contents as cell array
%        contents{get(hObject,'Value')} returns selected item from pmnu_batchScript
    contents = cellstr(get(hObject,'String'));
    selectedScript = contents{get(hObject,'Value')};
    setappdata(gcf,'selectedScript', selectedScript);

% --- Executes during object creation, after setting all properties.
function pmnu_batchScript_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pmnu_batchScript (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
    % set color
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end
    
    % scan scripts directory and populate this popup menu
    puppyScriptDir = getPuppyScriptDir();
    
    d = dir([puppyScriptDir,filesep,'*.m']);
    
    scriptsList = {d.name};
    
    set(hObject,'String',scriptsList);
    setappdata(gcf,'selectedScript', scriptsList{1});
    
function puppyScriptDir = getPuppyScriptDir()
    % locate batch processing scripts directory
    puppyPath = which('puppy');
    puppyScriptDir = [puppyPath(1:strfind(puppyPath,[filesep,'puppy.m']))];
    puppyScriptDir= [puppyScriptDir,'puppy_scripts'];



function ed_outDir_Callback(hObject, eventdata, handles)
% hObject    handle to ed_outDir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ed_outDir as text
%        str2double(get(hObject,'String')) returns contents of ed_outDir as a double
setappdata(ancestor(hObject,'figure'), 'outDir', get(hObject,'string'));
% --- Executes during object creation, after setting all properties.
function ed_outDir_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ed_outDir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
    setappdata(ancestor(hObject,'figure'), 'outDir', get(hObject,'string'));


% --- Executes on button press in pushbutton14.
function pushbutton14_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function edit5_Callback(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit5 as text
%        str2double(get(hObject,'String')) returns contents of edit5 as a double


% --- Executes during object creation, after setting all properties.
function edit5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit6_Callback(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit6 as text
%        str2double(get(hObject,'String')) returns contents of edit6 as a double


% --- Executes during object creation, after setting all properties.
function edit6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
