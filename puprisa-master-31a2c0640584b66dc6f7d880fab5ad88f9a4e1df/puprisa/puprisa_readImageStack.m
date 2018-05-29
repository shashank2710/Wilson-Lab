function [slices, header] = puprisa_readImageStack( fileName, varargin)
%PUPRISA_READIMAGESTACK Reads image stack generated by ScanImage
%   Jesse Wilson (jesse.wilson@duke.edu) 2011
%   (Warren Lab, Duke University)
%   Basic format: 
%   commented lines, #-delimited (or %-delimited) with information (key = value) pairs
%   followed by slicenumber=no., position = 1, 2,3
%
% optional argument hDisplayImages shows each image as it loads
%
% This function caches image stacks as a .mat file, and will load the
% cached file instead of the raw ASCII data if the .mat file exists.
%
% Supported Formats:
%
% .dat files from Image Scanner
%
% TIFF stacks
%  stack info (z position or time delay) stored in each page's name
%
% LSM (Zeiss)

    % for now, ignore the extra info, but come back to this later
    % todo: try/catch, return empty if failed
    % todo: return z-pos or t-delay depending on the type of stack taken.
    %       (Should be able to read this from the header)

% CHANGELOG:
% 3/7/2011:
%   - Added feature to cache loaded files in matlab format to speed up
%     loading times.
% 5/25/2011:
%   - fixed header parsing to accept % as a comment character
% 5/27/2011:
%   - fixed reading of date in header
% 6/01/2011:
%   - renamed to puprisa_readImageStack
%   - retain entire header text in header.fullHeaderText
% 6/13/2011:
%   - display progress bar if no hAxDisplayImages specified
% 8/30/13
%   - added ability to read .lsm files; FER
%     To do:`
%         - Allow for multiple images (e.g, Pump-probe and multiphoton) to
%         be loaded.  For now it only adds Pump-probe
%         - Add time stamp
%
% 7/10/14, MCF:
%   - added flag to read in all channels rather than just one channel
%     in case of a tif file input
%
% 8/6/15 FER:
% - added ability to lead MAT files directly.  Assumes a format of [x,y,t/lam/z]
%  initial intended used is for hyperspectral data. 

    readImageStackVersion = 20140218.1;    % increment this each time something changes

    % check if file exists
    fileID = fopen( fileName );
    
    if fileID == -1
        error(['Could not open file: ', fileName]);
    end
    
    fclose( fileID );
    
    % determine file type and load accordingly
    extension = fileName(end-2:end);
    
    switch( extension )
        case 'dat'
            
            [slices, header] = ...
                ReadImageStack_dat( ...
                    fileName, ...
                    readImageStackVersion );
            
        case 'lsm'
            
            [slices, header] = ...
                ReadImageStack_LSM( ...
                    fileName, ...
                    readImageStackVersion );
                
        case 'tif'

            if (~isempty(varargin)) 
                p = inputParser;
                addParameter(p,'TIFLoadAllChannels',true);
                parse(p,varargin{:})
                loadAllChannels = p.Results.TIFLoadAllChannels;
            else
                loadAllChannels = true;
            end
            
            [slices, header] = ...
                ReadImageStack_TIFF( ...
                    fileName, ...
                    readImageStackVersion, ...
                    loadAllChannels);
        case 'mat'      
               [slices, header] = ...
                ReadImageStack_MAT( ...
                    fileName, ...
                    readImageStackVersion ); 
                
    end
    
end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% ReadImageStack_dat()
%
% Read image stack from .dat file, recorded by old Image Scanner program
%
function [slices, header] = ReadImageStack_dat( fileName, readImageStackVersion )

    % look for matlab cached version, if it exists
    cacheReadSuccess = 0;
    
    cacheFileName = [fileName(1:(end-4)), '.mat'];
    
    if exist( cacheFileName )
        dat = load( cacheFileName );
        
        % check that the cache is up to date
        if( isfield(dat, 'readImageStackVersion') )
            if( dat.readImageStackVersion == readImageStackVersion )
                slices = dat.slices;
                header = dat.header;
                cacheReadSuccess = 1;
            else
                cacheReadSuccess = 0;     % if not up to date, flag to re-open
            end
        else
            cacheReadSuccess = 0;
        end
    end

    
    % if unable to read from cache, then open the file for reading
    if( ~cacheReadSuccess )

        fid = fopen( fileName );

        wb = waitbar(0,'Loading image stack...');
  
        % assume imageStack file has '.dat' appended
        baseFileName = fileName(1:(end-4));

        iSlice = 0;
        header = [];
        header.fullHeaderText = {};
        header.nSlices = 0;

        % loop through and determine no. slices
        tline = fgetl(fid);
        while( ischar(tline) )
            if strcmp(tline(1:11),'slicenumber')
                header.nSlices = header.nSlices + 1;
            end
            tline = fgetl(fid); 
        end

        % loop through file line by line
        frewind( fid );
        tline = fgetl(fid);
        while ischar(tline)
            if (tline(1) == '#' || tline(1) == '%')
                % header line -- save it in the header text
                header.fullHeaderText{end+1} = tline;

                % parse header line into:
                % '# key = value' or '% key = value'
                [~, tok] = ...
                    regexpi(tline,'[#|%] (.+) = (.*)','match','tokens');

                if isempty(tok)
                    % treat as a comment
                    disp(['comment: ', tline]);
                else
                    % treat as a key/value pair

                    key = tok{:}{1};
                    val = tok{:}{2};

                    % convert value to number, if possible
                    if ~(strcmp(key,'date') || strcmp(key,'time'))
                        valNum = str2num(val);
                        if ~isempty( valNum )
                            val = valNum;
                        end
                    end

                    % remove all non-alphanumeric an underscore characters
                    % from key, to make a valid key,
                    % and get rid of first character if it's a number
                    key = regexprep(key,'^\d+|[^\w]','');


                    % store in header structure
                    header = setfield(header, key, val);
                end

            else
                % check if header was formed correctly
                % just check presence of one field
                % (ideally we should check all fields)
                if ~isfield(header, 'linesperframe')
                    error('Header not read correctly.');
                end

                % move on to loading the slice
                iSlice = iSlice + 1;
                % break apart into position values
                A = sscanf(tline, 'slicenumber = %d , position = %f\t%f\t%f\t%s');

                slices(iSlice).slicenumber = A(1);
                slices(iSlice).posX = A(2);
                slices(iSlice).posY = A(3);
                slices(iSlice).posZ = A(4);
                slices(iSlice).delays = -A(4) * 2 / 299.792458;

                % generate file name corresponding to slice
                if length(A) == 4
                    imageFileName = ...
                        [baseFileName,'_',num2str(slices(iSlice).slicenumber)];
                elseif length(A) > 4
                    thePath = baseFileName( ...
                        1: max(strfind( baseFileName,filesep)));
                    imageFileName = [thePath, filesep, char(A(5:end).')];
                else
                    error(['slices line of file does not match format:',...
                         tline]);
                end

                slices(iSlice).imageData =...
                    readImage(imageFileName, ...
                              header.linesperframe, ...
                              header.pixelsperline );

                waitbar(iSlice / header.nSlices,wb)

            end
            tline = fgetl(fid);
        end
        
       
        fclose(fid);        % close the file
    
        close(wb);          % close the waitbar
    end
     
    % save to cache
    try
        save( cacheFileName, 'slices','header', 'readImageStackVersion' );
    catch
        warning('Could not write cache file. Remember that archive is write-protected.');
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% ReadImageStack_LSM()
%
% Read image stack from a Zeiss .lsm file
%
function [slices, header] = ReadImageStack_LSM( fileName, readImageStackVersion )
    % look for matlab cached version, if it exists
    cacheReadSuccess = 0;
    
    cacheFileName = [fileName(1:(end-4)), '.mat'];
    
    if exist( cacheFileName )
        dat = load( cacheFileName );
        
        % check that the cache is up to date
        if( isfield(dat, 'readImageStackVersion') )
            if( dat.readImageStackVersion == readImageStackVersion )
                slices = dat.slices;
                header = dat.header;
                cacheReadSuccess = 1;
            else
                cacheReadSuccess = 0;     % if not up to date, flag to re-open
            end
        else
            cacheReadSuccess = 0;
        end
    end

    
    % if unable to read from cache, then open the file for reading
    if( ~cacheReadSuccess )
        ImageStack  = tiffread32(fileName);

        header = [];
        header.fullHeaderText = {};


        header.nSlices = ImageStack(1).lsm.DimensionTime;
        header.filename = ImageStack(1).file_name;
        header.acqusitionnumber = 1;
        header.version = 1;
        header.pixelsperline = ImageStack(1).lsm.DimensionX;
        header.linesperframe = ImageStack(1).lsm.DimensionY;
        header.numofchannels = ImageStack(1).lsm.DimensionChannels;
        header.scanrangex = ImageStack(1).lsm.VoxelSizeX*ImageStack(1).lsm.DimensionX*1E3;
        header.scanrangey = ImageStack(1).lsm.VoxelSizeY*ImageStack(1).lsm.DimensionY*1E3;
        header.msline = ImageStack(1).lsm.TimeInterval/ImageStack(1).lsm.DimensionY*1E3;
        header.variableaxisZt = ImageStack(1).lsm.DimensionZ;
        header.date = fileName(end-11:(end-4));
        header.time = '00:00:00'; % get time
        header.usefile = 1;

        Time_Stack = 1; 
        if header.variableaxisZt> header.nSlices;
            header.nSlices = header.variableaxisZt;
            Time_Stack = 0;
        end

        header.fullHeaderText{1} = ['Num time pts = ',num2str(ImageStack(1).lsm.DimensionTime)];
        header.fullHeaderText{2} = ['file name:' ,ImageStack(1).file_name];
        header.fullHeaderText{3} = 1;
        header.fullHeaderText{4} = 1;
        header.fullHeaderText{5} = ['num x pixels = ',num2str(ImageStack(1).lsm.DimensionX)];
        header.fullHeaderText{6} = ['num y pixels = ',num2str(ImageStack(1).lsm.DimensionY)];
        header.fullHeaderText{7} = ['num z pixels= ',num2str(ImageStack(1).lsm.DimensionChannels)];
        header.fullHeaderText{8} = ['x range mm = ', num2str(ImageStack(1).lsm.VoxelSizeX*ImageStack(1).lsm.DimensionX*1E3)];
        header.fullHeaderText{9} = ['y range mm = ',num2str(ImageStack(1).lsm.VoxelSizeY*ImageStack(1).lsm.DimensionY*1E3)];
        header.fullHeaderText{11} = 1;
        header.fullHeaderText{10} = ['line scan time = ',num2str(ImageStack(1).lsm.TimeInterval/ImageStack(1).lsm.DimensionY*1E3)];
        header.fullHeaderText{12} = ['date: ',fileName(end-11:(end-4))];
        header.fullHeaderText{13} = ['time: ','00:00:00']; % get time
        header.fullHeaderText{14} = 1;

        if Time_Stack


            display('-------------------------------------------');
            display('Find file with time axis or do sequential?');
            display('Sequential : 0');
            display('Find file: 1');
            FindFile = input('Enter 0 or 1: ');
            display('-------------------------------------------');


            if FindFile
    %                         wd = pwd();
                wd = getpref('puprisa','imageWorkingDirectory');

                [DelayFileName,DelayFilePath] = uigetfile([wd,'/*.dat'],...
                    'Select corresponding delay file' );



                TimeDelays_ps = load([DelayFilePath,DelayFileName]);
            else
                TimeDelays_ps = 1:header.nSlices;
            end
        else
            TimeDelays_ps(1:header.nSlices) = 0; 
        end 

        if (length(TimeDelays_ps)~= header.nSlices);
            error('Incorrect time delay file: num of time pts does not match num of slices')
        end

        wb = waitbar(0,'Loading image stack...');

        for iSlice = 1: header.nSlices
            slices(iSlice).slicenumber = iSlice;

           if (length(ImageStack(1).data)==1) ; 
            slices(iSlice).imageData{1} = double(ImageStack(iSlice).data{1});
           else
               for iData = 1:length(ImageStack(1).data);
                slices(iSlice).imageData{iData} = double(ImageStack(iSlice).data{iData});
               end
           end


            % need to change for multiple types of images acquiered
            slices(iSlice).posX = ImageStack(iSlice).lsm.OriginX;
            slices(iSlice).posY = ImageStack(iSlice).lsm.OriginY;
           if Time_Stack; 
            slices(iSlice).posZ = ImageStack(iSlice).lsm.OriginZ;
           else
            slices(iSlice).posZ = iSlice;
           end
            slices(iSlice).delays = TimeDelays_ps(iSlice);
            waitbar(iSlice / header.nSlices,wb)
        end

         close( wb );       % close the waitbar
    end
    
    
    % save to cache
    save( cacheFileName, 'slices','header', 'readImageStackVersion' );
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% ReadImageStack_TIFF()
%
% Read image stack from a TIF file
%
function [slices, header] = ReadImageStack_TIFF( fileName, readImageStackVersion, loadAllChannels )
    header = [];

    if (loadAllChannels)
        % figure out how many channels there are
        underScores = strfind(fileName, '_');
        nameBase = fileName(1:underScores(end));
        if exist( [nameBase 'CH1.tif'], 'file')
            header.numofchannels = 1;
            fileList{1} = [nameBase 'CH1.tif'];
        else
            error(['Could not find first channel file: ', nameBase 'CH1.tif']);
        end
        if exist( [nameBase 'CH2.tif'], 'file')
            header.numofchannels = 2;
            fileList{2} = [nameBase 'CH2.tif'];
            if exist( [nameBase 'CH3.tif'], 'file')
                header.numofchannels = 3;
                fileList{3} = [nameBase 'CH3.tif'];
                if exist( [nameBase 'CH4.tif'], 'file')
                    header.numofchannels = 4;
                    fileList{4} = [nameBase 'CH4.tif'];
                end
            end
        end
    else % if(loadAllChannels)
        header.numofchannels = 1;
        fileList{1} = fileName;
    end % if(loadAllChannels)
    
    info  = imfinfo(fileName);
    
    % fill in header
    header.nSlices = length( info );
    header.linesperframe = info(1).Width;
    header.fullHeaderText = info(1).ImageDescription;
    
    % determine whether it's a t- or z-stack
    parsed = sscanf(info(1).PageName, '%c = %f %c%c');
    stackType = char(parsed(1));
    switch stackType
        case 't'
            header.scanaxis = 0;
        case 'z'
            header.scanaxis = 1;
        otherwise
            error(['Invalid page name: ', info(1).PageName]);
    end
    
    
    wb = waitbar(0,'Loading image stack...');

    
    % temporarily turn off warning for missing 'TIFF Photometric Tag'
    warning('off', 'MATLAB:imagesci:rtifc:missingPhotometricTag');
    
    for iSlice = 1:header.nSlices

        for iChannel = 1:header.numofchannels
            slices(iSlice).imageData{iChannel} = imread( fileList{iChannel}, ...
                                                         'index', iSlice );
        end
        
        % parse slice info
        parsed = sscanf(info(iSlice).PageName, '%c = %f %c%c');
        
        if header.scanaxis
            % z stack
            slices(iSlice).posZ = parsed(2);
            slices(iSlice).delays = 2*parsed(2) / 299.792458;
        else
            slices(iSlice).delays = parsed(2);
            slices(iSlice).posZ = 0.5 * parsed(2) * 299.792458;
            
        end
        
        % leave other positions empty; fill these in later
        slices(iSlice).posX = 0;
        slices(iSlice).posY = 0;
        
        waitbar(iSlice / header.nSlices,wb)
    end
    
    % turn warning for missing 'TIFF Photometric Tag' back on
    warning('on', 'MATLAB:imagesci:rtifc:missingPhotometricTag');

    
    close(wb);
    
    
    

end



function [slices, header] = ReadImageStack_MAT( fileName, readImageStackVersion )  
% for a matlab file that has a data cube. 
% No need to specify time axis
   
        dat = load(fileName);
        contents = fieldnames(dat);
        
        slices = [];
        header = [];
        header.fullHeaderText = {};
        
        if length(contents)==1;
            % if only one matrix,  it is the data cube
            % third axis used as time axis (slices)
            NoSlices = eval(['size(dat.',contents{1},',3)']);
            for iSlice = 1:NoSlices;
                slices(iSlice).slicenumber = iSlice;
                slices(iSlice).imageData{1} = eval(['dat.',contents{1},'(:,:,',num2str(iSlice),')']);
                slices(iSlice).delays = iSlices;
            end
        else     
            DataCubes = 0;
             for ContentNo = 1:length(contents);
                    fieldSize = eval(['size(dat.',contents{ContentNo},')']);
                    
                        if length(fieldSize)<3;
                            for iSlice = 1:max(fieldSize); 
                                 slices(iSlice).slicenumber = iSlice;
                                 slices(iSlice).delays = ...
                                     eval(['dat.',contents{ContentNo},'(',num2str(iSlice),')']);
                            end
                        else
                            DataCubes = DataCubes + 1;
                          for iSlice = 1:fieldSize(3);  
                            slices(iSlice).imageData{DataCubes} = ...
                                eval(['dat.',contents{ContentNo},'(:,:,',num2str(iSlice),')']); 
                          end
                        end
            
             end
        end
        
        header.nSlices = length(slices);
        header.filename = fileName;
        header.acqusitionnumber = 1;
        header.version = 1;
        header.pixelsperline = size(slices(1).imageData{1},2);
        header.linesperframe = size(slices(1).imageData{1},1);
        header.numofchannels = 1;
        header.scanrangex = 500;
        header.scanrangey = 500;
        header.msline = 1;
        header.variableaxisZt = 1;
        header.date = '0000-00-00';
        header.time = '00:00:00'; % get time
        header.usefile = 1;
        header.fullHeaderText{1} = 'Matlab matrix loaded directly';
        

        slices(iSlice).posX = 1;
        slices(iSlice).posY = 1;
        slices(iSlice).posZ = 1;

  
end



