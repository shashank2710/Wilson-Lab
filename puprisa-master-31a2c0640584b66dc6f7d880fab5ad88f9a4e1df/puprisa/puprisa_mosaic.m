function puprisa_mosaic( inputFileNames, xCoords, yCoords, outputFileName )
% puprisa_mosaic
% THIS NEEDS TO READ COORDS FROM FILE NAMES; CAN'T TAKE INTO ACCOUNT THE
% SCALE YET-- THAT IS THE JOB OF viewChannelMosaic, SINCE THIS FILE DOES
% NOT ACTUALLY STITCH RAW DATA TOGETHER. IT JUST MAKES A FILE THAT POINTS
% TO ALL THE RAW DATA.
%
% add several stacks together to compose a mosaic, outputting a single 
% .dat file
%
% The arguments passed to this function is a list of input filenames
% followed by the output filename.
%
% The original image slices remain unmodified.
%
% The new stack will have the header of the first stack in 
% the input list. It is assumed that all image stacks in the fileList have
% the same acquisition parameters (same field of view, objective, number of
% slices, etc.); the only difference between these files should be the x
% and y coordinates.
%
% example usage:
%
% Suppose we have a 4x4 mosaic of image stacks, each image having a 
% 1x1 mm  field of view.
%
% puprisa_mosaic({'topLeft.dat',...
%                 'topRight.dat',...
%                 'lowerLeft.dat',...
%                 'lowerRight.dat'},...
%                [0, 1000, 0, 1000], [0, 0, 1000, 1000], ...
%                'mosaic.dat');
%
% Once completed, using puprisa to open the output 'mosaic.dat' file will
% result in the display of a single delay stack covering the entire 2x2 mm
% wide field of view.
%
% Alternatively, the coordinates may be extracted from the file names of
% the individual tile stacks. To do this, the file must contain both x and
% y coordinates in microns, specified by underscore separated fields, using
% 'Pos' and 'Neg' to specifiy positive(+) or negative(-) values. For
% example, the file name 'tile_xPos0400_yNeg0215.dat' specifies that this
% tile is located at (x = +400, y = -215). To extract coordinates from the
% filename, leave the xCoords and yCoords arguments empty (that is, []):
%
% puprisa_mosaic({'tile_xPos0000_yPos0000.dat',...
%                 'tile_xPos1000_yPos0000.dat',...
%                 'tile_xPos0000_yPos1000.dat',...
%                 'tile_xPos0000_yPos1000.dat'},...
%                [], [], 'mosaic.dat' );
%
% Another example: stitch together all delay stacks starting with
%  the name 'HREM_mosaicDelayStacks_zPos2100'
%
% D = dir('HREM_mosaicDelayStacks_zPos2100*.dat');
% fileNames = {D.name};
% puprisa_mosaic(fileNames,[],[],'mosaic_zPos2100.dat');
%
%

    % === check for test condition ========================================
    if( (nargin == 1) && (strcmp(inputFileNames,'-test')) )
        testMe();
        return;
    end

    % === Set-up and parsing input arguments ==============================
    % first make sure coordinates are specified for /each/ file
    nFiles = length(inputFileNames);
    
    if( isempty( xCoords ) && isempty( yCoords ) )
        [xCoords, yCoords] = getCoordsFromFilenames( inputFileNames );
    end
    
    if( length( xCoords ) < length( inputFileNames ) )
        error('Too few x coordinates for the given list of files.');
    elseif( length( xCoords ) < length( inputFileNames ) )
        error('Too many x coordinates for the given list of files.');
    elseif( length( yCoords ) < length( inputFileNames ) )
        error('Too few y coordinates for the given list of files.');
    elseif( length( yCoords ) < length( inputFileNames ) )
        error('Too many y coordinates for the given list of files.');
    end
    
    % === Open the output file ============================================
    % open output file with write permission
    fidOut = fopen(outputFileName, 'w');
    
    % double-check that the file opened correctly
    if( fidOut == -1 )
        error(['Could not open ', outputFileName]);
    end

    % === Copy header from the first input file ===========================
    
    % open the input file with read permission
    fidIn = fopen(inputFileNames{1});   
    if( fidIn == -1 )
        fclose( fidOut );   % close output file before we crash
        error(['Could not open ', inputFileNames{1}]);
    end
    
    % iterate through the whole input file, looking for header lines,
    % which start with '#' or '%'

    lineIn = fgetl(fidIn);       % get an input line
    foundMosaicLine = 0;
    while ischar(lineIn)         % loop while there are still lines to read
        
        % check if this is a header line by examining the first character
        if (lineIn(1) == '#' || lineIn(1) == '%')
            % replace 'mosaic = 0' with 'mosaic = 1
            if strcmp(lineIn,'# mosaic = 0')
                lineIn = '# mosaic = 1';
                foundMosaicLine =1;
            end
            
            % copy header lines to the output
            fprintf(fidOut, '%s \n', lineIn);   
        end
        
        lineIn = fgetl(fidIn);      % get the next line
    end
    
    if ~foundMosaicLine
        fprintf(fidOut, '# mosaic = 1\n');
    end
    
    fclose( fidIn );            % we're done with the 1st file for now
    
    newSliceNumber = 0;    % keep track of which image slice we're on
    
    % === Loop through all input files, copying the slicenumber lines =====
    for ii = 1:length( inputFileNames )
        
        % open the input file with read permission
        fidIn = fopen(inputFileNames{ii});   
        if( fidIn == -1 )
            fclose( fidOut );   % close output file before we crash
            error(['Could not open ', inputFileNames{1}]);
        end
        
        % iterate through the whole input file, looking for slicenumber
        % lines which start with 'slicenumber'
        
        oldSliceNumber = 0; % keep track of slicenumber in input file
        lineIn = fgetl(fidIn);
        while ischar(lineIn)
            % check if this is a slicenumber line
            if ( strncmp( lineIn, 'slicenumber', 11 ) )
                
                % increment the slice number
                newSliceNumber = newSliceNumber + 1;
                oldSliceNumber = oldSliceNumber + 1;
                
                % parse z position from the slice line
                r = sscanf(lineIn,...
                    'slicenumber = %d , position = %f\t%f\t%f');
                zCoord = r(4);
                          
                % get the slice name to append to the slicenumber line
                % this way puprisa knows where to find the file that
                % contains imaging data for this slice.
                sliceName = inputFileNames{ii};     % start with input .dat
                sliceName = sliceName(1:(end-4));   % cut off '.dat'
                sliceName = [sliceName, '_', num2str(oldSliceNumber)];
                
                % add the slice filename to the end of the line
                lineOut = sprintf(...
                    'slicenumber = %d , position = %f\t%f\t%f\t%s',...
                    newSliceNumber, xCoords(ii), yCoords(ii), zCoord, sliceName );
                
                % then write this slice line to the output
                fprintf(fidOut, '%s \n', lineOut); 
                
            end

            lineIn = fgetl(fidIn);      % get the next line
        end

        fclose( fidIn );            % we're done with the 1st file for now

    end

    % close the output file
    fclose( fidOut );
end

function [xCoords, yCoords] = getCoordsFromFilenames( fileNames )
    % Parse list of filenames to extract x- and y- coordinates. See
    % comments at top of this file (puprisa_mosaic.m) for filename
    % specification.
    
    % convert single file to cell array, if necessary
    if ischar( fileNames )
        fileNames = {fileNames};
    end
    
    nFiles = length(fileNames);
    
    xCoords = zeros(1, nFiles);
    yCoords = zeros(1, nFiles);
    
    for iFile = 1:nFiles
        parsed = regexp(fileNames{iFile},...
            '_x(?<x>(Pos|Neg)\d+)_y(?<y>(Pos|Neg)\d+).dat','names','ignorecase');
        
        xCoords(iFile) = ...
            str2num(regexprep(parsed.x,{'Pos','Neg'},{'+','-'},...
                'ignorecase'));
            
        yCoords(iFile) = ...
            str2num(regexprep(parsed.y,{'Pos','Neg'},{'+','-'},...
                'ignorecase'));
    end
end

function test_getCoordsFromFilenames()
    % generate a few examples, and expected results
    disp('Testing puprisa_mosaic:getCoordsFromFilenames()...');
    pass = 1;
    
    % test a single file
    fprintf('  Testing a single filename...')
    fileNames = 'tile_xPos0400_yNeg0280.dat';
    xCoords_expected = 400;
    yCoords_expected = -280;
    [xCoords, yCoords] = getCoordsFromFilenames( fileNames );
    if( (xCoords ~= xCoords_expected) || (yCoords ~= yCoords_expected) )
        pass = pass*0;
        fprintf(' FAIL\n');
        fprintf('    Test fileName: %s\n',fileNames);
        fprintf('    Expected coordinates: (x = %f, y = %f)\n', ...
            xCoords_expected, yCoords_expected );
        fprintf('    Returned coordinates: (x = %f, y = %f)\n', ...
            xCoords, yCoords );
    else
        fprintf(' PASS\n');
    end
    
    fprintf('  Testing a multiple filenames...')
    fileNames = {'tile_xPos0400_yNeg0280.dat',...
                 'other_xNeg938_ypos321.dat',...
                 'dar_xpos0312_yNeg414.dat',...
                 'wii_xNeg0312_yNeg4215.dat',...
                 'foo_xPos0000_yNeg0000.dat'};
    xCoords_expected = [400, -938, 312, -312,0];
    yCoords_expected = [-280,321,-414,-4215,0];
    [xCoords, yCoords] = getCoordsFromFilenames( fileNames );
    if( ~isequal(xCoords,xCoords_expected) || ~isequal(yCoords, yCoords_expected) )
        pass = pass*0;
        fprintf(' FAIL\n');
        
        disp('fileNames:');
        fileNames
        
        disp('xCoords_expected:');
        xCoords_expected
        
        disp('xCoords:');
        xCoords
        
        
        disp('yCoords_expected:');
        yCoords_expected
        
        disp('yCoords:');
        yCoords
    else
        fprintf(' PASS\n');
    end
    
    if pass == 1
        fprintf('puprisa_mosaic:getCoordsFromFilenames() ... PASS\n');
    else
        fprintf('puprisa_mosaic:getCoordsFromFilenames() ... FAIL\n');
    end
    
end

%==================== TEST FUNCTIONS ======================================
function testMe()
    % test all the sub-functions in this m-file
    
    disp('Testing puprisa_mosaic.m...');
    
    test_getCoordsFromFilenames();
end
