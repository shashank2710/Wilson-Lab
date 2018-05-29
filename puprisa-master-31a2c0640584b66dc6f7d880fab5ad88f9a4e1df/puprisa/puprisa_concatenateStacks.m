function puprisa_concatenateStacks( varargin )
% puprisa_concatenateStacks
% add several stacks together, outputting a single .dat file
%
% The arguments passed to this function is a list of input filenames
% followed by the output filename.
%
% example usage:
% puprisa_concatenateStacks('lowerStack.dat',...
%                           'upperStack.dat',...
%                           'concatenated.dat');
%
% This will take two stacks, lower and upper, and combine them into
% a single stack, saving the new combined stack as 'concatenated.dat'
%
% This may be done for an arbitrary number of stacks, e.g.:
%
% puprisa_concatenateStacks('stack1.dat',...
%                           'stack2.dat',...
%                           'stack3.dat',...
%                           'stack4.dat',...
%                           'manyPutTogether.dat');
% 
% The original image slices are kept unmodified.
%
% The new stack will have the header of the first stack in 
% the input list.

    % === Set-up and parsing input arguments ==============================
    % first make sure at least 3 filenames have been specified
    if( nargin < 3 )
        error('Must specify at least 2 input files and 1 output file');
    end

    inputFileNames = {varargin{1:(end-1)}};
    outputFileName = varargin{end};
    
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
    while ischar(lineIn)         % loop while there are still lines to read
        
        % check if this is a header line by examining the first character
        if (lineIn(1) == '#' || lineIn(1) == '%') 
            
            % copy header lines to the output
            fprintf(fidOut, '%s \n', lineIn);   
        end
        
        lineIn = fgetl(fidIn);      % get the next line
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
                
                % Replace the slicenumber in this line
                % This replacement is done with regular expressions.
                lineOut = ...
                    regexprep(lineIn,...
                              '(slicenumber = )(\d+)',...
                              ['$1',num2str(newSliceNumber)] );
                          
                % get the slice name to append to the slicenumber line
                % this way puprisa knows where to find the file that
                % contains imaging data for this slice.
                sliceName = inputFileNames{ii};     % start with input .dat
                sliceName = sliceName(1:(end-4));   % cut off '.dat'
                sliceName = [sliceName, '_', num2str(oldSliceNumber)];
                
                % add the slice filename to the end of the line
                lineOut = sprintf('%s\t%s', lineOut, sliceName );
                
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