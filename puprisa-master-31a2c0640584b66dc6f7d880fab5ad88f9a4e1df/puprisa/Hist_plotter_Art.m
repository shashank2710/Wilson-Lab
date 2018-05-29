% Function to process output of 'Multiple files delayConst Histogram'


% prompt user for an image stack to load
clear all;
% recall working directory
%wd = pwd;
wd = 'D:\Art imaging\2011-08-31';
% or use MATLAB's working directory if none has been set yet
if isempty(wd)
    wd = 'C:\Users\Prathyush\Desktop\Work\WarrenlabRepos\Matlab\puprisa_PUmpPRobeImageStackAnalysis\puprisa';
end

[fileName, pathName] = uigetfile([wd,'/*.txt'],...
    'Select a delay set' );


fitErrorFlag = 0;

if(strcmp(fileName(end-8:end-4),'Short'))
    if(strcmp(fileName(end-11:end-9),'tau'))
        FileName_S = fileName;
        FileName_AS = [fileName(1:end-12),'ampShort.txt'];
        FileName_L = [fileName(1:end-12),'tauLong.txt'];
        FileName_AL = [fileName(1:end-12),'ampLong.txt'];
    elseif(strcmp(fileName(end-11:end-9),'amp'))
        FileName_AS = fileName;
        FileName_S = [fileName(1:end-12),'tauShort.txt'];
        FileName_L = [fileName(1:end-12),'tauLong.txt'];
        FileName_AL = [fileName(1:end-12),'ampLong.txt'];
    end
elseif(strcmp(fileName(end-7:end-4),'Long'))
    if(strcmp(fileName(end-10:end-8),'tau'))
        FileName_L = fileName;
        FileName_AS = [fileName(1:end-11),'ampShort.txt'];
        FileName_S = [fileName(1:end-11),'tauShort.txt'];
        FileName_AL = [fileName(1:end-11),'ampLong.txt'];
    elseif(strcmp(fileName(end-10:end-8),'amp'))
        FileName_AL = fileName;
        FileName_S = [fileName(1:end-11),'tauShort.txt'];
        FileName_L = [fileName(1:end-11),'tauLong.txt'];
        FileName_AS = [fileName(1:end-11),'ampShort.txt'];
    end
elseif(strcmp(fileName(end-8:end-4),'Error'))
    fitErrorFlag = 1;
    if(strcmp(fileName(end-13:end-9),'Short'))
      FileName_AS = [fileName(1:end-17),'ampShort.txt'];
      FileName_S = [fileName(1:end-17),'tauShort.txt'];
      FileName_AL = [fileName(1:end-17),'ampLong.txt'];
      FileName_L = [fileName(1:end-17),'tauLong.txt'];
      
      FileName_ASerr = [fileName(1:end-17),'ampShortError.txt'];
      FileName_Serr = [fileName(1:end-17),'tauShortError.txt'];
      FileName_ALerr = [fileName(1:end-17),'ampLongError.txt'];
      FileName_Lerr = [fileName(1:end-17),'tauLongError.txt'];
    elseif(strcmp(fileName(end-12:end-9),'Long'))
      FileName_AS = [fileName(1:end-16),'ampShort.txt'];
      FileName_S = [fileName(1:end-16),'tauShort.txt'];
      FileName_AL = [fileName(1:end-16),'ampLong.txt'];
      FileName_L = [fileName(1:end-16),'tauLong.txt'];  
      
      FileName_ASerr = [fileName(1:end-16),'ampShortError.txt'];
      FileName_Serr = [fileName(1:end-16),'tauShortError.txt'];
      FileName_ALerr = [fileName(1:end-16),'ampLongError.txt'];
      FileName_Lerr = [fileName(1:end-16),'tauLongError.txt'];
    end
end

if ~isequal(fileName,0)
    tauShort = load([pathName,'/',FileName_S]);
    tauLong = load([pathName,'/',FileName_L]);
    ampShort = load([pathName,'/',FileName_AS]);
    ampLong = load([pathName,'/',FileName_AL]);
    if(fitErrorFlag)
        tauShortError = load([pathName,'/',FileName_Serr]);
        tauLongError = load([pathName,'/',FileName_Lerr]);
        ampShortError = load([pathName,'/',FileName_ASerr]);
        ampLongError = load([pathName,'/',FileName_ALerr]);
    end
end


m = size(tauShort);
j = 1;

for i = 1:m(1)
    if(tauShort(i) > 0.1  && (tauShort(i) < 40) && tauLong(i) > 0.1  && (tauLong(i) < 40))%valid range
        if(fitErrorFlag)
            if(isnan(tauShortError(i)) || isnan(ampShortError(i)) || isnan(tauLongError(i)) || isnan(ampLongError(i)))
                continue
            end
            if((tauShortError(i) > tauShort(i)*0.5) || (tauLongError(i) > tauLong(i)*0.5))
               continue
            end
            if( (abs(ampShortError(i)) > 0.5*abs(ampShort(i))) || (abs(ampLongError(i)) > 0.5*abs(ampLong(i)) ))
                continue
            end
            tauShort2Error(j) = tauShortError(i);
            ampShort2Error(j) = ampShortError(i);
            tauLong2Error(j) = tauLongError(i);
            ampLong2Error(j) = ampLongError(i);
        end
        tauShort2(j) = tauShort(i);
        ampShort2(j) = ampShort(i);
        tauLong2(j) = tauLong(i);
        ampLong2(j) = ampLong(i);
        j = j+1;
    end
end
j
MeanTimeConstShort = mean(tauShort2);
sigmaShort = std(tauShort2);
MeanampShort = mean(ampShort2);
sigmaampShort = std(ampShort2);

MeanTimeConstLong = mean(tauLong2);
sigmaLong = std(tauLong2);
MeanampLong = mean(ampLong2);
sigmaampLong = std(ampLong2);

  figure (1)
  hist(tauShort2,200);
  figure (2)
  hist(tauLong2,200);

excelWrite = ({'MeanTimeConstShort','sigmaShort','MeanampShort','sigmaampShort','MeanTimeConstLong'...
    ,'sigmaLong','MeanampLong','sigmaampLong'; MeanTimeConstShort,sigmaShort,MeanampShort,sigmaampShort,MeanTimeConstLong...
    ,sigmaLong,MeanampLong,sigmaampLong})';

if(fitErrorFlag)
    MeanTimeConstShortError = mean(tauShort2Error);
    sigmaShortError = std(tauShort2Error);
    MeanampShortError = mean(ampShort2Error);
    sigmaampShortError = std(ampShort2Error);
    
    MeanTimeConstLongError = mean(tauLong2Error);
    sigmaLongError = std(tauLong2Error);
    MeanampLongError = mean(ampLong2Error);
    sigmaampLongError = std(ampLong2Error);
    
    excelWriteERR = ({'MeanTimeConstShortError','sigmaShortError','MeanampShortError','sigmaampShortError'...
    ,'MeanTimeConstLongError','sigmaLongError','MeanampLongError','sigmaampLongError'; MeanTimeConstShortError,...
    sigmaShortError,MeanampShortError,sigmaampShortError,MeanTimeConstLongError,sigmaLongError,MeanampLongError,sigmaampLongError})';
    excelWrite = [excelWrite excelWriteERR]
    
end


xlswrite('data.xlsx',excelWrite);





