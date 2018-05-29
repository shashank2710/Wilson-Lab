function filtersignal = denoise(ori_signal,lowband,highband)
% ori_signal=xlsread('CL-051415-748201_I.csv');
% lowband=0.5;
% highband=300;
% Author: Hui Yang
% Affiliation: 
       %Department of Industrial Engineering and Management 
       %Oklahoma State University,Stillwater, 74075
       %Email: yanghui@gmail.com
       %SAID: Signal Analysis and Integrated Diagnosis Research Center
% ori_signal: the input time series, should be a vector
% lowband: cutoff frequency low boundary, suggested value: 0.5
% highband: cutoff frequency high boundary, suggested value: 200

fs=100; % Physionet ecg signal were origanlly sampled at 1000 Hz. 
passband(1) = lowband;
passband(2) = highband;

N = length(ori_signal);
y = fft(ori_signal);

lowicut = round(passband(1)*N/fs);
lowmirror = N-lowicut+2;

highicut =  round(passband(2)*N/fs);
highmirror = N-highicut+2;

y([1:(lowicut-1) (lowmirror+1):end])=0;
y((highicut+1):(highmirror-1))=0;

filtersignal = ifft(y);
% filtersignal=2*filtersignal;
% figure, plot(filtersignal)

%clear all
% csvwrite('leadII.csv',filtersignal)