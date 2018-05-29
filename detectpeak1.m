% function [Rpeaks, RR, heart_rate, avg_HR]=detectpeak(vq7)

    file='Volt analysis10.xlsx';
    %file='Neckvolt3.xlsx';
    [status,sheets,xlFormat] = xlsfinfo(file);
    sheets
    %volt analysis8 -4102
    %volt analysis9 -2058
    %Volt analysis10 - 8198
    ecg1=xlsread(file,'F7:F8198');
    %y=xlsread(file,'F7:F2054');
    %ecg=detrend(y);
    ecg1(1:230)=[];
    y=detrend(ecg1);
    %ddata=detrend(y);
    %ecg1=detrend(y);
    %figure, plot(ecg1);
    ecg=denoise(y,0.8,50);
    %ecg2=denoise(y,0.8,34);
    %xq8=1:20:length(ddata);
    %vq8=interp1(1:length(ddata), ddata, xq8,'nearest');
    %time=(1:length(vq8))/256;
    %figure, plot(time, vq8);
    figure, plot(ecg);
    hold on;
    %figure, plot(ecg2);
    xlabel('Sample value');
    ylabel('Voltage');
    N = length(ecg);
    fs=100;
    m = 1;
for j = 1:fs:N-fs
    PeakMax(m) = max(ecg(j:j+fs));
    m = m+1;
end
    [hn,xout] = hist(PeakMax,16);
    [C,I] = max(hn);
    PeakM = xout(I);
    R_series =0;
    h= 1;
    k = 1;
    t=1;
    Rpeaks = 0;
    ecg1=0;
% fprintf('----------------------- Start of QRS detector  ---------------------\n');

for i = 1:1:N
       
    if (ecg(i)>=0.5*PeakM)
        if (i-1)<0.4*fs || (i+1)> N-0.4*fs %consider current signal position belong to QRS
            continue;
        end
        
        if ecg(i) == max(ecg(i-0.4*fs:i+0.4*fs))            
            ecg1(i)=ecg(i);
            Rpeaks(k) = i;
            k = k+1;    
        end    
    end
    
end

index=find(ecg1~=0);
plot(index,ecg1(index),'ro');
hold off;

i =1;
while(i<length(Rpeaks))
    if (Rpeaks(i+1)<Rpeaks(i)+10)
        Rpeaks(i+1)=[];
        i=i-1;
    end
    i=i+1;
end

for i=1:length(Rpeaks)
R_series(i)= ecg(Rpeaks(i));  %get the R amplitude  
end
        
RR = diff(Rpeaks)/100; %gives the difference in Rpeaks in sec (has precision of milliseconds)
%i.e. if it shows 0.800 that means 800 milliseconds or 0.800 seconds

for i=1:length(RR)
heart_rate(i)=60/RR(i);
end
% heart_rate
avg_HR=mean(heart_rate);
avg_HR
