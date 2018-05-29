%Script for calculating wavenumber differences given wavelengths and
%bandwidths in nanometers.

%Make a matrix to accomodate up to 50 combinations
A = zeros(50,18);
%Determine probe wavelength and bandwidth
promptpr = 'What is the probe wavelength? ';
prwave = input(promptpr);
A(1,6) = prwave;
promptprb = 'What is the probe bandwidth (FWHM)? ';
prband = input(promptprb);
A(1,7) = prband;
%If more than one pump wavelength, do a loop
promptsw = 'How many pump wavelengths?';
Number = input(promptsw);
if Number > 100;
    promptchk = 'Are you sure? How many pump wavelengths?';
    Number = input(promptchk);
end
if Number == 1;
    promptpu = 'What is the pump wavelength? ';
    A(1,1) = input(promptpu);
    promptpup = 'What is the pump bandwidth (FWHM)? ';
    A(1,2) = input(promptpup);
    %A(1,6) = input(promptpr);
    %A(1,7) = input(promptprb);
     %Throw error if pump wavelength not higher energy than probe wavelength
    if (A(1,1)>A(1,6))
        disp('Invalid Raman combination.'), break;
    end
    %Convert center wavelengths in nanometers to wavenumbers
    A(1,3) = 1/(A(1,1)*(10^-7));
    A(1,8) = 1/(A(1,6)*(10^-7));
    %Calculate wavenumber difference
    A(1,11) = abs(A(1,3)-A(1,8));
    %Calculate bandwidths
    A(1,15) = A(1,1)-(A(1,2)/2);
    A(1,16) = A(1,1)+(A(1,2)/2);
    A(1,17) = A(1,6)-(A(1,7)/2);
    A(1,18) = A(1,6)+(A(1,7)/2);
    %Convert bandwidths in nanometers to wavenumbers
    A(1,4) = 1/(A(1,15)*(10^-7));
    A(1,5) = 1/(A(1,16)*(10^-7));
    A(1,9) = 1/(A(1,17)*(10^-7));
    A(1,10) = 1/(A(1,18)*(10^-7));
    %Calculate bandwidth of wavenumber difference
    A(1,12) = abs(A(1,5)-A(1,9));
    A(1,13) = abs(A(1,4)-A(1,10));
    A(1,14) = abs(A(1,13)-A(1,12));
    %Print results
    s1 = 'The wavenumber difference is';
    sprintf('%s %f', s1, A(1,11))
    s2 = 'The FWHM bandwidth is';
    sprintf('%s %f', s2, A(1,14))
elseif Number > 1 ;
    %Loop through several pump wavelengths
    n = 1;
    while n<(Number+1)
        %Get conditions
        promptpu = 'What is the pump wavelength? ';
        A(n,1) = input(promptpu);
        promptpup = 'What is the pump bandwidth (FWHM)? ';
        A(n,2) = input(promptpup);
        %Fill in probe wavelengths
        A(n,6) = prwave;
        A(n,7) = prband;
        %Throw error if pump wavelength not higher energy than probe wavelength
        if (A(n,1)>A(n,6))
            disp('Invalid Raman combination.'), break;
        end
        %Convert center wavelengths in nanometers to wavenumbers
        A(n,3) = 1/(A(n,1)*(10^-7));
        A(n,8) = 1/(A(n,6)*(10^-7));
        %Calculate wavenumber difference
        A(n,11) = abs(A(n,3)-A(n,8));
        %Calculate bandwidths
        A(n,15) = A(n,1)-(A(n,2)/2);
        A(n,16) = A(n,1)+(A(n,2)/2);
        A(n,17) = A(n,6)-(A(n,7)/2);
        A(n,18) = A(n,6)+(A(n,7)/2);
        %Convert bandwidths in nanometers to wavenumbers
        A(n,4) = 1/(A(n,15)*(10^-7));
        A(n,5) = 1/(A(n,16)*(10^-7));
        A(n,9) = 1/(A(n,17)*(10^-7));
        A(n,10) = 1/(A(n,18)*(10^-7));
        %Calculate bandwidth of wavenumber difference
        A(n,12) = abs(A(n,5)-A(n,9));
        A(n,13) = abs(A(n,4)-A(n,10));
        A(n,14) = abs(A(n,13)-A(n,12));
        n = n+1;
        %Print results
        R = zeros(50,3);
        R(:,1) = A(:,1);
        R(:,2) = A(:,11);
        R(:,3) = A(:,14);
    end
    disp('Results in Table R')
end