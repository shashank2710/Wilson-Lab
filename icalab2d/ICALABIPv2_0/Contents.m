% ICALAB for Image Processing
% (abbreviation of Independent Component Analysis Laboratory)
%
% developed and tested under Matlab versions 6.0 and higher
% Version 2.0, March 20, 2004
% Copyright (c) Andrzej Cichocki, Shun-ichi Amari, Krzysztof Siwek, 
%               Toshihisa Tanaka, Tomasz Rutkowski, Sergio Cruces, Seungjin Choi, 
%               Pando Georgiev and Yasushi Terazono and others. 
%
% The package contains several algorithms for ICA, BSS, BSE, 
% PCA and whitening of images. In this demo version some algorithms are 
% disabled.
%
% Full version of the package is available on the web page of the authors
% http://www.bsp.brain.riken.jp
%
% Type icalab in command window to launch main program 
% with the graphical user interface
%
% Please refer to your Matlab documentation on how to add 
% ICALAB for Image Processing to your Matlab search path. 
%
% ICALAB for Image Processing program:
%   icalab  - Main programm with graphical user interface for ICALAB
%
% Most m-files used by ICALAB for Image Processing are p-coded
% and have .p extensions. 
%
% User algorithms
% 
% The user can integrate with the package his or her own algorithms
% or any suitable algorithm  available in matlab, by simply inserting its
% code in:
% 
%    user_alg1.m - user_alg10.m  files.
%
% The algorithm inserted be the user should return only one variable:
% demixing matrix W
% If your algorithm estimate the mixing matrix H you
% can use W=pinv(H).
%
