% PUPRISA_IMAGESTACK
% class definition for an image stack
% 
% Jesse Wilson (2011) syrex314@gmail.com
% Duke University
%
% This class contains an image stack, which will be capable of viewing and
% manipulating by a number of different GUIs. This class should hold the
% framework for accessing the underlying image data.
%
% Later on we may want to define views and workspaces containing views and
% overlays of images.
%
% Created: 5/9/2011
% Last modified: 5/9/2011
%
% CHANGELOG:
%
classdef puprisa_ImageStack < handle
    
    properties (SetAccess = protected)
        imageSlices;    % 3D array of size n x m x p, 
                        % containg p slices of n x m pixel data arrays
    end
    
    % need a constructor capable of reading from a file
end