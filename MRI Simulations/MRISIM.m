%This file written to "simulate" a 3D Cartesian MRI sequence
cat = imread('CatGrayscale.gif');

%Data is written in an odd way:

% (k1,k2) are the k-space coordinates of a voxel
% (x,y,z) are the cartesian components of the magnetization vector of a
%voxel

%Gradient strengths-- assume 30 mT/m
%assume voxel spacing is 1mm center to center (512 ~ .5 meters)