%This file loads the lung images into memory. Lung datasets are given from the EMPIRE10 challenge
%The ReadData3D Function was written by Dirk-Jan Kroon and is avaliable at
%http://www.mathworks.com/matlabcentral/fileexchange/29344-read-medical-data-3d

clear all
close all
clc

addpath('C:\Users\Jeff\Documents\Grad School\Courses\EECE 544 Medical Imaging\Final Project\ReadData3D');
%addpath('C:\Users\Jeff\Documents\Grad School\Courses\EECE 544 Medical Imaging\Final Project\imshow3D');
%addpath('C:\Users\Jeff\Documents\Grad School\Courses\EECE 544 Medical Imaging\Final Project\MIRT\mirt3d');

addpath('C:\Users\Jeff\Documents\Grad School\Courses\EECE 544 Medical Imaging\Final Project\01to06\Online\scans');
[rawFixed,infoFixed]=ReadData3D('02_Fixed.mhd',true); 
[rawMoving,infoMoving] = ReadData3D('02_Moving.mhd',true);
rmpath('C:\Users\Jeff\Documents\Grad School\Courses\EECE 544 Mediscal Imaging\Final Project\01to06\Online\scans');

addpath('C:\Users\Jeff\Documents\Grad School\Courses\EECE 544 Medical Imaging\Final Project\01to06\Online\lungMasks');
[rawFixedMask,infoFixedMask]=ReadData3D('02_Fixed.mhd',true); 
[rawMovingMask,infoMovingMask] = ReadData3D('02_Moving.mhd',true);
rmpath('C:\Users\Jeff\Documents\Grad School\Courses\EECE 544 Medical Imaging\Final Project\01to06\Online\lungMasks');
