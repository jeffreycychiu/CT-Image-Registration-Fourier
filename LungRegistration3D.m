%% 
%This file implements the Algorithm described in "Non-rigid registration of tomographic images with Fourier
%Transforms" by AR Osorio , RA Isoardi, and G Mato on the lung image datasets given in the EMPIRE10 Challenge
%NOTE: Run LoadLungImages.m first to get the data

close all;
clc;

fixed = mat2gray(rawFixed);
moving = mat2gray(rawMoving);
fixedMask = mat2gray(rawFixedMask);
movingMask = mat2gray(rawMovingMask);

%Use the provided lung masks to separate out the lungs only
fixedLung = fixed.*fixedMask;
movingLung = moving.*movingMask;

%% Show the slices on the centre of the volume sets
fixedLungCentre = centreSlice(fixedLung);
movingLungCentre = centreSlice(movingLung);
figure()
imshowpair(fixedLungCentre,movingLungCentre);
title('Unregistered Images, Centre Slice');

%% Zero pad the two 3d volumes so they are the same size
padXYZ = size(movingLung)-size(fixedLung);
fixedLungPad = fixedLung;
movingLungPad = movingLung;
if padXYZ(1) > 0
    fixedLungPad = padarray(fixedLungPad,[padXYZ(1)/2 0 0]);
elseif padXYZ(1) < 0
    movingLungPad = padarray(movingLungPad,[abs(padXYZ(1))/2 0 0]);
end

if padXYZ(2) > 0
    fixedLungPad = padarray(fixedLungPad,[0 padXYZ(2)/2 0]);
elseif padXYZ(2) < 0
    movingLungPad = padarray(movingLungPad,[0 abs(padXYZ(2))/2 0]);
end

if padXYZ(3) > 0
    fixedLungPad = padarray(fixedLungPad,[0 0 padXYZ(3)/2]);
elseif padXYZ(3) < 0
    movingLungPad = padarray(movingLungPad,[0 0 abs(padXYZ(3))/2]);
end

%% Affine Transformation Maximizing the Mutual Information
[optimizer,metric] = imregconfig('monomodal');
metric = registration.metric.MattesMutualInformation();
optimizer.MaximumIterations = 20;   %limit this for now to improve the speed
affineTransform = imregtform(movingLungPad, fixedLungPad, 'rigid', optimizer, metric,'DisplayOptimization',1); %Currently, this uses meansquares instead of correlation :/
movingLungTransformed = imwarp(movingLungPad, affineTransform);

%% zero pad again
padTransXYZ = size(movingLungTransformed)-size(fixedLung);
fixedLungTransformedPad = fixedLung;
movingLungTransformedPad = movingLungTransformed;
if padTransXYZ(1) > 0
    fixedLungTransformedPad = padarray(fixedLungTransformedPad,[floor(padTransXYZ(1)/2) 0 0], 0, 'post');
    fixedLungTransformedPad = padarray(fixedLungTransformedPad,[ceil(padTransXYZ(1)/2) 0 0], 0, 'pre');
elseif padTransXYZ(1) < 0
    movingLungTransformedPad = padarray(movingLungTransformedPad,[floor(abs(padTransXYZ(1))/2) 0 0], 0, 'post');
    movingLungTransformedPad = padarray(movingLungTransformedPad,[ceil(abs(padTransXYZ(1))/2) 0 0], 0, 'pre');
end

if padTransXYZ(2) > 0
    fixedLungTransformedPad = padarray(fixedLungTransformedPad,[0 floor(padTransXYZ(2)/2) 0], 0, 'post');
    fixedLungTransformedPad = padarray(fixedLungTransformedPad,[0 ceil(padTransXYZ(2)/2) 0], 0, 'pre');
elseif padTransXYZ(2) < 0
    movingLungTransformedPad = padarray(movingLungTransformedPad,[0 floor(abs(padTransXYZ(2))/2) 0], 0, 'post');
    movingLungTransformedPad = padarray(movingLungTransformedPad,[0 ceil(abs(padTransXYZ(2))/2) 0], 0, 'pre');
end

if padTransXYZ(3) > 0
    fixedLungTransformedPad = padarray(fixedLungTransformedPad,[0 0 floor(padTransXYZ(3)/2)], 0, 'post');
    fixedLungTransformedPad = padarray(fixedLungTransformedPad,[0 0 ceil(padTransXYZ(3)/2)], 0, 'pre');
elseif padTransXYZ(3) < 0
    movingLungTransformedPad = padarray(movingLungTransformedPad,[0 0 floor(abs(padTransXYZ(2))/2)], 0, 'post');
    movingLungTransformedPad = padarray(movingLungTransformedPad,[0 0 ceil(abs(padTransXYZ(2))/2)], 0, 'pre');
end

%% Display Images pre and post affine transform
centreFixedPad = round(size(fixedLungPad)/2);
centreMovingPad = round(size(movingLungPad)/2);
centreFixedTransformedPad = round(size(fixedLungTransformedPad)/2);
centreMovingTransformed=round(size(movingLungTransformed)/2);

fixedLungCentrePad = rot90(squeeze(fixedLungPad(:,centreFixedPad(2),:)));
movingLungCentrePad = rot90(squeeze(movingLungPad(:,centreMovingPad(2),:)));
fixedLungCentreTransformedPad = rot90(squeeze(fixedLungTransformedPad(:,centreFixedTransformedPad(2),:)));
movingLungTransformedCentre = rot90(squeeze(movingLungTransformed(:,centreMovingTransformed(2),:)));

figure()
imshowpair(fixedLungCentrePad,movingLungCentrePad);
title('Centre Slice zero Padded, no transform');
figure()
imshowpair(fixedLungCentreTransformedPad,movingLungTransformedCentre);
title('Centre Slice After Affine Transform (Mutual Information)');

%% Calculate average displacement error
avgDispErrorOriginal = (1/numel(fixedLungPad)) * sum(abs(movingLungPad(:) - fixedLungPad(:)))
avgDispErrorAffine = (1/numel(fixedLungTransformedPad)) * sum(abs(movingLungTransformedPad(:) - fixedLungTransformedPad(:)))

%% Divide the Fixed Image into 'k' sub-volumes, 'Fr) (k = 8)
k = 8;  %MUST CHANGE FUNCTION AS WELL IF THIS TO BE CHANGED
s = 2; %number of heiarchal subdivides
fixedLungSubVolume = subDivide8(fixedLungTransformedPad);
% for index = 1:k
%     fixedLungSubSubVolume = subDivide8(fixedLungSubVolume(:,:,:,index));
% end
   
%% Divide the Affine-Transformed Moving Images into 'k' sub-volumes
movingLungSubVolume = subDivide8(movingLungTransformedPad);


%% Loop - Rigid Registration on each sub volume, then non-rigid (fourier) registration

%Rigid Registration:
% [optimizer,metric] = imregconfig('monomodal');
% metric = registration.metric.MattesMutualInformation();
% optimizer.MaximumIterations = 20;   %limit this for now to limit the time taken
% rigidTransform = imregtform(movingLungSubVolume(:,:,:,1), fixedLungSubVolume(:,:,:,1), 'rigid', optimizer, metric,'DisplayOptimization',1); %Currently, this uses meansquares instead of correlation :/
% movingLungTransformedSub = imwarp(movingLungSubVolume(:,:,:,1), rigidTransform);
% 
% centreMovingTransformedSubVol_1 = centreSlice(movingLungTransformedSub);
% figure()
% imshowpair(centreFixedSubVol_1,centreMovingTransformedSubVol_1);
% title('subvolume 1, after global affine after subvolume rigid');


%% Independent Fourier Transformations applied to each sub volume Fr
n = 2;      %order of the system
sizeX = size(fixedLungSubVolume,1);
sizeY = size(fixedLungSubVolume,2);
sizeZ = size(fixedLungSubVolume,3);

%create coefficents up to order n. Number of coefficents=n^3;
[a_ijk,b_ijk,c_ijk,d_ijk,e_ijk,f_ijk] = deal(normrnd(0,1/n,n,n,n)); %gaussian distributed random for 1st guess at coeffs

%perform the calculation on each x-dimension of the subvolume
movingLungSubVolumeFourier = movingLungSubVolume;
xLocFourier = zeros(size(movingLungSubVolume));
yLocFourier = zeros(size(movingLungSubVolume));
zLocFourier = zeros(size(movingLungSubVolume));

for volumeNumber = 1:k     %for each subvolume, calculate the deformation fields
    for x = 1:size(movingLungSubVolume,1)  %for each volume element in the subvolume - this is a lot
        for y = 1: size(movingLungSubVolume,2)
            for z = 1:size(movingLungSubVolume,3)
                xLocFourier(x,y,z,k) = x;
                yLocFourier(x,y,z,k) = y;
                zLocFourier(x,y,z,k) = z;
                for i = 1:n
                    for j = 1:n
                        for k = 1:n
                            theta_ijk = sin(pi*i*x/sizeX) * sin(pi*j*y/sizeY) * sin(pi*k*z/sizeZ);
                            phi_ijk = cos(pi*i*x/sizeX) * cos(pi*j*y/sizeY) * cos(pi*k*z/sizeZ);
                            
                            xLocFourier(x,y,z,k) = xLocFourier(x,y,z,k) + a_ijk(i,j,k)*theta_ijk + b_ijk(i,j,k)*phi_ijk;
                            yLocFourier(x,y,z,k) = yLocFourier(x,y,z,k) + c_ijk(i,j,k)*theta_ijk + d_ijk(i,j,k)*phi_ijk;
                            zLocFourier(x,y,z,k) = zLocFourier(x,y,z,k) + e_ijk(i,j,k)*theta_ijk + f_ijk(i,j,k)*phi_ijk;
                            %fprintf('%d , %d, %d\n',x, y, z)
                            %x
                        end
                    end
                end
            end
        end
    end
end
%% Apply deformation fields
volumeDisplacement = cell(1,8);
volumeFourierRegistered = cell(1,8);
for index = 1:8
    volumeDisplacement{index} = cat(4,xLocFourier(:,:,:,index),yLocFourier(:,:,:,index),zLocFourier(:,:,:,index));
    volumeFourierRegistered{index} = imwarp(movingLungSubVolume(:,:,:,index),volumeDisplacement{index});
end

%% Show sub volume 1
centreFixedSubVol_1 = centreSlice(fixedLungSubVolume(:,:,:,1));
centreMovingSubVol_1 = centreSlice(movingLungSubVolume(:,:,:,1));
centreFourierSubVol_1 = centreSlice(volumeFourierRegistered{1}(:,:,:));

figure()
imshowpair(centreFixedSubVol_1,centreMovingSubVol_1);
title('subvolume 1, after global affine but before subvolume rigid');
figure()
imshowpair(centreFixedSubVol_1, centreFourierSubVol_1);
title('fixed to fourier');
figure()
imshowpair(centreMovingSubVol_1, centreFourierSubVol_1,'diff');
title('moving before and after fourier');


fixed1 = fixedLungSubVolume(:,:,:,1);
movingBefore1 = movingLungSubVolume(:,:,:,1);

avgDispError1Original = (1/numel(movingBefore1)) * sum(abs(movingBefore1(:) - fixed1(:)))
avgDispError1New = (1/numel(newImage)) * sum(abs(newImage(:) - fixed1(:)))

%Show the displacement field for subvolume 1
dispFieldx = squeeze(volumeDisplacement{1}(:,round(161/2),:,1));
dispFieldz = squeeze(volumeDisplacement{1}(:,round(161/2),:,3));
figure()
quiver(dispFieldx,dispFieldz);
figure()
histogram(volumeDisplacement{1}(:,round(161/2)));

outputImage = subCombine8(volumeFourierRegistered);
outputCentre = centreSlice(outputImage);
imshowpair(fixedLungCentreTransformedPad,outputCentre);

avgDispErrorOriginal = (1/numel(fixedLungPad)) * sum(abs(movingLungPad(:) - fixedLungPad(:)))
avgDispErrorAffine = (1/numel(fixedLungTransformedPad)) * sum(abs(movingLungTransformedPad(:) - fixedLungTransformedPad(:)))
avgDispErrorFourier = (1/numel(fixedLungTransformedPad)) * sum(abs(outputImage(:) - fixedLungTransformedPad(:)))
