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
centreFixed = round(size(fixedLung)/2);
centreMoving = round(size(movingLung)/2);

fixedLungCentre = rot90(squeeze(fixedLung(:,centreFixed(2),:)));
movingLungCentre = rot90(squeeze(movingLung(:,centreMoving(2),:)));

imshowpair(fixedLungCentre,movingLungCentre);
title('Unregistered Images, Centre Slice');

%% Zero pad the 2D centre so they are the same size
padXY = size(movingLungCentre)-size(fixedLungCentre);
fixedLungCentrePad = fixedLungCentre;
movingLungCentrePad = movingLungCentre;
if padXY(1) > 0
    fixedLungCentrePad = padarray(fixedLungCentrePad,[padXY(1)/2 0]);
elseif padXY(1) < 0
    movingLungCentrePad = padarray(movingLungCentrePad,[abs(padXY(1))/2 0]);
end

if padXY(2) > 0
    fixedLungCentrePad = padarray(fixedLungCentrePad,[0 padXY(2)/2]);
elseif padXY(2) < 0
    movingLungCentrePad = padarray(movingLungCentrePad,[0 abs(padXY(2))/2]);
end

figure()
imshowpair(fixedLungCentrePad,movingLungCentrePad);
title('Unregistered Images, Centre Slice, zero padded');
avgDispErrorOriginal = (1/numel(fixedLungCentrePad)) * sum(abs(movingLungCentrePad(:) - fixedLungCentrePad(:)))

%% Estimates geometric transformation that aligns two 2-D images using correlation coefficent
corrCoeff1 = corr2(fixedLungCentrePad,movingLungCentrePad);

%function movingTransform = affineTransformMatrix * movingTransform;
a = 1; b = 0; c = 0; d = 1; e = 0; f = 0;

func = @negCorr;
opt=optimset('LargeScale','off','Display','iter');
start = [a b 0; c d 0; e f 1];
XF = fminsearch(@(X) negCorr(X,movingLungCentrePad,fixedLungCentrePad), start, opt);
%% Optimize the affine transform
XF(1,3)=0;
XF(2,3)=0;
XF(3,3)=1;
affineTransform = affine2d(XF);
movingTransformed = imwarp(movingLungCentrePad,affineTransform);

%% zero pad again
padXY = size(movingTransformed)-size(fixedLungCentrePad);
if padXY(1) > 0
    fixedLungCentrePad = padarray(fixedLungCentrePad,[floor(padXY(1)/2) 0],0,'post');
    fixedLungCentrePad = padarray(fixedLungCentrePad,[ceil(padXY(1)/2) 0],0,'pre');
elseif padXY(1) < 0
    movingTransformed = padarray(movingTransformed,[floor(abs(padXY(1))/2) 0],0,'post');
    movingTransformed = padarray(movingTransformed,[ceil(abs(padXY(1))/2) 0],0,'pre');
end
if padXY(2) > 0
    fixedLungCentrePad = padarray(fixedLungCentrePad,[0 floor(padXY(2)/2)],0,'post');
    fixedLungCentrePad = padarray(fixedLungCentrePad,[0 ceil(padXY(2)/2)],0,'pre');
elseif padXY(2) < 0
    movingTransformed = padarray(movingTransformed,[0 floor(abs(padXY(2))/2)],0,'post');
    movingTransformed = padarray(movingTransformed,[0 ceil(abs(padXY(2))/2)],0,'pre');
end

%% Show Figure after affine/correlation
figure()
imshowpair(fixedLungCentrePad,movingTransformed);
title('After affine transform maximizing the correlation');
corrCoeff1
corrCoeff2 = corr2(fixedLungCentrePad,movingTransformed)
avgDispErrorAffine = (1/numel(fixedLungCentrePad) * sum(abs(movingTransformed(:) - fixedLungCentrePad(:))))

%% Divide image into s = 4 equal parts
s = 4;
subFixed = subDivide4(fixedLungCentrePad);
subMoving = subDivide4(movingTransformed);

%% Calculate order n = 3 fourier coeffs
n = 5;

startFourierCoeff = cell(1,s);
fourierCoeff = cell(1,s);
movingFourierTrans = cell(1,s);
for i = 1:s
    startFourierCoeff{i} = zeros(n,n,s);
    startFourierCoeff{i}(:,:,1) = normrnd(0,1/s,n,n);
    startFourierCoeff{i}(:,:,2) = normrnd(0,1/s,n,n);
    startFourierCoeff{i}(:,:,3) = normrnd(0,1/s,n,n);
    startFourierCoeff{i}(:,:,4) = normrnd(0,1/s,n,n);
    opt=optimset('LargeScale','off','Display','off','TolFun',5e-4,'MaxIter',750);
    fourierCoeff{i} = fminsearch(@(X) fourierMethod(X,subMoving{1},subFixed{1}, n), startFourierCoeff{i}, opt);

    movingFourierTrans{i} = fourierTransImage(fourierCoeff{i}, subMoving{i}, subFixed{i}, n);
end

%Rough combination
sizeVol = size(subFixed{1});
movingFourierCombined (1:sizeVol(1), 1:sizeVol(2)) = movingFourierTrans{1};
movingFourierCombined (sizeVol(1)+1:2*sizeVol(1), 1:sizeVol(2)) = movingFourierTrans{2};
movingFourierCombined (1:sizeVol(1), sizeVol(2)+1:2*sizeVol(2)) = movingFourierTrans{3};
movingFourierCombined (sizeVol(1)+1:2*sizeVol(1), sizeVol(2)+1:2*sizeVol(2)) = movingFourierTrans{4};

figure()
imshowpair(fixedLungCentrePad,movingFourierCombined)
title('Recombined Fourier')
corrCoeff3 = corr2(fixedLungCentrePad,movingFourierCombined)

figure()
bar(diag([corrCoeff1,corrCoeff2,corrCoeff3]));
legend('Unregistered','Affine','Affine+Fourier')
ylabel('Correlation Coefficient with Fixed')
title('Correlation Coefficent Before and After Transforms')

%% Calculate Average Displacement Error
avgDispErrorFourier = (1/numel(fixedLungCentrePad) * sum(abs(movingFourierCombined(:) - fixedLungCentrePad(:))))

figure()
bar(diag([avgDispErrorOriginal,avgDispErrorAffine,avgDispErrorFourier]));
legend('Unregistered','Affine','Affine+Fourier')
ylabel('Averagae Displacement Error')
title('Average Displacement Error Before and After Transforms')

%% Evaluate Lung Boundaries
centreMask = round(size(fixedMask)/2);
fixedLungMask = rot90(squeeze(fixedMask(:,centreMask(2),:)));
lungEdge = edge(fixedLungMask);
figure()
imshow(lungEdge)

add2mm = zeros(size(lungEdge));
for x = 1:size(lungEdge,1)
    for y = 1:size(lungEdge,2)
        if lungEdge(x,y) == 1;
            %draw 3 pixels out to each side. each pixel is 0.7mm *3 = 2.1mm
            for i = -3:3
               for j = -3:3
                   add2mm(x+i,y+j)=1;
               end 
            end
            
        end
    end
end
lungEdge = lungEdge + add2mm;
figure()
imshow(lungEdge)


%% zero pad again
padXY = size(movingFourierCombined)-size(lungEdge);
if padXY(1) > 0
    lungEdge = padarray(lungEdge,[floor(padXY(1)/2) 0],0,'post');
    lungEdge = padarray(lungEdge,[ceil(padXY(1)/2) 0],0,'pre');
elseif padXY(1) < 0
    movingFourierCombined = padarray(movingFourierCombined,[floor(abs(padXY(1))/2) 0],0,'post');
    movingFourierCombined = padarray(movingFourierCombined,[ceil(abs(padXY(1))/2) 0],0,'pre');
end
if padXY(2) > 0
    lungEdge = padarray(lungEdge,[0 floor(padXY(2)/2)],0,'post');
    lungEdge = padarray(lungEdge,[0 ceil(padXY(2)/2)],0,'pre');
elseif padXY(2) < 0
    movingFourierCombined = padarray(movingFourierCombined,[0 floor(abs(padXY(2))/2)],0,'post');
    movingFourierCombined = padarray(movingFourierCombined,[0 ceil(abs(padXY(2))/2)],0,'pre');
end

figure()
imshowpair(movingFourierCombined,lungEdge)

