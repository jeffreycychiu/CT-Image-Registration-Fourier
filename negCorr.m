function negCorrelation = negCorr(X,movingImage,fixedImage)
    
X(1,3)=0;
X(2,3)=0;
X(3,3)=1;
affineTransform = affine2d(X);

movingTransformed = imwarp(movingImage,affineTransform);
%zero pad
padXY = size(movingTransformed)-size(fixedImage);

if padXY(1) > 0
    fixedImage = padarray(fixedImage,[floor(padXY(1)/2) 0],0,'post');
    fixedImage = padarray(fixedImage,[ceil(padXY(1)/2) 0],0,'pre');
elseif padXY(1) < 0
    movingTransformed = padarray(movingTransformed,[floor(abs(padXY(1))/2) 0],0,'post');
    movingTransformed = padarray(movingTransformed,[ceil(abs(padXY(1))/2) 0],0,'pre');
end
if padXY(2) > 0
    fixedImage = padarray(fixedImage,[0 floor(padXY(2)/2)],0,'post');
    fixedImage = padarray(fixedImage,[0 ceil(padXY(2)/2)],0,'pre');
elseif padXY(2) < 0
    movingTransformed = padarray(movingTransformed,[0 floor(abs(padXY(2))/2)],0,'post');
    movingTransformed = padarray(movingTransformed,[0 ceil(abs(padXY(2))/2)],0,'pre');
end

negCorrelation = -1*corr2(fixedImage,movingTransformed);    
