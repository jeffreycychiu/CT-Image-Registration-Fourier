function movingFourierImage = fourierTransImage(fourierCoeffs, movingImage, fixedImage, n)

%n is the order. FourierCoeff
sizeX = size(movingImage,1);
sizeY = size(movingImage,2);

%initial rigid registration

%initial guess at fouriers, use rand
%a_ij = normrnd(0,1/s,n,n); b_ij = normrnd(0,1/s,n,n); c_ij = normrnd(0,1/s,n,n); d_ij = normrnd(0,1/s,n,n); 
a_ij = fourierCoeffs(:,:,1);
b_ij = fourierCoeffs(:,:,2);
c_ij = fourierCoeffs(:,:,3);
d_ij = fourierCoeffs(:,:,4);

xLocFourier = zeros(size(movingImage));
yLocFourier = zeros(size(movingImage));

for x = 1:size(movingImage,1)
    for y = 1:size(movingImage,2)
        xLocFourier(x,y) = 0;
        yLocFourier(x,y) = 0;
        for i = 1:n
            for j = 1:n
                theta_ij = sin(pi*i*x/sizeX)*sin(pi*j*y/sizeY);
                phi_ij = cos(pi*i*x/sizeX)*cos(pi*j*y/sizeY);
                
                xLocFourier(x,y) = xLocFourier(x,y) + a_ij(i,j)*theta_ij + b_ij(i,j)*phi_ij;
                yLocFourier(x,y) = yLocFourier(x,y) + c_ij(i,j)*theta_ij + d_ij(i,j)*phi_ij;                
            end
        end
    end
end

%now we have a matrix of xLocFourier and yLocFourier
%Tx(x,y) = x+xLocFourier
% movingFourierImage = zeros(size(movingImage));
% for xIndex = 1:size(movingImage,1)
%     for yIndex = 1:size(movingImage,2)
%         movingFourierImage(round(x+xLocFourier(x,y)),round(y+yLocFourier(x,y))) = movingImage(x,y);    
%     end
% end

displacement = cat(3,xLocFourier,yLocFourier);
movingFourierImage = imwarp(movingImage,displacement);

%% zero pad
padXY = size(movingFourierImage)-size(fixedImage);
if padXY(1) > 0
    fixedImage = padarray(fixedImage,[floor(padXY(1)/2) 0],0,'post');
    fixedImage = padarray(fixedImage,[ceil(padXY(1)/2) 0],0,'pre');
elseif padXY(1) < 0
    movingFourierImage = padarray(movingFourierImage,[floor(abs(padXY(1))/2) 0],0,'post');
    movingFourierImage = padarray(movingFourierImage,[ceil(abs(padXY(1))/2) 0],0,'pre');
end
if padXY(2) > 0
    fixedImage = padarray(fixedImage,[0 floor(padXY(2)/2)],0,'post');
    fixedImage = padarray(fixedImage,[0 ceil(padXY(2)/2)],0,'pre');
elseif padXY(2) < 0
    movingFourierImage = padarray(movingFourierImage,[0 floor(abs(padXY(2))/2)],0,'post');
    movingFourierImage = padarray(movingFourierImage,[0 ceil(abs(padXY(2))/2)],0,'pre');
end

%% Wcalculate the negative correlation to be minimized
fourierCorr = -1*corr2(fixedImage,movingFourierImage)
