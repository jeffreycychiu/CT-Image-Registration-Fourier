function sub = subDivide8(A)
%Divides the 3-d image A into k subvolumes

sizeVol = round(size(A)/2);
sub = zeros([sizeVol 8]);

sub(:,:,:,1) = A(1:sizeVol(1), 1:sizeVol(2), 1:sizeVol(3));
sub(:,:,:,2) = A(sizeVol(1)+1:2*sizeVol(1), 1:sizeVol(2), 1:sizeVol(3));
sub(:,:,:,3) = A(1:sizeVol(1), sizeVol(2)+1:2*sizeVol(2), 1:sizeVol(3));
sub(:,:,:,4) = A(sizeVol(1)+1:2*sizeVol(1), sizeVol(2)+1:2*sizeVol(2), 1:sizeVol(3));
sub(:,:,:,5) = A(1:sizeVol(1), 1:sizeVol(2), sizeVol(3)+1:2*sizeVol(3));
sub(:,:,:,6) = A(sizeVol(1)+1:2*sizeVol(1), 1:sizeVol(2), sizeVol(3)+1:2*sizeVol(3));
sub(:,:,:,7) = A(1:sizeVol(1), sizeVol(2)+1:2*sizeVol(2), sizeVol(3)+1:2*sizeVol(3));
sub(:,:,:,8) = A(sizeVol(1)+1:2*sizeVol(1), sizeVol(2)+1:2*sizeVol(2), sizeVol(3)+1:2*sizeVol(3));
