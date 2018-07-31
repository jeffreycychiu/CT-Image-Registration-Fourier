function sub = subDivide4(A)

sizeVol = floor(size(A)/2);
sub = cell(1,4);

sub{1}(:,:) = A(1:sizeVol(1), 1:sizeVol(2));
sub{2}(:,:) = A(sizeVol(1)+1:2*sizeVol(1), 1:sizeVol(2));
sub{3}(:,:) = A(1:sizeVol(1), sizeVol(2)+1:2*sizeVol(2));
sub{4}(:,:) = A(sizeVol(1)+1:2*sizeVol(1), sizeVol(2)+1:2*sizeVol(2));