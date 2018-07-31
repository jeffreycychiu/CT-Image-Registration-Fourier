function outputImage = subCombine8(subVolume)

sizeVol = size(subVolume{1});
outputImage = 2*size(subVolume,2);

outputImage(1:sizeVol(1), 1:sizeVol(2), 1:sizeVol(3)) = subVolume{1};
outputImage(sizeVol(1)+1:2*sizeVol(1), 1:sizeVol(2), 1:sizeVol(3)) = subVolume{2};
outputImage(1:sizeVol(1), sizeVol(2)+1:2*sizeVol(2), 1:sizeVol(3)) = subVolume{3};
outputImage(sizeVol(1)+1:2*sizeVol(1), sizeVol(2)+1:2*sizeVol(2), 1:sizeVol(3)) = subVolume{4};
outputImage(1:sizeVol(1), 1:sizeVol(2), sizeVol(3)+1:2*sizeVol(3)) = subVolume{5};
outputImage(sizeVol(1)+1:2*sizeVol(1), 1:sizeVol(2), sizeVol(3)+1:2*sizeVol(3)) = subVolume{6};
outputImage(1:sizeVol(1), sizeVol(2)+1:2*sizeVol(2), sizeVol(3)+1:2*sizeVol(3)) = subVolume{7};
outputImage(sizeVol(1)+1:2*sizeVol(1), sizeVol(2)+1:2*sizeVol(2), sizeVol(3)+1:2*sizeVol(3)) = subVolume{8};
