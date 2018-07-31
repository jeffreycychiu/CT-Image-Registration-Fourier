function centreSlice = centreSlice(A) %outputs the centre slice of image A

centreAIndex= round(size(A)/2);
centreSlice = rot90(squeeze(A(:,centreAIndex(2),:)));