function threshIm = areaThreshold(Im,minArea,maxArea)
%Get areas of each component for thresholding

threshIm = zeros(size(Im),class(Im));

for frameNum = 1:size(Im,3)
    CCPreThresh = bwconncomp(Im(:,:,frameNum));
    statsPreThresh = regionprops(CCPreThresh,'Area');

    %Get components above the area threshold
    keptComponents = find([statsPreThresh.Area]>minArea & [statsPreThresh.Area]<maxArea);

    %Make a new image with only components above threshold
    threshIm(:,:,frameNum) = ismember(labelmatrix(CCPreThresh),keptComponents);
end

end