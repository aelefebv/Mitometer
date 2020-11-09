function [sigmaOptimal,threshOptimal,costs] = optimizeSigmaThresh(ImBgRemoved,minArea,maxArea)
%

numFrames = min(size(ImBgRemoved,3),10);

%We use Otsu’s thresholding and verify its effectiveness metric to
%determine whether it is satisfactory (>80% effective) to use rather than
%iterating through intensity values. If it is not satisfactory (<80%
%effective) we calculate the entire stack’s 80th quantile of intensity and
%use that as the upper limit for thresholding.
[t,em] = graythresh(uint8(ImBgRemoved));
if em > 0.8 && round(t*max(ImBgRemoved(:))) ~= 0
    intensityQuantile = round(t*max(ImBgRemoved(:)));
    threshMatrix = intensityQuantile;
else
    intensityQuantile = quantile(ImBgRemoved(ImBgRemoved~=0),0.8);
    threshMatrix = (2:1:intensityQuantile);
end

%We also limit the gaussian filter’s standard deviation between 0.33 and
%0.5. We choose these values as we want to keep the gaussian filter’s size
%to a 3x3 kernel. We choose a minimum standard deviation of greater than
%0.32 as it has a center value of 0.97 and edge values of 0.0074,
%essentially making it and anything below it useless as a filter. We choose
%a maximum standard deviation of 0.5 as, in general, the rule of thumb is
%to limit the filter size to 3 times the standard deviation on either side
%of the center (making it 6*0.5 for a size of 3). This is because the tails
%of the Gaussian have values that are effectively 0 after three standard
%deviations, and anything greater than 0.5 with a filter size of 3 produces
%a non-Gaussian filter.
sigmaMatrix = (0.33:0.01:0.5);


numComponents = zeros([length(threshMatrix),length(sigmaMatrix),numFrames]);
medianArea = zeros([length(threshMatrix),length(sigmaMatrix),numFrames]);

loading = waitbar(0,'Please wait...','CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
setappdata(loading,'canceling',0);
pause(.0002)


for sigma = 1:length(sigmaMatrix)
    ImGaussFiltered = gaussFilter(ImBgRemoved,sigmaMatrix(sigma));
    for thresh = 1:length(threshMatrix)
        
        %
        if getappdata(loading,'canceling')
            delete(loading)
            break
        end
        waitbar((thresh+length(threshMatrix)*(sigma-1))/(length(sigmaMatrix)*length(threshMatrix)),loading,sprintf('Testing parameter %d of %d.',(thresh+length(threshMatrix)*(sigma-1)),(length(sigmaMatrix)*length(threshMatrix))));
        pause(.0002)
        %
        
        %         [~,ImMask] =
        %         segmentObjects(Im,ImBgRemoved,ImGaussFiltered,threshMatrix(thresh));
        ImMask = thresholdImage(ImGaussFiltered,threshMatrix(thresh));
        
        %We then iterate through the first 10 timeframes (or all frames if
        %the total number of timeframes is less than 10) for each standard
        %deviation and intensity value.
        for frameNum = 1:numFrames
            CCPreThresh = bwconncomp(ImMask(:,:,frameNum));
            statsPreThresh = regionprops(CCPreThresh,'Area');
            %We calculate the number of mitochondria and the median area of
            %all mitochondria that are within the area threshold in each
            %temporal frame.
            if nargin > 1
                keptComponents = find([statsPreThresh.Area]>minArea & [statsPreThresh.Area]<maxArea);
                numComponents(thresh,sigma,frameNum) = length(keptComponents);
                medianArea(thresh,sigma,frameNum) = nanmedian([statsPreThresh(keptComponents).Area]);
            else
                keptComponents = find([statsPreThresh.Area]>1);
                numComponents(thresh,sigma,frameNum) = length(keptComponents);
                medianArea(thresh,sigma,frameNum) = nanmedian([statsPreThresh(keptComponents).Area]);
            end
        end
    end
end
%We calculate the standard deviation of the median area of the mitochondria
%and the mean and standard deviation of the number of mitochondria across
%all temporal stacks for each Gaussian sigma and absolute threshold value.
meanArea = nanmean(medianArea,3);
stdNum = nanstd(numComponents,[],3);
stdArea = nanstd(medianArea,[],3);
stdArea(isnan(stdArea)) = nanmean(stdArea(:));
meanNum = nanmean(numComponents,3);

%We take the z-score of each of these mitochondrial parameters to build up
%our cost matrix.
normStdArea = zscore(stdArea,0,'all');
normStdNum = zscore(stdNum,0,'all');
normMeanNum = zscore(meanNum,0,'all');

% normStdArea(meanArea<minArea) = 1e4; normStdArea(meanArea>maxArea) = 1e4;

%To ensure choosing a stable value, we run the cost matrix through a
%symmetric 3x3 median filter, which acts to remove any outlying regions.
costMatrix = medfilt2(normStdArea + normStdNum - 0.5*normMeanNum,'symmetric');
costMatrix(costMatrix==0) = 1e4;

%We then select the minimum value of this cost matrix to set our sigma and
%threshold values.
minVal = min(costMatrix,[],'all');

[rowMin,colMin] = find(costMatrix == minVal);

threshOptimal = threshMatrix(floor(median(rowMin)));
sigmaOptimal = sigmaMatrix(floor(median(colMin)));

costs.stdArea = stdArea;
costs.stdNum = stdNum;
costs.meanNum = meanNum;
costs.normStdArea = normStdArea;
costs.normStdNum = normStdNum;
costs.normMeanNum = normMeanNum;
costs.costMatrix = costMatrix;
costs.numComponents = numComponents;
costs.medianArea = medianArea;
costs.rowMin = rowMin;
costs.colMin = colMin;


% figure imagesc(costMatrix) axis image colormap hot hold on
% scatter(colMin,rowMin,'r') figure imagesc(normStdArea) axis image
% colormap hot hold on scatter(colMin,rowMin,'r') figure
% imagesc(normStdNum) axis image colormap hot hold on
% scatter(colMin,rowMin,'r') figure imagesc(-normMeanNum) axis image
% colormap hot hold on scatter(colMin,rowMin,'r')

delete(loading)

end
