function [weights,allPerfectTracks,perfectTrackIdx] = getTrackingWeights(track)

perfectTrack = zeros(1,length(track));

for trackNum = 1:length(track)
    if sum(track(trackNum).confident)==length(track(trackNum).confident) && length(track(trackNum).frame)>1
        perfectTrack(trackNum) = 1;
    end
end
perfectTrackIdx = find(perfectTrack);
allPerfectTracks = track(find(perfectTrack));

%Initialize
meanArea = zeros(1,length(allPerfectTracks));
meanMajorAxis = zeros(1,length(allPerfectTracks));
meanMinorAxis = zeros(1,length(allPerfectTracks));
meanSolidity = zeros(1,length(allPerfectTracks));
meanPeri = zeros(1,length(allPerfectTracks));
meanMeanIntensity = zeros(1,length(allPerfectTracks));

intraAreaCoV = zeros(1,length(allPerfectTracks));
intraMajAxCoV = zeros(1,length(allPerfectTracks));
intraMinAxCoV = zeros(1,length(allPerfectTracks));
intraSolCoV = zeros(1,length(allPerfectTracks));
intraPeriCoV = zeros(1,length(allPerfectTracks));
intraMeanIntCoV = zeros(1,length(allPerfectTracks));


for trackNum = 1:length(allPerfectTracks)
    meanArea(trackNum) = mean(allPerfectTracks(trackNum).Area);
    meanMajorAxis(trackNum) = mean(allPerfectTracks(trackNum).MajorAxisLength);
    meanMinorAxis(trackNum) = mean(allPerfectTracks(trackNum).MinorAxisLength);
    meanSolidity(trackNum) = mean(allPerfectTracks(trackNum).Solidity);
    meanPeri(trackNum) = mean(allPerfectTracks(trackNum).Perimeter);
    meanMeanIntensity(trackNum) = mean(allPerfectTracks(trackNum).MeanIntensity);
    
    intraAreaCoV(trackNum) = std(allPerfectTracks(trackNum).Area)/mean(allPerfectTracks(trackNum).Area);
    intraMajAxCoV(trackNum) = std(allPerfectTracks(trackNum).MajorAxisLength)/mean(allPerfectTracks(trackNum).MajorAxisLength);
    intraMinAxCoV(trackNum) = std(allPerfectTracks(trackNum).MinorAxisLength)/mean(allPerfectTracks(trackNum).MinorAxisLength);
    intraSolCoV(trackNum) = std(allPerfectTracks(trackNum).Solidity)/mean(allPerfectTracks(trackNum).Solidity);
    intraPeriCoV(trackNum) = std(allPerfectTracks(trackNum).Perimeter)/mean(allPerfectTracks(trackNum).Perimeter);
    intraMeanIntCoV(trackNum) = std(allPerfectTracks(trackNum).MeanIntensity)/mean(allPerfectTracks(trackNum).MeanIntensity);
end

interAreaCoV = nanstd(meanArea(:))/nanmean(meanArea(:));
interMajAxCoV = nanstd(meanMajorAxis(:))/nanmean(meanMajorAxis(:));
interMinAxCoV = nanstd(meanMinorAxis(:))/nanmean(meanMinorAxis(:));
interSolCoV = nanstd(meanSolidity(:))/nanmean(meanSolidity(:));
interPeriCoV = nanstd(meanPeri(:))/nanmean(meanPeri(:));
interMeanIntCoV = nanstd(meanMeanIntensity(:))/nanmean(meanMeanIntensity(:));

avIntraAreaCoV = nanmean(intraAreaCoV);
avIntraMajAxCoV = nanmean(intraMajAxCoV);
avIntraMinAxCoV = nanmean(intraMinAxCoV);
avIntraSolCoV = nanmean(intraSolCoV);
avIntraPeriCoV = nanmean(intraPeriCoV);
avIntraMeanIntCoV = nanmean(intraMeanIntCoV);


areaWeight = interAreaCoV/avIntraAreaCoV;
majAxWeight = interMajAxCoV/avIntraMajAxCoV;
minAxWeight = interMinAxCoV/avIntraMinAxCoV;
solWeight = interSolCoV/avIntraSolCoV;
periWeight = interPeriCoV/avIntraPeriCoV;
meanIntWeight = interMeanIntCoV/avIntraMeanIntCoV;

sumWeights = areaWeight+majAxWeight+minAxWeight+solWeight+periWeight+meanIntWeight;

areaWeightNorm = areaWeight/sumWeights;
majAxWeightNorm = majAxWeight/sumWeights;
minAxWeightNorm = minAxWeight/sumWeights;
solWeightNorm = solWeight/sumWeights;
periWeightNorm = periWeight/sumWeights;
meanIntWeightNorm = meanIntWeight/sumWeights;

weights = [areaWeightNorm,majAxWeightNorm,minAxWeightNorm,solWeightNorm,periWeightNorm,meanIntWeightNorm];
% 
% if nnz(~weights)==6
%     weights(:) = 1/6;
% end


end