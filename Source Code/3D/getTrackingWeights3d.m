function [weights,allPerfectTracks,perfectTrackIdx] = getTrackingWeights3d(track)
%The weighting of each mitochondrial morphological feature depends on the
%coefficient of variation (CoV) within optimal tracks and between optimal
%tracks. If a morphological feature changes greatly within the optimal
%tracks, it will have a lesser weighting, whereas if it varies greatly
%between optimal tracks, it will have a greater weighting.

perfectTrack = zeros(1,length(track));

for trackNum = 1:length(track)
    if sum(track(trackNum).confident)==length(track(trackNum).confident) && length(track(trackNum).frame)>1
        perfectTrack(trackNum) = 1;
    end
end
perfectTrackIdx = find(perfectTrack);
allPerfectTracks = track(find(perfectTrack));

%Initialize
meanVolume = zeros(1,length(allPerfectTracks));
meanMajorAxis = zeros(1,length(allPerfectTracks));
meanMinorAxis = zeros(1,length(allPerfectTracks));
meanZAxis = zeros(1,length(allPerfectTracks));
meanSolidity = zeros(1,length(allPerfectTracks));
meanSurf = zeros(1,length(allPerfectTracks));
meanMeanIntensity = zeros(1,length(allPerfectTracks));

intraVolumeCoV = zeros(1,length(allPerfectTracks));
intraMajAxCoV = zeros(1,length(allPerfectTracks));
intraMinAxCoV = zeros(1,length(allPerfectTracks));
intraZAxCoV = zeros(1,length(allPerfectTracks));
intraSolCoV = zeros(1,length(allPerfectTracks));
intraSurfCoV = zeros(1,length(allPerfectTracks));
intraMeanIntCoV = zeros(1,length(allPerfectTracks));


for trackNum = 1:length(allPerfectTracks)
    %We build a vector of mean values for each feature of each optimal
    %track.
    meanVolume(trackNum) = mean(allPerfectTracks(trackNum).Volume);
    meanMajorAxis(trackNum) = mean(allPerfectTracks(trackNum).MajorAxisLength);
    meanMinorAxis(trackNum) = mean(allPerfectTracks(trackNum).MinorAxisLength);
    meanZAxis(trackNum) = mean(allPerfectTracks(trackNum).ZAxisLength);
    meanSolidity(trackNum) = mean(allPerfectTracks(trackNum).Solidity);
    meanSurf(trackNum) = mean(allPerfectTracks(trackNum).SurfaceArea);
    meanMeanIntensity(trackNum) = mean(allPerfectTracks(trackNum).MeanIntensity);
    
    % For each individual optimal track, we calculate the CoV of each
    % feature within that track, the intraCoV.
    intraVolumeCoV(trackNum) = std(allPerfectTracks(trackNum).Volume)/mean(allPerfectTracks(trackNum).Volume);
    intraMajAxCoV(trackNum) = std(allPerfectTracks(trackNum).MajorAxisLength)/mean(allPerfectTracks(trackNum).MajorAxisLength);
    intraMinAxCoV(trackNum) = std(allPerfectTracks(trackNum).MinorAxisLength)/mean(allPerfectTracks(trackNum).MinorAxisLength);
    intraZAxCoV(trackNum) = std(allPerfectTracks(trackNum).ZAxisLength)/mean(allPerfectTracks(trackNum).ZAxisLength);
    intraSolCoV(trackNum) = std(allPerfectTracks(trackNum).Solidity)/mean(allPerfectTracks(trackNum).Solidity);
    intraSurfCoV(trackNum) = std(allPerfectTracks(trackNum).SurfaceArea)/mean(allPerfectTracks(trackNum).SurfaceArea);
    intraMeanIntCoV(trackNum) = std(allPerfectTracks(trackNum).MeanIntensity)/mean(allPerfectTracks(trackNum).MeanIntensity);
end

%We then calculate the CoV of each feature for all tracks, the interCoV
interVolumeCoV = std(meanVolume(:))/mean(meanVolume(:));
interMajAxCoV = std(meanMajorAxis(:))/mean(meanMajorAxis(:));
interMinAxCoV = std(meanMinorAxis(:))/mean(meanMinorAxis(:));
interZAxCoV = std(meanZAxis(:))/mean(meanZAxis(:));
interSolCoV = std(meanSolidity(:))/mean(meanSolidity(:));
interSurfCoV = std(meanSurf(:))/mean(meanSurf(:));
interMeanIntCoV = std(meanMeanIntensity(:))/mean(meanMeanIntensity(:));

avIntraVolumeCoV = mean(intraVolumeCoV);
avIntraMajAxCoV = mean(intraMajAxCoV);
avIntraMinAxCoV = mean(intraMinAxCoV);
avIntraZAxCoV = mean(intraZAxCoV);
avIntraSolCoV = mean(intraSolCoV);
avIntraSurfCoV = mean(intraSurfCoV);
avIntraMeanIntCoV = mean(intraMeanIntCoV);

%We then weight each feature as the ratio of the interCoV to the mean
%intraCoV and normalize each weight by the sum of the weights to produce a
%final sum weighting of 1.

volumeWeight = interVolumeCoV/avIntraVolumeCoV;
majAxWeight = interMajAxCoV/avIntraMajAxCoV;
minAxWeight = interMinAxCoV/avIntraMinAxCoV;
zAxWeight = interZAxCoV/avIntraZAxCoV;
solWeight = interSolCoV/avIntraSolCoV;
surfWeight = interSurfCoV/avIntraSurfCoV;
meanIntWeight = interMeanIntCoV/avIntraMeanIntCoV;

sumWeights = volumeWeight+majAxWeight+minAxWeight+zAxWeight+solWeight+surfWeight+meanIntWeight;

volumeWeightNorm = volumeWeight/sumWeights;
majAxWeightNorm = majAxWeight/sumWeights;
minAxWeightNorm = minAxWeight/sumWeights;
zAxWeightNorm = minAxWeight/sumWeights;
solWeightNorm = solWeight/sumWeights;
surfWeightNorm = surfWeight/sumWeights;
meanIntWeightNorm = meanIntWeight/sumWeights;

weights = [volumeWeightNorm,majAxWeightNorm,minAxWeightNorm,zAxWeightNorm,solWeightNorm,surfWeightNorm,meanIntWeightNorm];

end