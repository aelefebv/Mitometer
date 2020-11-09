function [fusionMatrix,extra] = checkFusionVolume(lostTrackNum,track,stdRange,CVperc)
numTracks = length(track);

%Initialize
preVolume = inf(1,numTracks);
postVolume = inf(1,numTracks);
lostVolume = inf(1,numTracks);
preStd = inf(1,numTracks);
postStd = inf(1,numTracks);
lostStd = inf(1,numTracks);

lostTrack = track(lostTrackNum);

for trackNum = 1:numTracks
    frameMin = lostTrack.frame(end);
    if sum(track(trackNum).frame>frameMin) && sum(track(trackNum).frame<=frameMin)
        index = find(track(trackNum).frame>frameMin);
    else
        continue
    end
    
    prefusion = track(trackNum).Volume(1:index(1)-1);
    postfusion = track(trackNum).Volume(index(1):end);
    
    if nnz(~track(trackNum).confident(1:index(1)-1)) || nnz(~track(trackNum).confident(index(1):end)) || length(prefusion) < 2 || length(postfusion) < 2
        continue
    end
    preVolume(trackNum) = nanmean(prefusion);
    postVolume(trackNum) = nanmean(postfusion);
    lostVolume(trackNum) = nanmean(lostTrack.Volume);
    
    preStd(trackNum) = nanstd(track(trackNum).Volume(1:index(1)-1));
    postStd(trackNum) = nanstd(track(trackNum).Volume(index(1):end));
    lostStd(trackNum) = nanstd(lostTrack.Volume);
end

volumeDiff = abs(postVolume-(preVolume+lostVolume));
stdTotal = sqrt(preStd.^2+postStd.^2+lostStd.^2);
% coeffVarMatrix = stdTotal./volumeDiff;
coeffVarMatrixPost = postStd./postVolume;
coeffVarMatrixPre = preStd./preVolume;
coeffVarMatrixLost = lostStd./lostVolume;
coeffVarMatrix = coeffVarMatrixPost.*coeffVarMatrixPre.*coeffVarMatrixLost;
coeffVarMatrix(isnan(coeffVarMatrix)) = Inf;
% coeffVarMatrix(coeffVarMatrix == 0) = Inf;
possibleFusionMatrix = volumeDiff<(stdTotal*stdRange);
fusionMatrix = possibleFusionMatrix.*(coeffVarMatrix<CVperc);
fusionMatrix(isnan(fusionMatrix)) = 0;

extra.preVolume = preVolume;
extra.postVolume = postVolume;
extra.lostVolume = lostVolume;
extra.preStd = preStd;
extra.postStd = postStd;
extra.lostStd = lostStd;
extra.volumeDiff = volumeDiff;
extra.stdTotal = stdTotal;
extra.coeffVarMatrixPost = coeffVarMatrixPost;
extra.coeffVarMatrixPre = coeffVarMatrixPre;
extra.coeffVarMatrixLost = coeffVarMatrixLost;
extra.coeffVarMatrix = coeffVarMatrix;
extra.possibleFusionMatrix = possibleFusionMatrix;
extra.fusionMatrix = fusionMatrix;

end