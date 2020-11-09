function [fissionMatrix, extra] = checkFissionVolume(newTrackNum,track,stdRange,CVperc)
numTracks = length(track);

%Initialize
preVolume = inf(1,numTracks);
postVolume = inf(1,numTracks);
newVolume = inf(1,numTracks);
preStd = inf(1,numTracks);
postStd = inf(1,numTracks);
newStd = inf(1,numTracks);

newTrack = track(newTrackNum);

for trackNum = 1:numTracks
    frameMax = newTrack.frame(1);
    if sum(track(trackNum).frame<frameMax) && sum(track(trackNum).frame>=frameMax)
        index = find(track(trackNum).frame<frameMax);
    else
        continue
    end
    
    prefission = track(trackNum).Volume(1:index(end));
    postfission = track(trackNum).Volume(index(end)+1:end);
    
    if nnz(~track(trackNum).confident(1:index(end))) || nnz(~track(trackNum).confident(index(end)+1:end)) || length(prefission) < 2 || length(postfission) < 2
        continue
    end
    preVolume(trackNum) = nanmean(track(trackNum).Volume(1:index(end)));
    postVolume(trackNum) = nanmean(track(trackNum).Volume(index(end)+1:end));
    newVolume(trackNum) = nanmean(newTrack.Volume);
    
    preStd(trackNum) = nanstd(track(trackNum).Volume(1:index(end)));
    postStd(trackNum) = nanstd(track(trackNum).Volume(index(end)+1:end));
    newStd(trackNum) = nanstd(newTrack.Volume);
end

volumeDiff = abs(preVolume-(postVolume+newVolume));
stdTotal = sqrt(preStd.^2+postStd.^2+newStd.^2);
% coeffVarMatrix = stdTotal./volumeDiff;
coeffVarMatrixPost = postStd./postVolume;
coeffVarMatrixPre = preStd./preVolume;
coeffVarMatrixNew = newStd./newVolume;
coeffVarMatrix = coeffVarMatrixPost.*coeffVarMatrixPre.*coeffVarMatrixNew;
coeffVarMatrix(isnan(coeffVarMatrix)) = Inf;
% coeffVarMatrix(coeffVarMatrix == 0) = Inf;
possibleFissionMatrix = volumeDiff<(stdTotal*stdRange);
fissionMatrix = possibleFissionMatrix.*(coeffVarMatrix<CVperc);
fissionMatrix(isnan(fissionMatrix)) = 0;

extra.preVolume = preVolume;
extra.postVolume = postVolume;
extra.newVolume = newVolume;
extra.preStd = preStd;
extra.postStd = postStd;
extra.newStd = newStd;
extra.volumeDiff = volumeDiff;
extra.stdTotal = stdTotal;
extra.coeffVarMatrixPost = coeffVarMatrixPost;
extra.coeffVarMatrixPre = coeffVarMatrixPre;
extra.coeffVarMatrixNew = coeffVarMatrixNew;
extra.coeffVarMatrix = coeffVarMatrix;
extra.possibleFissionMatrix = possibleFissionMatrix;
extra.fissionMatrix = fissionMatrix;


end