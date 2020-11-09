function [fissionMatrix, extra] = checkFissionArea(newTrackNum,track,stdRange,CVperc)
numTracks = length(track);

%Initialize
preArea = inf(1,numTracks);
postArea = inf(1,numTracks);
newArea = inf(1,numTracks);
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
    
    prefission = track(trackNum).Area(1:index(end));
    postfission = track(trackNum).Area(index(end)+1:end);
    
    if nnz(~track(trackNum).confident(1:index(end))) || nnz(~track(trackNum).confident(index(end)+1:end)) || length(prefission) < 2 || length(postfission) < 2
        continue
    end
    
    preArea(trackNum) = nanmean(track(trackNum).Area(1:index(end)));
    postArea(trackNum) = nanmean(track(trackNum).Area(index(end)+1:end));
    newArea(trackNum) = nanmean(newTrack.Area);
    
    preStd(trackNum) = nanstd(track(trackNum).Area(1:index(end)));
    postStd(trackNum) = nanstd(track(trackNum).Area(index(end)+1:end));
    newStd(trackNum) = nanstd(newTrack.Area);
end

areaDiff = abs(preArea-(postArea+newArea));
stdTotal = sqrt(preStd.^2+postStd.^2+newStd.^2);
% coeffVarMatrix = stdTotal./areaDiff;
coeffVarMatrixPost = postStd./postArea;
coeffVarMatrixPre = preStd./preArea;
coeffVarMatrixNew = newStd./newArea;
coeffVarMatrix = coeffVarMatrixPost.*coeffVarMatrixPre.*coeffVarMatrixNew;
coeffVarMatrix(isnan(coeffVarMatrix)) = Inf;
% coeffVarMatrix(coeffVarMatrix == 0) = Inf;
possibleFissionMatrix = areaDiff<(stdTotal*stdRange);
fissionMatrix = possibleFissionMatrix.*(coeffVarMatrix<CVperc);
fissionMatrix(isnan(fissionMatrix)) = 0;

extra.preArea = preArea;
extra.postArea = postArea;
extra.newArea = newArea;
extra.preStd = preStd;
extra.postStd = postStd;
extra.newStd = newStd;
extra.areaDiff = areaDiff;
extra.stdTotal = stdTotal;
extra.coeffVarMatrixPost = coeffVarMatrixPost;
extra.coeffVarMatrixPre = coeffVarMatrixPre;
extra.coeffVarMatrixNew = coeffVarMatrixNew;
extra.coeffVarMatrix = coeffVarMatrix;
extra.possibleFissionMatrix = possibleFissionMatrix;
extra.fissionMatrix = fissionMatrix;


end