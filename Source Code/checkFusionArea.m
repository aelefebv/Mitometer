function [fusionMatrix,extra] = checkFusionArea(lostTrackNum,track,stdRange,CVperc)
numTracks = length(track);

%Initialize
preArea = inf(1,numTracks);
postArea = inf(1,numTracks);
lostArea = inf(1,numTracks);
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
    
    prefusion = track(trackNum).Area(1:index(1)-1);
    postfusion = track(trackNum).Area(index(1):end);
    
    if nnz(~track(trackNum).confident(1:index(1)-1)) || nnz(~track(trackNum).confident(index(1):end)) || length(prefusion) < 2 || length(postfusion) < 2
        continue
    end
    preArea(trackNum) = nanmean(prefusion);
    postArea(trackNum) = nanmean(postfusion);
    lostArea(trackNum) = nanmean(lostTrack.Area);
    
    preStd(trackNum) = nanstd(track(trackNum).Area(1:index(1)-1));
    postStd(trackNum) = nanstd(track(trackNum).Area(index(1):end));
    lostStd(trackNum) = nanstd(lostTrack.Area);
end

areaDiff = abs(postArea-(preArea+lostArea));
stdTotal = sqrt(preStd.^2+postStd.^2+lostStd.^2);
% coeffVarMatrix = stdTotal./areaDiff;
coeffVarMatrixPost = postStd./postArea;
coeffVarMatrixPre = preStd./preArea;
coeffVarMatrixLost = lostStd./lostArea;
coeffVarMatrix = coeffVarMatrixPost.*coeffVarMatrixPre.*coeffVarMatrixLost;
coeffVarMatrix(isnan(coeffVarMatrix)) = Inf;
% coeffVarMatrix(coeffVarMatrix == 0) = Inf;
possibleFusionMatrix = areaDiff<(stdTotal*stdRange);
fusionMatrix = possibleFusionMatrix.*(coeffVarMatrix<CVperc);
fusionMatrix(isnan(fusionMatrix)) = 0;

extra.preArea = preArea;
extra.postArea = postArea;
extra.lostArea = lostArea;
extra.preStd = preStd;
extra.postStd = postStd;
extra.lostStd = lostStd;
extra.areaDiff = areaDiff;
extra.stdTotal = stdTotal;
extra.coeffVarMatrixPost = coeffVarMatrixPost;
extra.coeffVarMatrixPre = coeffVarMatrixPre;
extra.coeffVarMatrixLost = coeffVarMatrixLost;
extra.coeffVarMatrix = coeffVarMatrix;
extra.possibleFusionMatrix = possibleFusionMatrix;
extra.fusionMatrix = fusionMatrix;

end