function NNExtrema = getNNExtrema(mito)
numMito = length(mito);

%initialize
extremaMatrix = zeros(numMito,numMito);

for mitoNum = 1:numMito
    for mitoNum2 = mitoNum:numMito
        extremaMatrix(mitoNum,mitoNum2) = min(sqrt( (mito(mitoNum).Extrema(:,1)-mito(mitoNum2).Extrema(:,1)).^2 + (mito(mitoNum).Extrema(:,2)-mito(mitoNum2).Extrema(:,2)).^2 ));
    end
end

extremaMatrix(extremaMatrix==0) = Inf;

NNExtremaRow = min(extremaMatrix,[],2);
NNExtremaCol = min(extremaMatrix,[],1);

NNExtrema = min([NNExtremaRow,NNExtremaCol'],[],2);

end