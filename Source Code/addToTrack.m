function outtrack = addToTrack(mito,track)
outtrack = struct;  %final structure
for field = fieldnames(mito)'
    fname = field{1};
    if strcmp((field),'PixelIdxList') && ~iscell(mito.(fname))
        outtrack.(fname) = horzcat(track.(fname),{mito.(fname)});
    else
        outtrack.(fname) = horzcat(track.(fname),mito.(fname));
    end
end
end