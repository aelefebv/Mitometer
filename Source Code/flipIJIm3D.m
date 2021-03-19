function imFlip = flipIJIm3D(im)

imFlip = permute(im,[2,1,4,3]);

end