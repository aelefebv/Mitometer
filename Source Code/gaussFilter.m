function ImGaussFiltered = gaussFilter(Im,sigma)
%parameters sigma to vary

if sigma <= 0
    sigma = 1E-9;
end

% filterSize = 2*ceil(2*sigma)+1;

gaussFilterMatrix = fspecial('Gaussian',[3 3],sigma);
ImGaussFiltered = imfilter(Im,gaussFilterMatrix);

end