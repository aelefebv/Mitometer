function ImMask = thresholdImage(ImGaussianFiltered,thr)

ImMask = zeros(size(ImGaussianFiltered),class(ImGaussianFiltered));
ImMask(ImGaussianFiltered>thr) = 1;
ImMask(ImGaussianFiltered<=thr) = 0;

end