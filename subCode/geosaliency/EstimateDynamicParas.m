function [clipVal, geoSigma, neiSigma] = EstimateDynamicParas(adjcMatrix, colDistM)
[meanMin1, meanTop, meanMin2] = GetMeanMinAndMeanTop(adjcMatrix, colDistM, 0.01);
clipVal = meanMin2;
% Emperically choose adaptive sigma for converting geodesic distance to
% weight
geoSigma = min([10, meanMin1 * 3, meanTop / 10]);
geoSigma = max(geoSigma, 5);

% Emperically choose adaptive sigma for smoothness term in Equa(9) of our
% paper.
neiSigma = min([3 * meanMin1, 0.2 * meanTop, 20]);