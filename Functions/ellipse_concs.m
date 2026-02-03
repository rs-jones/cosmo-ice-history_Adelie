function [sd1,sd2] = ellipse_concs(NBe, delNBe, NAl, delNAl, plotSelect, plotString)
% Produces uncertainty ellipses for [Be-10] vs [Al-26] concentration plots.
%
% Syntax:
%   [h1,h2] = ellipse_concentration(NBe,delNBe,NAl,delNAl,2,plotString)
%   [h1] = ellipse_concentration(NBe,delNBe,NAl,delNAl,1,plotString)
%   [x,y] = ellipse_concentration(NBe,delNBe,NAl,delNAl,0,plotString)
%
% Inputs:
%   NBe       - 10Be concentration
%   delNBe    - Uncertainty in 10Be
%   NAl       - 26Al concentration
%   delNAl    - Uncertainty in 26Al
%   plotSelect - 0 = return (x,y) only, 1 = plot 68%, 2 = plot 68% & 95%
%   plotString - Line style/colour string (e.g., 'r', '--k')
%
% Outputs:
%   sd1 - Handle to 68% ellipse (or x-values if plotSelect==0)
%   sd2 - Handle to 95% ellipse (or y-values if plotSelect==0)
%
% Author: Adapted from Balco's Lal-Klein-Nishiizumi ellipse function

% Check arguments
if nargin < 5
    plotSelect = 1;
end
if nargin < 6
    plotString = 'k';
end
if plotSelect ~= 0 && plotSelect ~= 1 && plotSelect ~= 2
    plotSelect = 1;
end

% Estimate range and create meshgrid
[x, y] = meshgrid(...
    (NBe - 4 * delNBe):0.1 * delNBe:(NBe + 4 * delNBe), ...
    (NAl - 4 * delNAl):0.1 * delNAl:(NAl + 4 * delNAl));

% Compute 2D uncorrelated Gaussian PDF
Prob = exp(-0.5 * (((x - NBe) ./ delNBe).^2 + ((y - NAl) ./ delNAl).^2));

% Normalize PDF to total volume = 1
normP = Prob ./ sum(Prob(:));

% Convert to integer representation for cumulative probability
normP = normP * 10000;
intP = round(normP);

maxval = max(intP(:));
cumprob = zeros(1, maxval);
for a = 1:maxval
    cumprob(a) = sum(intP(intP >= a)) / 10000;
end

% Find thresholds for 68% and 95% confidence levels
sigma1 = find(abs(cumprob - 0.68) == min(abs(cumprob - 0.68)), 1);
sigma2 = find(abs(cumprob - 0.95) == min(abs(cumprob - 0.95)), 1);

% Compute contours
cmat = contourc(x(1,:), y(:,1), normP, [sigma1 sigma2]);

% Extract 68% ellipse (first level)
contourStarts = find(cmat(1,:) == sigma1);
contourSizes = cmat(2, contourStarts);
contourToPlot = find(contourSizes == max(contourSizes));

x1 = cmat(1, (contourStarts(contourToPlot)+1):(contourStarts(contourToPlot) + contourSizes(contourToPlot)));
y1 = cmat(2, (contourStarts(contourToPlot)+1):(contourStarts(contourToPlot) + contourSizes(contourToPlot)));

% Output
if plotSelect == 0
    sd1 = x1;
    sd2 = y1;
elseif plotSelect == 1
    sd1 = plot(x1, y1, plotString);
    sd2 = [];
elseif plotSelect == 2
    hold on
    sd1 = plot(x1, y1, plotString);
    
    % Extract 95% ellipse
    contourStarts = find(cmat(1,:) == sigma2);
    contourSizes = cmat(2, contourStarts);
    contourToPlot = find(contourSizes == max(contourSizes));

    x2 = cmat(1, (contourStarts(contourToPlot)+1):(contourStarts(contourToPlot) + contourSizes(contourToPlot)));
    y2 = cmat(2, (contourStarts(contourToPlot)+1):(contourStarts(contourToPlot) + contourSizes(contourToPlot)));

    sd2 = plot(x2, y2, plotString);
end
end
