
% In this script, four values are closely approximated via golden section search
%  - λ value for which correlation coefficient ρ is strongest between S^o and ABC_λ
%  - λ value for which ρ is strongest between S^o and H_λ
%  - λ value for which ρ is strongest between C_p and ABC_λ
%  - λ value for which ρ is strongest between C_p and H_λ
% Additionally, curves for ρ against λ near these values are plotted in 4 figs

close all; % Close any figures already opened
clear;     % and clear all variables
format long; % More significant figures printed in the console
pkg load statistics;

% Config: Adjust lines' width, font size, whether or not to save to file
lineWidth = 2;
fontSize = 16;
saveToFile = true;

% Utility: Below is a function to round-off to 4 decimal places | returns string
as_4_dp_str = @(x) sprintf('%.04f', round(x*(10^4))/(10^4));

% Cell containing Entropy and Heat Capacity of 30 lower benzenoids
expData = {reshape([ % Entropy
    269.722 334.155 389.475 395.882 444.724 447.437
    457.958 455.839 450.418 399.491 499.831 513.857
    508.537 507.395 506.076 512.523 500.734 510.307
    509.210 513.879 511.770 509.611 461.545 463.738
    468.712 555.409 472.295 554.784 468.796 551.708
  ]', 30, 1), reshape([ % Heat Capacity
    83.019 133.325 184.194 183.654 235.165 233.497
   234.568 234.638 233.558 200.815 286.182 285.056
   284.037 284.088 285.148 284.595 284.870 284.503
   284.785 284.740 284.233 284.552 251.175 250.568
   251.973 336.098 267.543 337.204 285.041 368.518
  ]', 30, 1);
  "S^o", "C_p" % <-- CHANGED: was "E", "ΔH"
};

% 30 by 3 array, number of edges in each dudv partition: (2,2), (2,3) then (3,3)
d_f = [
  6 6 6 7  6 8  7 8 9 6  6  7  8  8  7 10 9  9  9  8  9  8 8 8  7 10  7  6  6  6
  0 4 8 6 12 8 10 8 6 8 16 14 12 12 14 8 10 10 10 12 10 12 8 8 10 12 10 20 12 16
  0 1 2 3  3 5  4 5 6 5  4  5  6  6  5 8  7  7  7  6  7  6 8 8  7  9 10  5 12 19
]';

% Cell containing the index-computing functions
% Usage: getIndexFns{n}(λ) | n=1:ABC_λ, n=2:H_λ, l = lambda
getIndexFns = {
  @(lambda) (sum(d_f.*[(1/2),(1/2),(4/9)].^lambda,2)); % Atom-Bond Connectivity (ABC)
  @(lambda) (sum(d_f.*[(1/2),(2/5),(1/3)].^lambda,2)); % Harmonic (H)
}';
indexName = {"ABC" "H"};

% Variables for indexing arrays and iterating for-loops
numData    = size(expData, 2);       % two (S^o & C_p)
numIndices = size(getIndexFns, 2);   % two (ABC & H)
numCurves  = numData*numIndices;     % four

% All x in visible ranges (both plots - near and far)
xnear = [linspace(-0.63, 0.2, 800); linspace(-3.30, 0.25, 800)];
xfar  = linspace(-20, 20, 800);

% Do the same procedure for each experimental data i.e., S^o and C_p
for ii = 1:numData
  figure(numData+ii); hold on;
  ymin_all = inf;
  ymax_all = -inf;

  % WARNING: these xmeet1-ymeet2 values are hardcoded, computed separately
  xmeet1 = [
    -0.46299; % for S^o
    -2.2887   % for C_p
  ](ii);
  ymeet1 = [
    0.99641; % for S^o
    0.97216  % for C_p
  ](ii);
  xmeet2 = [
    -0.00045;    % for S^o
    -0.00010081  % for C_p
  ](ii);
  ymeet2 = [
    0.99574; % for S^o
    0.94512  % for C_p
  ](ii);

  ybox = [
    0.9968; % for figure 3 (blue dotted line)
    0.98    % for figure 4 (blue dotted line)
  ](ii);

  plot([xmeet1 xmeet1 0 0], [0 ybox ybox 0], '--b', 'LineWidth', lineWidth);

  yend = 0;

  % Draw correlation curves for each index-property pair
  for n = 1:numIndices
    get_indices_vals = @(lambda) getIndexFns{n}(lambda)(!isnan(expData{1,ii}));
    ccFn = @(lambda) corrcoef(
      get_indices_vals(lambda),
      expData{1,ii}(!isnan(expData{1,ii}))
    )(1,2);

    ynear = arrayfun(ccFn, xnear(ii,:));
    ymin_all = min(ymin_all, min(ynear));
    ymax_all = max(ymax_all, max(ynear));
    yfar = arrayfun(ccFn, xfar);

    % Compute peak values and display in console
    disp(sprintf("%s against general_%s", expData{2,ii}, indexName{n}));
    peakLambda = mean(GoldenSectionSearch_Maximum(ccFn, xnear(ii,1), xnear(ii,end), 1e-15))
    peakCorrCoeff = ccFn(peakLambda)

    X_lm = getIndexFns{n}(peakLambda);
    Y_lm = expData{1,ii};
    mdl = fitlm(X_lm, Y_lm);
    close

    grad_coef            = mdl{3,2}
    grad_lower           = mdl{3,4};
    CI_grad              = grad_coef - grad_lower
    intercept_coef       = mdl{2,2}
    intercept_lower_CI   = mdl{2,4};
    CI_intercept         = intercept_coef - intercept_lower_CI
    disp("")

    % Generate curve label
    curveLabels{n} = sprintf("%s and %s_l", expData{2,ii}, indexName{n});

    figure(ii);
    hold on;
    curvesFar(n) = plot(xfar, yfar, '-', 'LineWidth', lineWidth);
    drawnow;

    figure(numData+ii);
    hold on;
    curves(n) = plot(xnear(ii,:), ynear, '-', 'LineWidth', lineWidth);

    plot([peakLambda peakLambda xnear(ii,1)], [0 peakCorrCoeff peakCorrCoeff],
         '--k', 'LineWidth', lineWidth/2);
    text(peakLambda, peakCorrCoeff,
        {'', sprintf("(−%s, %s)", as_4_dp_str(abs(peakLambda)), as_4_dp_str(peakCorrCoeff))},
        'VerticalAlignment', 'bottom');

    yend = max(yend, ynear(end));
  end

  % Mark intersection points where ABC_λ is better than H_λ
  plot([xmeet1 xmeet2], [ymeet1 ymeet2], '*b',
       'MarkerSize', 16, 'LineWidth', lineWidth/1.5);
  text(xmeet1, ymeet1, {'', sprintf(" (−%s, %s)", as_4_dp_str(abs(xmeet1)), as_4_dp_str(ymeet1))},
       'VerticalAlignment', 'top', 'Color', [0, 0, 0.8]);
  text(xmeet2, ymeet2, {'', sprintf("(0, %s) ", as_4_dp_str(ymeet2))},
       'VerticalAlignment', 'top', 'HorizontalAlignment', 'right', 'Color', [0, 0, 0.8]);

  % Label zoomed-in plot
  xlabel('λ');
  ylabel('ρ');
  leg = legend(curves, curveLabels);
  set(leg, 'location', "southwest");

  ybox_space = [
    0.000005;
    0.0005
  ](ii);
  axis([xnear(ii,1) xnear(ii,end) ymin_all-0.0005 ymax_all+0.0005]);
  drawnow;

  % Label zoomed-out plot
  figure(ii);
  xlabel('λ');
  ylabel('ρ');
  leg = legend(curvesFar, curveLabels);
  set(leg, 'location', "southeast");

  if ii == 2
    set(leg, 'location', "southeast");
  end

  hold off;
end

for ii = 1:4
  figure(ii);
  xticklabels(strrep(xticklabels, '-', '−'));
  yticklabels(strrep(yticklabels, '-', '−'));
  set(findall(figure(ii), '-property', 'FontSize'), 'FontSize', fontSize)
end

if saveToFile
  saveas(figure(1), "01_comb_ccurves_So_indices_FAR.png");
  saveas(figure(2), "01_comb_ccurves_Cp_indices_FAR.png");
  saveas(figure(3), "01_comb_ccurves_So_indices_NEAR.png");
  saveas(figure(4), "01_comb_ccurves_Cp_indices_NEAR.png");
end
