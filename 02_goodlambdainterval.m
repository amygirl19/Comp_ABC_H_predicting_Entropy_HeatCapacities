% In this script, four plots are generated which shows λ intervals for
%  - good correlation coefficient ρ between S^o and ABC_λ
%  - good correlation coefficient ρ between S^o and H_λ
%  - good correlation coefficient ρ between Cp and ABC_λ
%  - good correlation coefficient ρ between Cp and H_λ
% The limits of the intervals are also printed to console

close all; % Close any figures already opened
clear;     % and clear all variables
format long; % More significant figures printed in the console

% Config: Adjust lines' width, font size, whether or not to save to file
lineWidth = 3;
fontSize = 24;
saveToFile = false; % Set to true to auto-save plots

% Utility: Below is a function to round-off to 4 decimal places | returns string
as_4_dp_str = @(x) sprintf('%.04f', round(x*(10^4))/(10^4));

% Cell containing S^o and Cp of 30 benzenoids
entropyData = reshape([
  269.722 334.155 389.475 395.882 444.724 447.437
  457.958 455.839 450.418 399.491 499.831 513.857
  508.537 507.395 506.076 512.523 500.734 510.307
  509.210 513.879 511.770 509.611 461.545 463.738
  468.712 555.409 472.295 554.784 468.796 551.708
]', 30, 1);

enthalpyData = reshape([
  83.019  133.325 184.194 183.654 235.165 233.497
  234.568 234.638 233.558 200.815 286.182 285.056
  284.037 284.088 285.148 284.595 284.870 284.503
  284.785 284.740 284.233 284.552 251.175 250.568
  251.973 336.098 267.543 337.204 285.041 368.518
]', 30, 1);

expData = {entropyData, enthalpyData, "S^o", "C_p"};

% 30 by 3 array, number of edges in each dudv partition: (2,2), (2,3) then (3,3)
d_f = [
  6  6  6  7  6  8  7  8  9  6  6  7  8  8  7  10 9  9  9  8  9  8  8  8  7  10 7  6  6  6
  0  4  8  6  12 8  10 8  6  8  16 14 12 12 14 8  10 10 10 12 10 12 8  8  10 12 10 20 12 16
  0  1  2  3  3  5  4  5  6  5  4  5  6  6  5  8  7  7  7  6  7  6  8  8  7  9  10 5  12 19
]'; % Used for computing indices based on edge endpoint degree partitions

% Cell containing the index-computing functions
% Usage: getIndexFns{n}(lambda) | n=1:ABC_lambda, n=2:H_lambda
getIndexFns = {
  @(lambda) (sum(d_f .* [(1/2),(1/2),(4/9)].^lambda, 2));  % General ABC index
  @(lambda) (sum(d_f .* [(1/2),(2/5),(1/3)].^lambda, 2));  % General H index
}';
% Cell containing their respective labels
indexName = {"ABC" "H"};

% Variables for indexing arrays and iterating for-loops
numData    = 2;                   % S^o and Cp
numIndices = size(getIndexFns,2); % two
numCurves  = numData*numIndices;  % four

% Boundaries for visible intervals, for each index-property pair
%              ABC_λ   H_λ
xstart = [    -15     -5        % S^o
              -10     -5    ];  % Cp
xend = [       40      8        % S^o
               25      8    ];  % Cp
ystart = [     1.0     1.0;     % S^o
               1.0     1.0];    % Cp
yend = [       0.95    0.90;    % S^o
               0.90    0.90 ];  % Cp

% Exact rho value considered good for
%            S^o    Cp
a_goodrho = [0.965; 0.9675];

% Colors (different shades of cyan and green)
colShaded        = {[0.85, 1, 1];      [0.85, 1, 0.85]};
colIndicatorVert = {[0.2, 0.55, 0.55]; [0, 0.5, 0]};
colIndicatorHorz = {[0.35, 0.75, 0.75]; [0.45, 0.7, 0.45]};
colCurve         = {[0, 0.75, 0.75];   [0, 0.5, 0]};

% Do the same procedure for each experimental data i.e., S^o and Cp
for ii = 1:numData
  expVals = expData{1, ii};
  validIdx = !isnan(expVals);

  % Draw correlation curves for each index-property pair
  for n = 1:numIndices
    % Function to get corrcoef ρ between S^o/Cp (depending on ii) with specified λ
    %                                and ABC_λ/H_λ (depending on n)
    ccFn = @(lambda) corrcoef(
      getIndexFns{n}(lambda)(validIdx),
      expVals(validIdx)
    )(1,2);

    this_fig = figure((ii-1)*numIndices+n); % Figures (1) to (4)
    hold on;

    % Get lambda for highest rho first
    peakLambda = mean(
      GoldenSectionSearch_Maximum(ccFn, -15, 40, 1e-15));

    % Prepare a function to calc |rho(lambda)-goodrho|
    ccFn_good = @(lambda)(-abs(ccFn(lambda)-a_goodrho(ii)));

    % Get the limit i.e., value of lambda where rho equals goodrho
    getLimitFromInterval = @(lb, ub) mean(
      GoldenSectionSearch_Maximum(ccFn_good, lb, ub, 1e-15));

    lambda_lb = getLimitFromInterval(peakLambda-25, peakLambda); % Search to the left
    lambda_ub = getLimitFromInterval(peakLambda, peakLambda+30); % Search to the right

    % Write the intervals in console
    disp(sprintf("ρ(%s,%s_λ) ≥ %.02f when λ ∈ [%.08f, %.08f]",
         expData{2+ii}, indexName{n}, a_goodrho(ii), lambda_lb, lambda_ub));

    % STEP 1: Shade the area inside the good lambda interval FIRST
    u0 = lambda_lb;    u_width  = lambda_ub - lambda_lb;
    v0 = ystart(ii,n); v_height = yend(ii,n) - ystart(ii,n);
    rectangle('Position', [u0, v0, u_width, v_height], 'FaceColor', colShaded{ii}, 'LineStyle', 'none');

    % STEP 2: Plot the full blue curve across the entire x range (on top of shading)
    x_full = linspace(xstart(ii,n), xend(ii,n), 600);
    y_full = arrayfun(ccFn, x_full);
    plot(x_full, y_full, '-', 'LineWidth', lineWidth);

    % STEP 3: Draw the indicator lines for the good lambda interval
    % Vertical dashed lines:
    plot([lambda_lb lambda_lb], [ystart(ii,n) yend(ii,n)],
         '--', 'LineWidth', lineWidth/1.75, 'Color', colIndicatorVert{ii});
    plot([lambda_ub lambda_ub], [ystart(ii,n) yend(ii,n)],
         '--', 'LineWidth', lineWidth/1.75, 'Color', colIndicatorVert{ii});
    % Horizontal black dashed line, horizontal colored dashed line:
    plot([xstart(ii,n) lambda_lb], [a_goodrho(ii) a_goodrho(ii)],
         '--k', 'LineWidth', lineWidth/1.75);
    plot([lambda_lb lambda_ub], [a_goodrho(ii) a_goodrho(ii)],
         '--', 'LineWidth', lineWidth/1.75, 'Color', colIndicatorHorz{ii});

    % Write on the plot the interval limits
    text(lambda_lb, yend(ii,n), {'', sprintf("  l=-%s", as_4_dp_str(abs(lambda_lb)))}, 'VerticalAlignment', 'top', 'Rotation', 90);
    text(lambda_ub, yend(ii,n), {'', sprintf("  l=-%s", as_4_dp_str(abs(lambda_ub)))}, 'VerticalAlignment', 'top', 'Rotation', 90);

    % STEP 4: Plot the colored curve within the good lambda range (on top of everything)
    x_in = linspace(lambda_lb, lambda_ub, 400);
    y_in = arrayfun(ccFn, x_in);
    plot(x_in, y_in, '-', 'LineWidth', lineWidth, 'Color', colCurve{ii});
    % Also highlight the good interval on the x-axis
    plot([lambda_lb, lambda_ub], [yend(ii,n), yend(ii,n)], '-',
         'LineWidth', lineWidth, 'Color', colCurve{ii});

    % Label the plot
    title(sprintf('between %s and %s_l', expData{2+ii}, indexName{n}));
    xlabel('λ');
    ylabel('ρ');
    drawnow;
    axis([xstart(ii,n) xend(ii,n) yend(ii,n) ystart(ii,n)]);

    % Replace all hyphens with minuses
    xticklabels(strrep(xticklabels, '-', '−'));
    yticklabels(strrep(yticklabels, '-', '−'));
    % Set all fontsizes to size specified early in the script
    set(findall(this_fig, '-property', 'FontSize'), 'FontSize', fontSize);
    drawnow;

    hold off;
  end
end

if saveToFile
  saveas(figure(1), "02_good_lambda_intervals_S0_ABC.png");
  saveas(figure(2), "02_good_lambda_intervals_S0_H.png");
  saveas(figure(3), "02_good_lambda_intervals_Cp_ABC.png");
  saveas(figure(4), "02_good_lambda_intervals_Cp_H.png");
end
