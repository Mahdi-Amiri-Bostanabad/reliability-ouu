function [R, sampleInputs] = buildRandomInputs(makePlots)
% STEP3_UNCERTAINTY  Build minimal stochastic model (Cd0, sigma)
%   [R, sampleInputs] = step3_uncertainty();          % plots PDFs
%   [...] = step3_uncertainty(false);                 % silent, no plots
%
% Outputs pushed to caller workspace too, so scripts can see them.
if nargin < 1,  makePlots = true;  end

% ---------- 1. PDFs -------------------------------------------------------
mu_Cd0  = 0.012;  cov_Cd0 = 0.20;
sigma_ln = sqrt(log(1+cov_Cd0^2));
mu_ln    = log(mu_Cd0) - 0.5*sigma_ln^2;
pdf_Cd0  = makedist('Lognormal','mu',mu_ln,'sigma',sigma_ln);

mu_sig = 0.13;   cov_sig = 0.12;
pdf_sig = makedist('Normal','mu',mu_sig,'sigma',cov_sig*mu_sig);

R = struct( ...
    'name',  {'Cd0','sigma'}, ...
    'units', {'–','–'}, ...
    'dist',  {pdf_Cd0 , pdf_sig}, ...
    'mu',    num2cell([mu_Cd0  , mu_sig ]), ...
    'cov',   num2cell([cov_Cd0 , cov_sig]) );


% ---------- 2. Sampler handle -------------------------------------------
sampleInputs = @(n) sampler(n,R);

% ---------- 3. Push to caller workspace for easy use ---------------------
assignin('base','R',R);
assignin('base','sampleInputs',sampleInputs);

% ---------- 4. Optional sanity-check plots ------------------------------
if makePlots
    figure('Name','Step-3 PDFs'), clf
    subplot(1,2,1)
        x = linspace(0,0.025,400);
        plot(x, pdf(R(1).dist,x)), grid on
        title('Lognormal PDF  C_d_0'), xlabel('C_d_0')
    subplot(1,2,2)
        x2 = linspace(0.08,0.18,400);
        plot(x2, pdf(R(2).dist,x2)), grid on
        title('Normal PDF  \sigma'), xlabel('\sigma')
end
end  %------------- end main function -------------------------------------

% ============ nested sampler ============================================
function X = sampler(n,R)
    Z = randn(n,2);                      % independent normals
    X(:,1) = icdf(R(1).dist, normcdf(Z(:,1)));   % Cd0
    X(:,2) = icdf(R(2).dist, normcdf(Z(:,2)));   % sigma
end
