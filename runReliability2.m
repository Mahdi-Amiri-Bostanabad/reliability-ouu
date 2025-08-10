function [betaF , pfF , pfS] = runReliability2(x_opt , Ndir , Nis)
% runReliability  –  FORM / SORM  + optional DS  + optional Importance Sampling
%
%  [betaF ,pfF ,pfS]             = runReliability(x_opt)
%  [...] = runReliability(x_opt , Ndir)               % + DS
%  [...] = runReliability(x_opt , Ndir , Nis)         % + DS + IS
%
%  Cd0  ~ log-normal ;   sigma ~ log-normal   (positivity enforced)

if nargin<2,  Ndir = 0; end          % default: skip DS
if nargin<3,  Nis  = 0; end          % default: skip IS

clip = @(z) max(min(z,4),-4);        % identical truncation everywhere

%% 0.  Constants & baseline
C            = load('baseline_constants.mat','baseline').baseline;
[r0 , V0]    = deal(x_opt(1), x_opt(2));

%% 1.  Random-variable parameters  (both log-normal)
[R,~]   = buildRandomInputs(false);
mu_ln   = log(R(1).mu) - 0.5*log(1+R(1).cov^2);
sig_ln  = sqrt(log(1+R(1).cov^2));               % Cd0

mu_sig  = log(R(2).mu) - 0.5*log(1+R(2).cov^2);
sig_sig = sqrt(log(1+R(2).cov^2));               % sigma

%% 2.  FORM (HL-RF)
betaF = zeros(3,1);  pfF = zeros(3,1);
uStar = zeros(2,3);  gradU=zeros(2,3);

for j = 1:3
    [betaF(j),pfF(j),uStar(:,j),gradU(:,j)] = ...
        formHLRF(@(u) gfun(u,j), 1e-4, 50);
end

%% 3.  SORM (Breitung)
pfS = zeros(3,1);
for j = 1:3
    kappa  = curvature(uStar(:,j),@(u) gfun(u,j),gradU(:,j));
    pfS(j) = pfF(j) ./ sqrt(1 + betaF(j)*kappa);
end

%% 4.  Directional Sampling  (optional)
if Ndir > 0
    [pfDS , cvDS] = directionalSampling(Ndir);
    fprintf('\nDS  (N=%d)  pf ≈ [%5.2e  %5.2e  %5.2e]   C.o.V ≈ [%4.2f  %4.2f  %4.2f]\n',...
            Ndir , pfDS , cvDS );
end

%% 5.  Importance Sampling around the FORM MPP  (optional)
if Nis > 0
    [pfIS , cvIS] = importanceSampling(Nis);
    fprintf('IS  (N=%d)  pf ≈ [%5.2e  %5.2e  %5.2e]   C.o.V ≈ [%4.2f  %4.2f  %4.2f]\n',...
            Nis , pfIS , cvIS );
end


%% 6.  Pretty print
nm = {'Disk load','Blade load','Energy'};
fprintf('\n--- FORM / SORM @ r = %.3f m ,  V = %.1f m/s ---\n',r0,V0);
fprintf('%-12s  %-6s  %-10s  %-10s\n','Limit','β','pf_FORM','pf_SORM');
for j = 1:3
    fprintf('%-12s  %6.2f  %10.2e  %10.2e\n',...
            nm{j},betaF(j),pfF(j),pfS(j));
end

% ========================= nested helpers ===============================

   function g = limitStates(Cd0,sigma)
       mu   = V0/(C.Omega0*r0);
       Pc   = (sigma*Cd0/8)*(1+4.65*mu^2)*1.225*C.A*C.Omega0^3*r0^3;
       Ereq = C.P_hover*4*C.thover + 4*Pc*(C.R/V0);
       
       DL   = C.T / C.A;
       CT   = C.T /(1.225*C.A*(C.Omega0*r0)^2);

       g(1) = C.DL_max - DL;
       g(2) = C.BL_max - CT/sigma;
       g(3) = C.E_use  - Ereq;
   end

   function [g,grad_u] = gfun(u,idx)
       Cd0   = exp(mu_ln  + sig_ln  * clip(u(1)));
       sigma = exp(mu_sig + sig_sig * clip(u(2)));

       gvec = limitStates(Cd0,sigma);   g = gvec(idx);

       rel = 1e-2;  x0=[Cd0 sigma];  dgx=zeros(1,2);
       for k = 1:2
           d  = rel*x0(k);
           gp = limitStates(x0(1)+(k==1)*d , x0(2)+(k==2)*d);
           gm = limitStates(x0(1)-(k==1)*d , x0(2)-(k==2)*d);
           dgx(k) = (gp(idx)-gm(idx))/(2*d);
       end
       grad_u = (dgx .* [sig_ln*Cd0 , sig_sig*sigma]).';
   end

   function [beta,pf,u_fin,grad_fin] = formHLRF(fun,tol,Nmax)
       u=[0;0];
       for k=1:Nmax
           [g,grad]=fun(u); ng=norm(grad);
           if ng<1e-5
               beta = sign(g)*inf; pf=(g<=0); u_fin=u; grad_fin=grad; return
           end
           alpha=grad/ng; beta=-g/ng; u_new=beta*alpha;
           if norm(u_new-u)<tol, u=u_new; break, end, u=u_new;
       end
       beta=norm(u); pf=normcdf(-beta); u_fin=u; grad_fin=grad;
   end

   function k=curvature(u0,fun,grad0)
       h=1e-3; H=zeros(2);
       for i=1:2
           e=zeros(2,1); e(i)=1;
           [~,gp]=fun(u0+h*e); [~,gm]=fun(u0-h*e);
           H(:,i)=(gp-gm)/(2*h);
       end
       t=[-grad0(2); grad0(1)]; t=t/norm(t);
       k=(t'*H*t)/max(norm(grad0),1e-8);
   end

  function [pf , cv] = directionalSampling(N)
% Directional Sampling with full trace saved to ds_debug.mat

    maxGrow  = 1e3;         % β cap if we never reach failure
    tolBisec = 1e-3;
    maxIter  = 40;

    fail   = zeros(1,3);
    b_enter  =  inf( N ,1);          %#ok<*NASGU>
    g_enter  =  nan( N ,3);
    g_design =  nan( N ,3);
    status   =  zeros( N ,1,'int8');

    for ii = 1:N
        d = randn(2,1);  d = d/norm(d);

        % -- grow bracket until any g≤0 ---------------------------------
        bL = 0.07;  bH = 0.2;
        while true
            uH = bH*d;
            CD0_e=exp(mu_ln  + sig_ln  * clip(uH(1)));
            sigma_e=exp(mu_sig + sig_sig * clip(uH(2)));

            gH = limitStates( ...
                    exp(mu_ln  + sig_ln  * clip(uH(1))) , ...
                    exp(mu_sig + sig_sig * clip(uH(2))) );
            if any(gH<=0) || bH>maxGrow, break, end
            bH = bH*2;
        end

        if bH>maxGrow           % never crossed → safe in this direction
            status(ii) = 0;
            continue
        end
        status(ii) = 1;
        b_enter(ii)= bH;
        g_enter(ii,:)= gH;

        % -- bisection ---------------------------------------------------
        for iter = 1:maxIter
            bM = 0.5*(bL+bH);   uM = bM*d;
            gM = limitStates( ...
                    exp(mu_ln  + sig_ln  * clip(uM(1))) , ...
                    exp(mu_sig + sig_sig * clip(uM(2))) );
            if all(gM>0), bL=bM; else, bH=bM; end
            if bH-bL<tolBisec, break, end
        end

        uF = bH*d;
                Cd0_f=exp(mu_ln  + sig_ln  * uF(1));
                sigma_f=exp(mu_sig + sig_sig * uF(2));
         % fprintf('DS  Cd0=%g  sig=%g  E_use=%g\n', Cd0_f, sigma_f, C.E_use);

        gF = limitStates(Cd0_f,sigma_f);
        g_design(ii,:) = gF;

        fail = fail + (gF<=0);
    end

    pf = fail / N;
    var = max( pf .* (1-pf) / N , 0 );
    cv  = sqrt(var) ./ max(pf,eps);

    % save ds_debug.mat  b_enter g_enter g_design status
end

   % ---------- Importance Sampling -------------------------------------
   function [pfIS,cvIS] = importanceSampling(N)
       pfIS = zeros(1,3);   varIS = zeros(1,3);

       % precompute MPPs and Cholesky of the shifted covariance (identity)
       mpp = uStar;                      % each column is u* for g_j
       eye2 = eye(2);

       for j = 1:3
           u0 = mpp(:,j);                % centre IS at FORM design point
           w_sum = 0;  w2_sum = 0;  nf = 0;

           for n = 1:N
               u_shift = u0 + randn(2,1);          % sample from N(u0,I)
               % importance weight  w = φ(u)/φ_shifted(u)
               % φ_shifted  = φ(u-u0)   because Σ = I
               w = exp( 0.5*(u_shift-u0).'*(u_shift-u0) ...
                         -0.5*u_shift.'*u_shift );

               Cd0   = exp(mu_ln  + sig_ln  * clip(u_shift(1)));
               sigma = exp(mu_sig + sig_sig * clip(u_shift(2)));
               gj    = limitStates(Cd0,sigma); gj = gj(j);

               I_f   = (gj <= 0);       % indicator of failure
               w_sum = w_sum + w*I_f;
               w2_sum= w2_sum+ w^2*I_f;
           end
           pfIS(j) = w_sum / N;
           varIS(j)= max( w2_sum/N - pfIS(j)^2 , 0 );
       end
       cvIS = sqrt(varIS) ./ max(pfIS,eps);
   end

dd=limitStates(0.0237,4.1222e-04);
end





