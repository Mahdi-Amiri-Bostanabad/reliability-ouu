%% Deterministic baseline sizing & optimisation 
% Author: Mahdi Amiribostanabad & Mohammad Bagheri
% Date: Jul‑2025 
% -------------------------------------------------------------------------
clear; clc;

%% CONSTANTS & INPUT DATA --------------------------------------------------
rho  = 1.225;         % [kg/m^3] air density
g    = 9.81;          % [m/s^2]  gravity

% Mission / payload
m_payload = 2.0;      % [kg]
R         = 50e3;     % [m]   one‑way range
thover    = 60;       % [s]   hover per leg (quad has 4 running legs)

% Aerodynamic / rotor parameters (deterministic)
sigma  = 0.13;        % [–]    rotor solidity 
Cd0    = 0.012;       % [–]    zero‑lift drag coef
V_tip  = 200;         % [m/s]  fixed tip speed 

% Battery data
rho_b_Wh  = 158;      % [Wh/kg] specific energy
eta_batt  = 0.85;     % charge–discharge efficiency
margin    = 0.30;     % keep 30 % SOC unused
rho_b     = rho_b_Wh*3600;  % [J/kg]

% Structural limits
DL_max = 250;         % [N/m^2]
BL_max = 0.14;        % [–]

%% DESIGN VARIABLES -------------------------------------------------------
% x = [r; V_inf] – rotor radius [m] & cruise speed [m/s]
x0 = [0.2; 25];
lb = [0.15 ; 10];        % r ≥ 0.15 m,  V ≥ 10 m/s
ub = [0.4 ; 80];         % r ≤ 0.4 m,  V ≤ 80 m/s


%% OPTIMISER --------------------------------------------------------------
opts = optimoptions('fmincon','Algorithm','sqp', ...
                    'OptimalityTolerance',1e-8, ...
                    'ConstraintTolerance',1e-8, ...
                    'StepTolerance',1e-8, ...
                    'Display','iter', ...
                    'MaxFunctionEvaluations',2e3);
problem = createOptimProblem('fmincon', 'objective', @objFun, ...
                             'nonlcon', @nlcon, 'x0', x0, ...
                             'lb', lb, 'ub', ub, 'options', opts);
[x_opt, m_opt] = fmincon(problem);

fprintf('\nOptimal radius   : %.3f  m', x_opt(1));
fprintf('\nOptimal speed    : %.1f  m/s', x_opt(2));
fprintf('\nTake‑off weight  : %.2f N (%.2f kg)\n', m_opt*g, m_opt);

ggg=limitStates(x_opt)

%% Output-------------------------------------------------------------------
[~, st_opt] = sizingModel(x_opt);
baseline = struct( ...
   'A',        pi*x_opt(1)^2 , ...
   'T',        m_opt*g/4     , ...
   'Omega0',   st_opt.Omega  , ...
   'P_hover',  (1/0.75)*(m_opt*g/4)^1.5/sqrt(2*rho*pi*x_opt(1)^2) , ...
   'E_use',    st_opt.E_usable , ...
   'DL_max',   DL_max , ...
   'BL_max',   BL_max , ...
   'thover',   thover , ...
   'R',        R );
save baseline_constants.mat baseline

%% OBJECTIVE FUNCTION -----------------------------------------------------
function m = objFun(x)
    m = sizingModel(x);  % returns total **mass** [kg]
end

%% NON-LINEAR CONSTRAINTS -------------------------------------------------
function [c,ceq] = nlcon(x)
    g   = limitStates(x);
    c   = -(g);      
    ceq = [];
end

%% SIZING MODEL -----------------------------------------------------------
function [m_total, st] = sizingModel(x)
    % Pull shared variables from caller
    g          = evalin('base','g');
    rho        = evalin('base','rho');
    sigma      = evalin('base','sigma');
    Cd0        = evalin('base','Cd0');
    rho_b      = evalin('base','rho_b');
    eta_batt   = evalin('base','eta_batt');
    margin     = evalin('base','margin');
    V_tip      = evalin('base','V_tip');
    thover     = evalin('base','thover');
    R          = evalin('base','R');
    m_payload  = evalin('base','m_payload');
    DL_max     = evalin('base','DL_max');
    BL_max     = evalin('base','BL_max');

    % Design vars
    r    = x(1);                     % [m]
    Vinf = x(2);                     % [m/s]

    % Initial mass guess [kg]
    m_total = m_payload + 2.0;

    for k = 1:50
        W_total = m_total * g;       % [N]
        A = pi*r^2;                  % disk area per rotor
        T = W_total/4;              % thrust per rotor [N]

        % Hover power (Eq.6)
        P_hover = (1/0.75)*T^1.5/sqrt(2*rho*A);

        % ---- fast numeric solve for Ω that satisfies  T·Vinf = P_cruise  ---------
        % K = (sigma*Cd0/8)*rho*A*r^3;          % common factor in P_cruise
        % a3 = K;                               % Ω³ coefficient
        % a1 = K*4.65*Vinf^2/r^2;               % Ω¹ coefficient
        % a0 = -T*Vinf;                         % constant term  (move to LHS)
        % 
        % Omega_candidates = roots([a3 0 a1 a0]);           % cubic: a3Ω³ + a1Ω + a0 = 0
        % Omega = Omega_candidates( ...
        %           imag(Omega_candidates)==0 & real(Omega_candidates)>0);  % pick real +
        % Omega = Omega(1);                    % there’s only one positive real root
        
        % Omega = sqrt( T / ( 0.3 * (BL_max * sigma) * rho * A * r^2 ) );

        Omega=1100;

        mu       = Vinf/(Omega*r);           % advance ratio with solved Ω
        P_cruise = (sigma*Cd0/8)*(1 + 4.65*mu^2)*rho*A*Omega^3*r^3;  % now TV = Pcruise
        % ---------------------------------------------------------;

        % Mission energy requirement [J]
        E_req = P_hover*4*thover + 4*P_cruise*(R/Vinf);

        % Battery mass sized to leave 30 % SOC unused
        extra_m = 0.10;
        m_batt = E_req / ((1-margin-extra_m)*eta_batt*rho_b); % [kg]
        E_usable = (1-margin)*eta_batt*rho_b*m_batt;

        % Propulsion system empirical weights 
        P_installed = 4*max(P_hover,P_cruise);  % total installed power [W]
        m_motors = 2.506e-4 * P_installed;      % [kg]
        m_ESC    = 3.594e-4 * P_installed;      % [kg]
        m_rotors = 4*(0.7484*r^2 - 0.0403*r);   % [kg]

        % Frame mass scales with total mass (empirical)
        m_empty_guess = m_motors + m_ESC + m_rotors;
        m_total_tmp   = m_payload + m_batt + m_empty_guess;
        m_frame       = 0.5 + 0.2*m_total_tmp;   % [kg]
        m_empty       = m_empty_guess + m_frame; % [kg]

        m_new = m_payload + m_batt + m_empty;    % [kg]
        if abs(m_new - m_total) < 1e-4
            m_total = m_new; break; end
        m_total = m_new;
     end

    % Constraint‑related quantities
    W_total = m_total * g;           % [N]
    T = W_total/4;
    A = pi*r^2;
    DL = T/A;
    CT = T/(rho*A*(Omega*r)^2);

    % Package
    st = struct('DL',DL,'DL_max',DL_max, ...
                'BL',CT/sigma,'BL_max',BL_max, ...
                'E_req',E_req,'E_usable',E_usable, ...
                'diskArea',A,'Omega',Omega);
end




%% LIMIT-STATE HELPER ------------------------------------------------------
function [g, st] = limitStates(x)
    % g(1) : disk-loading margin      DL_max − DL
    % g(2) : blade-loading margin     BL_max − CT/σ
    % g(3) : energy-reserve margin    E_usable − E_req
    [~, st] = sizingModel(x);            % physics & book-keeping
    g = [ st.DL_max - st.DL ; ...
          st.BL_max - st.BL ; ...
          st.E_usable - st.E_req ];
end









