function ex_9_6_3()
% Calculations for Example 9.6.3 from REB, The Book
    % global constants available to all functions
    % given
    yA_in = 0.76;
    yI_in = 0.25;
    Vdot_in = 100; % cm^3 /s
    mDot_in = 0.44; % g /s
    P_in = 3.0; % atm
    T_in = 400 + 273.15; % K
    D = 2.5; % cm
    L = 800.; % cm
    T_ex = 375 + 273.15; % K
    U = 187E-1/3600; % J /s /cm^2 /K
    Dp = 0.25; % cm
    phi = 0.7;
    eps = 0.6;
    k0f = 9E17; % mol /cm^3 /s /atm
    Ef = 285E3; % J /mol
    k0r = 4.09E-4; % mol /cm^3 /s /atm^4
    Er = 85E3; % J /mol
    dH = 200E3; % J /mol
    CpA = 11.7*4.184; % J /mol /K
    CpY = 8.3*4.184; % J /mol /K
    CpZ = 4.2*4.184; % J /mol /K
    CpI = 5.8*4.184; % J /mol /K
    mu = 0.027E-2; % g /cm /s
    % known
    Re = 8.314; % J /mol /K
    Rw = 82.06; % cm^3 atm /mol /K
    P_conv = 9.872E-7; % atm cm^2 /dyne

    % PFR reactor function
    function [z, nA, nY, nZ, nI, T, P] = pfr_model_variables()
        % set the initial values
        ind_0 = 0.0;
        nA_in = yA_in*Vdot_in*P_in/Rw/T_in;
        nI_in = yI_in*Vdot_in*P_in/Rw/T_in;
        dep_0 = [nA_in, 0.0, 0.0, nI_in, T_in, P_in];

        % define the stopping criterion
        stopVar = 0;
        stopVal = L;
        
        % solve the IVODEs
        odes_are_stiff = false;
        [z, dep, success, message] = solve_ivodes(ind_0, dep_0, stopVar...
            , stopVal, @pfr_derivatives, odes_are_stiff);

        % check that a solution was found
        if ~success
            disp(' ')
            disp(['PFR model function issue: ', message])
            disp('Press any key to continue or CTRL-C to exit.')
            pause
        end

        % extract the dependent variables to be returned
        nA = dep(:,1);
        nY = dep(:,2);
        nZ = dep(:,3);
        nI = dep(:,4);
        T = dep(:,5);
        P = dep(:,6);
    end

    % PFR derivatives function
    function ddz = pfr_derivatives(~, dep)
        % extract the dependent variables
        nA = dep(1);
        nY = dep(2);
        nZ = dep(3);
        nI = dep(4);
        T = dep(5);
        P = dep(6);

        % calculate the additional unknowns
        kf = k0f*exp(-Ef/Re/T);
        kr = k0r*exp(-Er/Re/T);
        ntot = nA + nY + nZ + nI;
        PA = nA/ntot*P;
        PY = nY/ntot*P;
        PZ = nZ/ntot*P;
        r = kf*PA - kr*PY*PZ^3;
        Vdot = ntot*Rw*T/P;
        rho = mDot_in/Vdot;
        G = mDot_in/(pi*D^2/4);

        % evaluate the derivatives
        dnAdz = -pi*D^2/4*r;
        dnYdz = pi*D^2/4*r;
        dnZdz = 3*pi*D^2/4*r;
        dnIdz = 0.0;
        dTdz = (pi*D*U*(T_ex - T) - pi*D^2/4*r*dH)/(nA*CpA + nY*CpY + nZ*CpZ...
            + nI*CpI);
        dPdz = -P_conv*(1-eps)/eps^3*G^2/rho/phi/Dp*(150*(1-eps)*mu/phi/Dp/G...
            + 1.75);

        % return the derivatives as a column vector
        ddz = [dnAdz; dnYdz; dnZdz; dnIdz; dTdz; dPdz];
    end

    % items of interest function
    function items_of_interest()
        % solve the reactor design equations
        [~, nA, ~, ~, ~, T, P] = pfr_model_variables();

        % calculate the quantities of interest
        T_out = T(end) - 273.15;
        P_out = P(end);
        nA_in = yA_in*Vdot_in*P_in/Rw/T_in;
        fA = 100*(nA_in - nA(end))/nA_in;

        % tabulate, show and save the results
        item = ["Outlet T"; "Outlet P"; "Conversion"];
        value = [T_out; P_out; fA];
        units = ["°C"; "atm"; "%"];
        resultsTable = table(item, value, units);
        disp(' ')
        disp(resultsTable)
        disp(' ')
        writetable(resultsTable, 'results_ex_9_6-3.csv');
    end

    % perform the calculations
    items_of_interest();
end