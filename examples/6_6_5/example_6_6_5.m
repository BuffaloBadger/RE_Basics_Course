function example_6_6_5()
% Calculations for Example 6.6.5 from REB, The Book
    % global constants available to all functions
    % given
    V = 10; % gal
    Ve = 1.25; % gal
    CB0 = 3; % mol/gal
    T0 = 20 + 273.15; % K
    VFRin = 0.5; % gal/min
    CAin = 7.5; % mol/gal
    CBin = 3; % mol/gal
    Tin = 50 + 273.15; % K
    Texin = 20 + 273.15; % K
    Tex0 = 20 + 273.15; % K
    mex_in = 250; % g/min
    U = 190; % cal/ft^2/min/K
    A = 4; % ft^2
    rhoex = 1; % g/cc
    Cpex = 1; % cal/g/K
    Cp = 1600; % cal/gal/K
    tf = 30; % min
    k01 = 4.3E8; % gal/mol/min
    E1 = 14200; % cal/mol
    dH1 = -11000; % cal/mol
    k02 = 2.7E8; % gal/mol/min
    E2 = 16100; % cal/mol
    dH2 = -11600; % cal/mol
    k03 = 3.9E8; % gal/mol/min
    E3 = 14800; % cal/mol
    dH3 = -12100; % cal/mol
    % known
    R = 1.987; % cal/mol/K
    % calculated
    VFR = VFRin;
    nB0 = VFR*CB0;
    nAin = CAin*VFRin;
    nBin = CBin*VFRin;

    % CSTR model function
    function [t, nA, nB, nW, nX, nY, nZ, T, Tex] = cstr_model_variables()
        % initial values
        ind0 = 0.0;
        dep0 = [0.0, nB0, 0.0, 0.0, 0.0, 0.0, T0, Tex0];

        % stopping criterion
        stopVar = 0;
        stopVal = tf;

        % solve the IVODEs
        odesAreStiff = true;
        [ind, dep, flag, message] = solve_ivodes(ind0, dep0, stopVar...
            , stopVal, @cstr_derivatives, odesAreStiff);
        
        % check for solver issues
        if flag <= 0
            disp(' ')
            disp(['     CSTR model function issue: ' message    ])
            disp(' ')
        end

        % return the individual CSTR model variables
        t = ind;
        nA = dep(:, 1);
        nB = dep(:, 2);
        nW = dep(:, 3);
        nX = dep(:, 4);
        nY = dep(:, 5);
        nZ = dep(:, 6);
        T = dep(:, 7);
        Tex = dep(:, 8);
    end

    % CSTR derivatives function
    function ddt = cstr_derivatives(~, dep)
        % unpack the dependent variables
        nA = dep(1);
        nB = dep(2);
        nW = dep(3);
        nX = dep(4);
        nY = dep(5);
        nZ = dep(6);
        T = dep(7);
        Tex = dep(8);

        % calculate the additional unknowns
        k1 = k01*exp(-E1/(R*T));
        k2 = k02*exp(-E2/(R*T));
        k3 = k03*exp(-E3/(R*T));
        CA = nA/VFR;
        CB = nB/VFR;
        r1 = k1*CA*CB;
        r2 = k2*CA*CB;
        r3 = k3*CA*CB;
        Q = U*A*(Tex - T);

        % calculate the derivatives
        dndt_A = VFR/V*(nAin - nA - V*(r1 + r2 + r3));
        dndt_B = VFR/V*(nBin - nB - V*(r1 + r2 + r3));
        dndt_W = VFR/V*(-nW + r1*V);
        dndt_X = VFR/V*(-nX + r2*V);
        dndt_Y = VFR/V*(-nY + r3*V);
        dndt_Z = VFR/V*(-nZ + V*(r1 + r2 + r3));
        dTdt = (Q - VFR*Cp*(T-Tin) - V*(r1*dH1 + r2*dH2 + r3*dH3))/V/Cp;
        dTexdt = (-Q-mex_in*Cpex*(Tex-Texin))/(rhoex*Ve*Cpex);

        % return the derivatives as a column vector
        ddt = [dndt_A; dndt_B; dndt_W; dndt_X; dndt_Y; dndt_Z; dTdt; dTexdt];
    end

    % deliverables function
    function deliverables()
        % solve the design equations
        [t, ~, nB, ~, ~, ~, ~, T, Tex] = cstr_model_variables();

        % calculate the quantities of interest
        CB = nB/VFR;
        T = T - 273.15;
        Tex = Tex - 273.15;

        % display and save the graphs
        figure;
        plot(t, T, 'k', 'LineWidth', 2);
        xlabel('Time (min)');
        ylabel('Reacting Fluid Temperature (C)');
        set(gca, 'FontSize', 14);
        saveas(gcf, 'T_vs_t.png');

        figure;
        plot(t, Tex, 'k', 'LineWidth', 2);
        xlabel('Time (min)');
        ylabel('Exchange Fluid Temperature (C)');
        set(gca, 'FontSize', 14);
        saveas(gcf, 'Tex_vs_t.png');

        figure;
        plot(t, CB, 'k', 'LineWidth', 2);
        xlabel('Time (min)');
        ylabel('Concentration of B (mol/gal)');
        set(gca, 'FontSize', 14);
        saveas(gcf, 'CB_vs_t.png');
    end

    % perform the calculations
    deliverables();
end
