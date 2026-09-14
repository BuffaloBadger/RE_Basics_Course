function example_6_6_6()
% Calculations for Example 6.6.6 from REB, The Book
    % global constants available to all functions
    % given
    V = 500.0; % cm^3
    Vdot_in = 1.0; % cm^3 /s
    CA_in = 0.015; % mol /cm^3
    CB_in = 0.015; % mol /cm^3
    T_in = 50 + 273.15; % K
    Tss = [264.75, 50.06, 137.96] + 273.15; % K
    fAss = [0.9749, 0.03, 0.3993];
    dT = [20, 20, 1]; % K
    % constants
    Cp = 0.35*4.184; % J /g /K
    rho = 0.93; % g /cm^3
    dH = -20000.0; % J /mol
    k0 = 3.24E12; % cm^3 /mol /s
    E = 105000.0; % J /mol
    % known
    R = 8.314; % J /mol /K
    % calculated
    nA_in = CA_in*Vdot_in;
    nB_in = CB_in*Vdot_in;
    Vdot = Vdot_in;

    % CSTR model function
    function [t, nA, nB, nY, nZ, T] = cstr_model_variables(dep0, tf)
        % initial value
        ind0 = 0.0;

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
        nY = dep(:, 3);
        nZ = dep(:, 4);
        T = dep(:, 5);
    end

    % CSTR derivatives function
    function ddt = cstr_derivatives(~, dep)
        % unpack the dependent variables
        nA = dep(1);
        nB = dep(2);
        nY = dep(3);
        nZ = dep(4);
        T = dep(5);

        % calculate the additional unknowns
        k = k0*exp(-E/(R*T));
        CA = nA/Vdot;
        CB = nB/Vdot;
        r = k*CA*CB;

        % calculate the derivatives
        dnAdt = Vdot/V*(nA_in - nA - V*r);
        dnBdt = Vdot/V*(nB_in - nB - V*r);
        dnYdt = Vdot/V*( -nY + V*r);
        dnZdt = Vdot/V*( -nZ + V*r);
        dTdt = -(Cp*Vdot*rho*(T-T_in) + V*r*dH)/(V*Cp*rho);

        % return the derivatives as a column vector
        ddt = [dnAdt; dnBdt; dnYdt; dnZdt; dTdt];
    end

    % deliverables function
    function deliverables()
        % loop through the steady states
        nSS = length(Tss);
        for iSS = 1:nSS
            % set the initial values for a positive perturbation
            nA0 = nA_in*(1-fAss(iSS));
            nB0 = nB_in - fAss(iSS)*nA_in;
            nY0 = fAss(iSS)*nA_in;
            nZ0 = fAss(iSS)*nA_in;
            T0 = Tss(iSS) + dT(iSS);
            dep0 = [nA0; nB0; nY0; nZ0; T0];

            % set the final value
            tf = 2500.0; % s

            % solve the design equations
            [tPos, ~, ~, ~, ~, TPos] = cstr_model_variables(dep0, tf);

            % repeat for a negative perturbation
            T0 = Tss(iSS) - dT(iSS);
            dep0 = [nA0; nB0; nY0; nZ0; T0];
            [tNeg, ~, ~, ~, ~, TNeg] = cstr_model_variables(dep0, tf);

            % convert units
            tPos = tPos/60.0; % min
            tNeg = tNeg/60.0; % min
            TPos = TPos - 273.15; % C
            TNeg = TNeg - 273.15; % C

            % plot, show, and save the results
            figure;
            plot(tPos, TPos, 'r', tNeg, TNeg, 'b', 'LineWidth', 2);
            xlabel('Time (min)');
            ylabel('Temperature (°C)');
            if iSS == 3
                yline(Tss(1) - 273.15, 'k-', 'LineWidth', 2);
                yline(Tss(2) - 273.15, 'k-', 'LineWidth', 2);
                yline(Tss(3) - 273.15, 'k--', 'LineWidth', 2);
            else
                yline(Tss(iSS) - 273.15, 'k-', 'LineWidth', 2);
            end
            saveas(gcf, ['perturbation_response_ss_' num2str(iSS) '.png']);
        end
    end

    % perform the calculations
    deliverables();
end