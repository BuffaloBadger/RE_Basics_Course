function activity_12()
% Calculations for the Class 12 Learning Activity from REB, The Course
    % global constants available to all functions
    % given
    V = 500.0; % cm^3
    Vdot_in = 1.0; % cm^3 /s
    CA_in = 0.015; % mol /cm^3
    CB_in = 0.015; % mol /cm^3
    T_in = 50 + 273.15; % K
    Cp = 0.35*4.184; % J /g /K
    rho = 0.93; % g /cm^3
    dH = -20000.0; % J /mol
    k0 = 3.24E12; % cm^3 /mol /s
    E = 105000.0; % J /mol
    T0_values = [180, 190] + 273.15; % K
    % known
    R = 8.314; % J /mol /K
    % calculated
    Vdot = Vdot_in;
    nA_in = CA_in*Vdot_in;
    nB_in = CB_in*Vdot_in;

    % CSTR model function
    function [t, nA, nB, nY, nZ, T] = cstr_model_variables(T0, tf)
        % set the initial values
        ind_0 = 0.0;
        dep_0 = [0, 0, 0, 0, T0];

        % define the stopping criterion
        StopVar = 0;
        StopVal = tf;
        
        % solve the IVODEs
        odes_are_stiff = false;
        [t, dep, success, message] = solve_ivodes(ind_0, dep_0, StopVar...
            , StopVal, @cstr_derivatives, odes_are_stiff);

        % check that a solution was found
        if ~success
            disp('')
            disp(["CSTR model function issue: ", message])
            disp('')
        end

        % extract the individual cstr_model_variables
        nA = dep(:,1);
        nB = dep(:,2);
        nY = dep(:,3);
        nZ = dep(:,4);
        T = dep(:,5);
    end

    % CSTR derivatives function
    function ddt = cstr_derivatives(~, dep)
        % extract the dependent variables for this integration step
        nA = dep(1);
        nB = dep(2);
        nY = dep(3);
        nZ = dep(4);
        T = dep(5);

        % calculate the additional unknowns
        k = k0*exp(-E/R/T);
        CA = nA/Vdot;
        CB = nB/Vdot;
        r = k*CA*CB;
        
        % evaluate the derivatives
        dnAdt = Vdot/V*(nA_in - nA - V*r);
        dnBdt = Vdot/V*(nB_in - nB - V*r);
        dnYdt = Vdot/V*( -nY + V*r);
        dnZdt = Vdot/V*( -nZ + V*r);
        dTdt = -(Cp*Vdot_in*rho*(T-T_in) + V*r*dH)/(V*Cp*rho);

        % return the derivatives as a column vector
        ddt = [dnAdt; dnBdt; dnYdt; dnZdt; dTdt];
    end

    % deliverables function
    function deliverables()
        % set the final value of t
        tf = 2500; % s

        % solve the design equations for the first initial temperature
        [t1, ~, ~, ~, ~, T1] = cstr_model_variables(T0_values(1), tf);

        % solve the design equations for the second initial temperature
        [t2, ~, ~, ~, ~, T2] = cstr_model_variables(T0_values(2), tf);

        % convert to °C
        T1 = T1 - 273.15;
        T2 = T2 - 273.15;

        % generate, show, and save as a graph
        figure;
        plot(t1, T1, 'r', t2, T2, 'b', 'LineWidth', 2)
        set(gca, 'FontSize', 14);
        xlabel('Time (min)','FontSize', 14)
        ylabel('Reacting Fluid Temperature (°C)','FontSize', 14)
        legend({'180 °C', '190 °C'}, 'FontSize', 14)
        saveas(gcf,"T_profile.png")
    end

    % perform the calculations
    deliverables();
end