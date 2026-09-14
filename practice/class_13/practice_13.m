function practice_13()
% Calculations for the Class 13 Practice Assignment from REB, The Course
    % global constants available to all functions
    % given
    V = 17.7; % ft^3
    CA_in = 0.125; % lbmol/ft^3
    Vdot_in = 0.15; % ft^3/min
    T_in = 300 + 458.67; % °R
    Cp = 73.1; % BTU/ft^3/°F
    k01 = 3.04e8; % ft^3/lbmol/min
    E1 = 29100; % BTU/lbmol
    dH1 = 32800; % BTU/lbmol
    rho_ex = 81.2; % lb/ft^3
    Cp_ex = 1.06; % BTU/lb/°F
    mDotex = 5.0; % lb/min
    Tex_in = 410 + 458.67; % °R
    Tex_0 = 355 + + 458.67; % °R
    Vex = 7.0; % ft^3
    A = 21.5; % ft^2
    U = 1.71; % BTU/ft^2/min/°F
    fA_ss = 0.853;
    T_0 = 332 + 458.67; % °R
    tf = 10*60; % min
    % known
    R = 1.986; % BTU/lbmol/°R
    % calculated constants
    Vdot = Vdot_in;
    nA_in = CA_in*Vdot_in;
    nA_0 = nA_in*(1-fA_ss);
    nZ_0 = nA_in*fA_ss;

    % CSTR model function
    function [t, nA, nZ, T, Tex] = cstr_model_variables()
        % set the initial values
        ind_0 = 0;
        dep_0 = [nA_0, nZ_0, T_0, Tex_0];

        % define the stopping criterion
        StopVar = 0;
        StopVal = tf;

        % solve the cstr design equations
        odes_are_stiff = false;
        %odes_are_stiff = true
        [t, dep, success, message] = solve_ivodes(ind_0, dep_0, StopVar...
            , StopVal, @cstr_derivatives, odes_are_stiff);
        
        % check for solver issues
        if ~success
            disp('')
            disp(["CSTR model function issue: ",message])
            disp('')
        end

        % extract the cstr model variables
        nA = dep(:,1);
        nZ = dep(:,2);
        T = dep(:,3);
        Tex = dep(:,4);
    end

    % cstr derivatives function
    function ddt = cstr_derivatives(~, dep)
        % extract the dependent variables
        nA = dep(1);
        nZ = dep(2);
        T = dep(3);
        Tex = dep(4);

        % calculate the additional unknowns
        k1 = k01*exp(-E1/(R*T));
        CA = nA/Vdot;
        r1 = k1*CA^2;
        Q = U*A*(Tex - T);

        % evaluate the derivatives
        dnAdt = Vdot/V*(nA_in - nA - r1*V);
        dnZdt = Vdot/V*(-nZ + r1*V);
        dTdt = 1/V/Cp*(-Vdot*Cp*(T - T_in) - r1*V*dH1 + Q);
        dTexdt = 1/rho_ex/Vex/Cp_ex*(-Q - mDotex*Cp_ex*(Tex - Tex_in));

        % return the derivatives as a vector
        ddt = [dnAdt; dnZdt; dTdt; dTexdt];
    end

    % deliverables function
    function deliverables()
        % get the cstr model variables
        [t, ~, nZ, T, ~] = cstr_model_variables();

        % calculate the concentration of Z
        CZ = nZ/Vdot;

        % generate, show, and save the graphs
        figure
        plot(t, CZ, 'LineWidth', 2)
        xlabel('Time (min)', 'FontSize', 14)
        ylabel('Concentration of Z (lbmol ft^-^3)', 'FontSize', 14)
        set(gca, 'FontSize', 14)
        saveas(gcf, 'CZ_vs_t.png')

        figure
        plot(t,T - 458.67, 'LineWidth', 2)
        xlabel('Time (min)', 'FontSize', 14)
        ylabel('Outlet Temperature (°F)', 'FontSize', 14)
        set(gca, 'FontSize', 14)
        saveas(gcf, 'T_vs_t.png')
    end

    %execution command
    deliverables()
end