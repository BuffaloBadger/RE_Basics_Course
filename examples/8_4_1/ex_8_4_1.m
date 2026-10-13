function ex_8_4_1()
% Calculations for Example 8.4.1 in "Reaction Engineering Basics, The Book"
    % global constants available to all functions
    % given
    CA_0 = 15.0; % mol/L
    T_0 = 20 + 273.15; % K
    CB_in = 15.0; % mol/L
    T_in = T_0;
    V_0 = 5.0; % L
    Vdot_in = 0.25/60; % L/s
    V_B = 5.0; % L
    V_ex = 2.5; % L
    Tex_in = 20 + 273.15; % K
    Vdot_ex = 2.5/60; % L/s
    A = 2150; % cm^2
    U = 73/60/60; % cal/cm^2/s/K
    Tex_0 = 20 + 273.15; % K
    rho = 1000.0; % g/L
    Cp = 1.0; % cal/g/K
    dH_1 = -13700.0; % cal/mol
    k0_1 = 8.11E12; % L/mol/s
    E_1 = 17700.0; % cal/mol
    P = 1.0; % atm
    % known
    Re = 1.987; % cal/mol/K
    Rw = 0.08206; % L-atm/mol/K
    % calculated
    nDotB_in = Vdot_in*CB_in;
    nA_0 = CA_0*V_0;
    mDot_ex = Vdot_ex*rho;
    V_f = V_0 + V_B;

    % SBSTR reactor function
    function [t, nA, nB, nS, T, Tex, V] = sbstr_model_variables()
        % set the initial values
        ind_0 = 0.0;
        dep_0 = [nA_0, 0.0, 0.0, T_0, Tex_0, V_0];

        % define the stopping criterion
        stopVar = 6;
        stopVal = V_f;
        
        % solve the IVODEs
        odes_are_stiff = false;
        [t, dep, success, message] = solve_ivodes(ind_0, dep_0, stopVar...
            , stopVal, @sbstr_derivatives, odes_are_stiff);

        % check for solver issues
        if ~success
            disp('')
            disp(["SBSTR function issue: ", message])
            disp(' Press any key to continue or Ctrl-C to exit')
            pause
        end

        % extract the dependent variables
        nA = dep(:,1);
        nB = dep(:,2);
        nS = dep(:,3);
        T = dep(:,4);
        Tex = dep(:,5);
        V = dep(:,6);
    end

    % SBSTR derivatives function
    function ddt = sbstr_derivatives(~, dep)
        % extract the dependent variables
        nA = dep(1);
        nB = dep(2);
        T = dep(4);
        Tex = dep(5);
        V = dep(6);

        % calculate the additional unknowns
        k_1 = k0_1*exp(-E_1/Re/T);
        CA = nA/V;
        CB = nB/V;
        r = k_1*CA*CB;
        Qdot = U*A*(Tex - T);

        % evaluate the derivatives
        dnAdt = -V*r;
        dnBdt = nDotB_in -V*r;
        dnSdt = V*r;
        dTdt = (Qdot - rho*Vdot_in*Cp*(T - T_in) - V*r*dH_1 ...
            + P*Vdot_in*Re/Rw) /(rho*V*Cp);
        dTexdt = (-Qdot - mDot_ex*Cp*(Tex - Tex_in))/(rho*V_ex*Cp);
        dVdt = Vdot_in;

        % return the derivatives as a column vector
        ddt = [dnAdt; dnBdt; dnSdt; dTdt; dTexdt; dVdt];
    end

    % items_of_interest function
    function items_of_interest()
        % solve the reactor design equations
        [t, nA, ~, ~, T, ~, V] = sbstr_model_variables();

        % calculate the items of interest for plotting and convert units
        CA = nA./V;
        T_C = T - 273.15;
        t_min = t/60;

        % tabulate, show and save the final concentration and temperature
        item = ["Final Concentration of A"; "Final Temperature"];
        value = [CA(end); T_C(end)];
        units = ["M"; "°C"];
        resultsTable = table(item, value, units);
        disp('')
        disp(resultsTable)
        disp('')
        writetable(resultsTable, 'ex_8_4_1_results.csv')
        
        % display and save the graphs
        figure
        plot(t_min, CA, 'LineWidth', 2)
        xlabel("Time (min)", "FontSize", 14)
        ylabel("Concentration of A (mol/L)", "FontSize", 14)
        set(gca, 'FontSize', 14)
        saveas(gcf, 'ex_8_4_1_CA_vs_t.png')

        figure
        plot(t_min, T_C, "LineWidth", 2)
        xlabel("Time (min)", "FontSize", 14)
        ylabel("Temperature (°C)", "FontSize", 14)
        set(gca, 'FontSize', 14)
        saveas(gcf, 'ex_8_4_1_T_vs_t.png')
    end

    % perform the calculations
    items_of_interest()
end
        