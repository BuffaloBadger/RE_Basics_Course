function ex_8_4_3()
% Calculations for Example 8.4.3 from REB, The Book
    % global constants available to all functions
    % given
    k0_1 = 1.192E15; % /min
    E_1 = 97600; % J/mol
    dH_1 =-58615; % J/mol
    V_ex = 300; % cm^3
    Tex_in = 60 + 273.15; % K
    Vdot_ex = 250; % cm^3/min
    UA = 260*4.184; % J/min/K
    P = 1.0; % atm
    VW_0 = 67; % cm^3
    VZ_0 = 283; % cm^3
    Vacid = 0.3; % cm^3
    T_0 = 60 + 273.15; % K
    Tex_0 = 60 + 273.15; % K
    Cp = 2.68; % J/cm^3/K
    V_A = 350; % cm^3
    T_in = 21 + 273.15; % K
    T_f = 65 + 273.15; % K
    T_max = 95 + 273.15; %K
    rho_A = 1.082; % g/cm^3
    rho_W = 1.0; % g/cm^3
    rho_Z = 1.049; % g/cm^3
    M_A = 102; % g/mol
    M_W = 18; % g/mol
    M_Z = 60; % g/mol
    Cp_A = 168.2; % J/mol/K
    rho_ex = 1.0; % g/cm^3
    Cp_ex = 1.0*4.184; % cal/g/K
    % known
    Re = 1.987*4.184; % J/mol/K
    Rw = 82.06; % cm^3-atm/mol/K
    % calculated
    V_0 = VW_0 + VZ_0 + Vacid;
    m_ex = Vdot_ex*rho_ex;
    nW_0 = VW_0*rho_W/M_W;
    nZ_0 = VZ_0*rho_Z*M_Z;

    % global variable for the current inlet volumetric flow rate
    g_Vdot_in = nan;

    % SBSTR reactor function
    function [t, nA, nW, nZ, T, Tex, V] = sbstr_model_variables(Vdot_in)
        % set the current value of the inlet volumetric flow rate
        g_Vdot_in = Vdot_in;

        % set the initial values for the first stage
        ind_0 = 0.0;
        dep_0 = [0.0, nW_0, nZ_0, T_0, Tex_0, V_0];

        % set the stopping criterion for the first stage
        t_1 = V_A/Vdot_in;
        stopVar = 0;
        stopVal = t_1;

        % solve the design equations for the first stage
        odes_are_stiff = false;
        [t, dep, success, message] = solve_ivodes(ind_0, dep_0, stopVar...
            , stopVal, @sbstr_derivatives, odes_are_stiff);
        
        % check for solver issues
        if ~success
            disp('')
            disp(["SBSTR model first stage issue: ", message])
            disp('Press any key to continue or CTRL-C to exit')
            pause
        end
        
        % extract the dependent variables
        nA = dep(:,1);
        nW = dep(:,2);
        nZ = dep(:,3);
        T = dep(:,4);
        Tex = dep(:,5);
        V = dep(:,6);

        % repeat for the second stage
        g_Vdot_in = 0;
        ind_0 = t(end);
        dep_0 = [dep(end,1), dep(end,2), dep(end,3), dep(end,4), dep(end,5)...
            , dep(end,6)];
        stopVar = 4;
        stopVal = T_f;
        [t2, dep2, success, message] = solve_ivodes(ind_0, dep_0, stopVar...
            , stopVal, @sbstr_derivatives, odes_are_stiff);
        if ~success
            disp('')
            disp(["SBSTR model second stage issue: ", message])
            disp('Press any key to continue or CTRL-C to exit')
            pause
        end
        
        % combine the results for the two stages
        t = [t; t2];
        nA = [nA; dep2(:,1)];
        nW = [nW; dep2(:,2)];
        nZ = [nZ; dep2(:,3)];
        T = [T; dep2(:,4)];
        Tex = [Tex; dep2(:,5)];
        V = [V; dep2(:,6)];
    end

    % SBSTR derivatives function
    function ddt = sbstr_derivatives(~, dep)
        % extract the dependent variables
        nA = dep(1);
        T = dep(4);
        Tex = dep(5);
        V = dep(6);

        % calculate the additional unknowns
        k_1 = k0_1*exp(-E_1/Re/T);
        CA = nA/V;
        r_1 = k_1*CA;
        Qdot = UA*(Tex-T);
        nA_in = g_Vdot_in*rho_A/M_A;
        
        % evaluate the derivatives
        dnAdt = nA_in - V*r_1;
        dnWdt = - V*r_1;
        dnZdt = V*r_1;
        dTdt = (Qdot - nA_in*Cp_A*(T - T_in) - V*r_1*dH_1 ...
                + P*g_Vdot_in*Re/Rw)/(Cp*V);
        dTexdt = (-Qdot - m_ex*Cp_ex*(Tex-Tex_in))/(rho_ex*V_ex*Cp_ex);
        dVdt = g_Vdot_in;

        % return the derivatives as a column vector
        ddt = [dnAdt; dnWdt; dnZdt; dTdt; dTexdt; dVdt];
    end

    % coupled unknown residual function
    function epsilon = coupled_unknown_residual(guess)
        % solve the SBSTR design equations using the guess
        [~, ~, ~, ~, T, ~, ~] = sbstr_model_variables(guess);

        % evaluate and return the residual
        epsilon = max(T) - T_max;
    end

    % items of interest function
    function items_of_interest
        % define a guess for the volumetric feed rate
        VdotGuess = 35;

        % solve the implicit equation for the Vdot_opt
        [soln, success, message] = solve_ates(@coupled_unknown_residual...
            , VdotGuess);

        % check for solver issues
        if ~success
            disp('')
            disp("Coupled unknown solver issue: ")
            disp(message)
            disp('Press any key to continue or CTRL-C to exit')
            pause
        end

        % extract the result
        Vdot_opt = soln(1);

        % solve the reactor design equations using the optimum feed rate
        [t, nA, ~, ~, T, ~, V] = sbstr_model_variables(Vdot_opt);

        % calculate the conversion
        f_A = 100*(V_A*rho_A/M_A - nA(end))/(V_A*rho_A/M_A);

        % tabulate, show and save the results
        item = ["Optimum Feed Rate"; "Processing Time"...
            ; "Acetic Anhydride Conversion"];
        value = [Vdot_opt; t(end); f_A];
        units = ["cm^3^ min^-1^"; "min"; "%"];
        resultsTable = table(item, value, units);
        disp('')
        disp(resultsTable)
        disp('')
        writetable(resultsTable, 'ex_8_4_3_results.csv')

        % Generate plot data
        CA = nA./V;

        % generate graphs
        figure
        plot(t, T - 273.15, 'LineWidth', 2)
        xlabel("Processing Time (min)", "FontSize", 16)
        ylabel("Reacting Fluid Temperature (°C)", "FontSize", 16)
        yline(T_max - 273.15, ':r', 'LineWidth', 2)
        set(gca,'FontSize', 16)
        saveas(gcf,'ex_8_4_3_T_vs_t.png')

        figure
        plot(t, CA, 'LineWidth', 2)
        xlabel("Processing Time (min)", "FontSize", 16)
        ylabel("Acetic Acid Concentration (mol cm^-^3)", "FontSize", 16)
        set(gca,'FontSize', 16)
        saveas(gcf,'ex_8_4_3_CA_vs_t.png')

        % as a check, calculate the maximum temperature and processing time for 
        % volumetric flows slightly smaller and larger than the optimum

        Vdot = [0.95*Vdot_opt; 1.05*Vdot_opt];
        Tmax = [nan; nan];
        tproc = [nan; nan];

        for i = 1:2
            [t, ~, ~, ~, T, ~, ~] = sbstr_model_variables(Vdot(i));
            tproc(i) = t(end);
            Tmax(i) = max(T) - 273.15;
        end

        % tabulate, show, and save the results
        checkTable = table(Vdot, Tmax, tproc);
        checkTable.Properties.VariableNames = ["Feed Rate (cm^3 min^-1)"...
            ,"Maximum Temperature (°C)", "Processing Time (min)"];
        disp('')
        disp(checkTable)
        disp('')
        writetable(checkTable, 'ex_8_4_3_check.csv')
    end

    % perform the calculations
    items_of_interest()
end
        