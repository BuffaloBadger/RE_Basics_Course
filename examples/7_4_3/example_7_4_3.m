function example_7_4_3()
%Calculations for Example 7.4.3 from REB, The Book
    % global constants available in all functions
    % given
    k0_1 = 2.59E9; % /min
    E_1 = 16500.; % cal /mol
    dH_1 = -22200.; % cal /mol
    CA_0 = 2.; % mol /L
    T_0 = 23  + 273.15; % K
    Cp = 440.; % cal /L /K
    V = 4.0; % L
    V_shell = 0.5; % L
    A_shell = 0.6; % ft^2
    U_shell = 1.13E4/60; % cal /ft^2 /min /K
    Tex_in = 20 + 273.15; % K
    rho_ex = 1.0; % g /cm^3
    Cp_ex = 1.0; % cal /g /K
    U_coil = 3.8E4/60; % cal /ft^2 /min /K
    A_coil = 0.23; % ft^2
    T_coil = 120 + 273.15; % K
    Tex_0 = 23 + 273.15; % K
    T_1 = 50 + 273.15; % K
    T_f = 25 + 273.15; % K
    t_turn = 25; % min
    % known
    R = 1.987; % cal /mol /K
    % calculated
    nA_0 = CA_0*V;

    % allocate global storage
    g_Vdot_ex = nan; % for the current coolant flow rate
    g_stage = 0; % for the current operational stage

    % BSTR reactor model
    function [t, nA, nZ, T, Tex] = bstr_model_variables(Vdot_ex)
        % set the current exchange fluid flow rate
        g_Vdot_ex = Vdot_ex;

        % set the current operational stage to 1
        g_stage = 1;

        % set the initial values for the first operational stage
        ind_0 = 0.0;
        dep_0 = [nA_0, 0.0, T_0, Tex_0];

        % define the stopping criterion for the first operational stage
        stopVar = 3;
        stopVal = T_1;
        
        % solve the design equations for the first operational stage
        odes_are_stiff = false;
        [t, dep, success, message] = solve_ivodes(ind_0, dep_0, stopVar...
            , stopVal, @bstr_derivatives, odes_are_stiff);
        
        % check for solver issues
        if ~success
            disp('')
            disp(["Stage 1 BSTR model issue: ", message])
            disp('')
        end

        % extract the dependent variable profiles
        nA = dep(:,1);
        nZ = dep(:,2);
        T = dep(:,3);
        Tex = dep(:,4);

        % set the current operational stage to 2
        g_stage = 2;

        % set the initial values for the second operational stage
        ind_0 = t(end);
        dep_0 = [nA(end), nZ(end), T(end), Tex(end)];

        % define the stopping criterion for the second operational stage
        stopVar = 3;
        stopVal = T_f;
        
        % solve the design equations for the second operational stage
        odes_are_stiff = false;
        [t2, dep, success, message] = solve_ivodes(ind_0, dep_0, stopVar...
            , stopVal, @bstr_derivatives, odes_are_stiff);
        
        % check for solver issues
        if ~success
            disp('')
            disp(["Stage 2 BSTR model issue: ", message])
            disp('')
        end

        % extract the dependent variable profiles
        nA2 = dep(:,1);
        nZ2 = dep(:,2);
        T2 = dep(:,3);
        Tex2 = dep(:,4);

        % concatenate the bstr model variables for the two stages
        t = [t; t2];
        nA = [nA; nA2];
        nZ = [nZ; nZ2];
        T = [T; T2];
        Tex = [Tex; Tex2];

    end

    % derivatives function
    function ddt = bstr_derivatives(~, dep)
        % extract the dependent variables
        nA = dep(1);
        T = dep(3);
        Tex = dep(4);

        % calculate the additional unknowns
        CA = nA/V;
        k_1 = k0_1*exp(-E_1/R/T);
        r_1 = k_1*CA;
        Qdot_shell = U_shell*A_shell*(Tex - T);
        Qdot_coil = U_coil*A_coil*(T_coil - T);

        % evaluate the derivatives
        dnAdt = -V*r_1;
        dnZdt = V*r_1;
        if g_stage == 1
            dTdt = (Qdot_shell + Qdot_coil - V*r_1*dH_1)/V/Cp;
            dTexdt = -Qdot_shell/rho_ex/V_shell/Cp_ex;
        else
            dTdt = (Qdot_shell - V*r_1*dH_1)/V/Cp;
            dTexdt = -(Qdot_shell + g_Vdot_ex*rho_ex*Cp_ex*(Tex-Tex_in))...
                /rho_ex/V_shell/Cp_ex;
        end
        
        % return the derivatives as a column vector
        ddt = [dnAdt; dnZdt; dTdt; dTexdt];
    end

    % deliverables function
    function deliverables()
        % choose a range of coolant flow rates
        %coolant_flows = np.linspace(100.0, 500.0, 100)
        coolant_flows = linspace(175.0, 200.0, 100);

        % allocate storage for the corresponding net rates
        net_rates = ones(1,100) * nan;

        % calculate the net rate for each coolant flow rate
        for i = 1:100
            % solve the reactor design equations
            [t, ~, nZ, ~, ~] = bstr_model_variables(coolant_flows(i));
            
            % calculate the net rate
            net_rates(i) = nZ(end)/(t(end) + t_turn);
        end

        % find the coolant flow where the net rate is maximized
        [rMax, iMax] = max(net_rates);
        Vdot_max = coolant_flows(iMax);

        % solve the reactor design equations using the optimum coolant flow
        [t, nA, nZ, T, ~] = bstr_model_variables(Vdot_max);
        
        % calculate the conversion vs time at the optimum coolant flow rate
        pct_conversion = 100*(nA_0 - nA)/nA_0;
        
        % tabulate, show, and save the results
        max_net_rate = nZ(end)/(t(end) + t_turn);
        item = ["Optimum Coolant Flow"; "Maximum Net Rate"];
        value = [Vdot_max; max_net_rate];
        units = ["cm^3 min^-1^"; "mol min^-1^"];
        resultsTable = table(item, value, units);
        disp('')
        disp(resultsTable)
        disp('')
        writetable(resultsTable, 'results.csv');

        % generate, show, and save the graphs
        figure; % net rate vs coolant flow
        plot(coolant_flows, net_rates, 'k', 'LineWidth', 2)
        xlabel("Coolant Flow (cm$^3$/min)", "FontSize", 14)
        ylabel("Net Rate (mol/min)","FontSize", 14)
        set(gca, "FontSize", 14)
        saveas(gcf, 'net_rate_vs_coolant_flow.png')

        figure(); % conversion profile
        plot(t, pct_conversion, 'k', 'LineWidth', 2)
        xlabel("Reaction Time (min)", "FontSize", 14)
        ylabel("Conversion (%)", "FontSize", 14)
        set(gca, "FontSize", 14)
        saveas(gcf, 'conversion_profile.png')

        figure(); % temperature profile
        plot(t,T - 273.15, 'k', 'LineWidth', 2)
        xlabel("Reaction Time (min)", "FontSize", 14)
        ylabel("Temperature (°C)", "FontSize", 14)
        set(gca, "FontSize", 14)
        saveas(gcf, 'temperature_profile.pdf')

    end

    % perform the calculations
    deliverables()
end