function example_7_4_4()
% Calculations for Example 7.4.4 from REB, The Book"""
    % global constants available to all functions
    % given
    PA_0 = 1.; % atm
    PB_0 = 2.; % atm
    T_0 = 25 + 273.14; % K
    V = 2.; % L
    Aex = 600.; % cm^2
    Uex = 0.6; % cal /cm^2 /min /K
    Tex = 30 + 273.15; % K
    k_0_1 =  3.34E9; % mol /cm^3 /min /atm^2
    k_0_2 = 4.99E9; % mol /cm^3 /min /atm^2
    E_1 = 20.5E3; % cal /mol
    E_2 = 21.8E3; % cal /mol
    dH_1 = -6300.; % cal /mol
    dH_2 = -6900.; % cal /mol
    Cp_A = 7.4; % cal /mol /K
    Cp_B = 8.6; % cal /mol /K
    Cp_D = 10.7; % cal /mol /K
    Cp_Z = 5.2; % cal /mol /K
    Cp_U = 10.3; % cal /mol /K
    % known
    Re = 1.987; % cal /mol /K
    Rw = 82.057; % cm^3 atm /mol /K
    % calculated
    nA_0 = PA_0*V/Rw/T_0;
    nB_0 = PB_0*V/Rw/T_0;
    P_0 = PA_0 + PB_0;

    % BSTR model function
    function [t, nA, nB, nD, nZ, nU, T, P] = bstr_model_variables(t_f)
        % set the initial values
        ind_0 = 0.0;
        dep_0 = [nA_0, nB_0, 0.0, 0.0, 0.0, T_0, P_0];

        % define the stopping criterion
        stopVar = 0;
        stopVal = t_f;
        
        % solve the IVODEs
        odes_are_stiff = true;
        [t, dep, success, message] = solve_ivodes(ind_0, dep_0, stopVar...
            , stopVal, @bstr_derivatives, odes_are_stiff);

        % check for solver issues
        if ~success
            disp('')
            disp(["BSTR model function issue: ", message])
            disp('')
        end

        % extract the dependent variable profiles to be returned
        nA = dep(:,1);
        nB = dep(:,2);
        nD = dep(:,3);
        nZ = dep(:,4);
        nU = dep(:,5);
        T = dep(:,6);
        P = dep(:,7);
    end

    % derivatives function
    function ddt = bstr_derivatives(~, dep)
        % extract the dependent variables
        nA = dep(1);
        nB = dep(2);
        nD = dep(3);
        nZ = dep(4);
        nU = dep(5);
        T = dep(6);

        % calculate the additional unknowns
        PA = nA*Rw*T/V;
        PB = nB*Rw*T/V;
        PD = nD*Rw*T/V;
        k_1 = k_0_1*exp(-E_1/Re/T);
        k_2 = k_0_2*exp(-E_2/Re/T);
        r_1 = k_1*PA*PB;
        r_2 = k_2*PD*PB;
        Qdot = Uex*Aex*(Tex - T);

        % Create mass matrix, setting all elements to zero
        mass_matrix = zeros(7,7);

        % Add 1 on the diagonal for the first 5 rows
        mass_matrix(1,1) = 1.0;
        mass_matrix(2,2) = 1.0;
        mass_matrix(3,3) = 1.0;
        mass_matrix(4,4) = 1.0;
        mass_matrix(5,5) = 1.0;

        % Add the elements for the energy balance
        mass_matrix(6,6) = nA*Cp_A + nB*Cp_B + nD*Cp_D + nZ*Cp_Z + nU*Cp_U;
        mass_matrix(6,7) = -V*Re/Rw;

        % Add the elements for the ideal gas law equation
        mass_matrix(7,1) = Rw*T;
        mass_matrix(7,2) = Rw*T;
        mass_matrix(7,3) = Rw*T;
        mass_matrix(7,4) = Rw*T;
        mass_matrix(7,5) = Rw*T;
        mass_matrix(7,6) = Rw*(nA + nB + nD + nZ +nU);
        mass_matrix(7,7) = -V;

        % Create right side vector
        rhs1 = -r_1*V;
        rhs2 = (-r_1 -r_2)*V;
        rhs3 = (r_1-r_2)*V;
        rhs4 = (r_1+r_1)*V;
        rhs5 = r_2*V;
        rhs6 = Qdot -(r_1*dH_1 + r_2*dH_2)*V;
        rhs7 = 0.0;
        rhs = [rhs1; rhs2; rhs3; rhs4; rhs5; rhs6; rhs7];

        % Evaluate the derivatives to be returned
        ddt = mass_matrix\rhs;
    end

    % deliverables function
    function deliverables()
        % choose a large final time
        tf = 60; % min

        % solve the BSTR design equations
        [t, nA, ~, nD, ~, ~, ~, ~] = bstr_model_variables(tf);

        % calculate the corresponding yield
        yield_DA = nD/nA_0;

        % find the reaction time where the yield is maximized
        [yield_opt, i_opt] = max(yield_DA);

        % find the corresponding reaction time
        t_opt = t(i_opt);

        % check that the optimum time is between 0 and t_f
        if (t_opt == 0) || (t_opt == tf)
            disp('')
            disp('The range of reaction times needs to be wider')
            disp('')
        end

        % calculate the conversion
        f_opt = 100*(nA_0 - nA(i_opt))/nA_0;

        % tabulate the results
        item = ["t opt"; "yield D/A"; "conversion"];
        value = [t_opt; yield_opt; f_opt];
        units = ["min"; "mol D per initial mol A"; "%"];
        resultsTable = table(item, value, units);

        % display the results
        disp(' ')
        disp(resultsTable)
        disp(' ')

        % save the results
        writetable(resultsTable, 'results.csv');
        
        % plot, show, and save the yield vs. the reaction time
        figure
        plot(t, yield_DA, 'k', 'LineWidth', 2)
        xlabel("Reaction Time (min)", "FontSize", 14)
        ylabel("Yield (mol D per initial mol A)", "FontSize", 14)
        set(gca, 'FontSize', 14)
        saveas(gcf, 'yield_vs_t.png')
    end

    % perform the calculations
    deliverables()
end