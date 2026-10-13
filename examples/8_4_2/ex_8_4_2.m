function ex_8_4_2()
% Calculations for Example 8.4.2 from REB, The Book
    % global constants available to all functions
    % given
    k0_1 = 1.83E12; % L/mol/h
    E_1 = 18000.0; % cal/mol
    k0_2 = 5.08E13; % L/mol/h
    E_2 = 20500.0; % cal/mol
    dH_1 = -9000.0; % cal/mol
    dH_2 = -7800.0; % cal/mol
    Cp = 863.0; % cal/L/K
    CA_0 = 2.0; % mol/L
    CB_in = 0.5; % mol/L
    T_0 = 40 + 273.15; % K
    T_in = T_0;
    P = 1.0; % atm
    V_A = 2000.0; % L
    V_B = 8000.0; % L
    tf = 8.0; % h
    t1_range = [1, 3, 5, 7];
    % known
    Re = 1.987; % cal/mol/K
    Rw = 0.08206; % L-atm/mol/K
    % calculated
    nA_0 = CA_0*V_A;

    % define a global variables for the current value of t1 and the current protocol stage
    g_t1 = nan;
    g_stage = -1;

    % SBSTR reactor function
    function [t, nA, nB, nD, nU, T, V] = sbstr_model_variables(t1)
        % make t1 available to the derivatives function
        g_t1 = t1;

        % set the current operational protocol stage;
        g_stage = 1;

        % set the initial values
        ind_0 = 0;
        dep_0 = [nA_0, 0, 0, 0, T_0, V_A];

        % set the stopping criterion
        stopVar = 0;
        stopVal = t1;

        % solve the IVODE design equations
        odes_are_stiff = false;
        [t1, dep1, success, message] = solve_ivodes(ind_0, dep_0, stopVar...
            , stopVal, @sbstr_derivatives, odes_are_stiff);
        
        % check for solver issues
        if ~success
            disp('')
            disp(["SBSTR model stage 1 issue: ", message])
            disp('Press any key to continue or CTRL-C to exit')
            pause
        end
        
        % repeat for the second stage
        g_stage = 2;
        ind_0 = t1(end);
        dep_0 = [dep1(end,1), dep1(end,2), dep1(end,3), dep1(end,4)...
            , dep1(end,5), dep1(end,6)];
        stopVal = tf;
        [t2, dep2, success, message] = solve_ivodes(ind_0, dep_0, stopVar...
            , stopVal, @sbstr_derivatives, odes_are_stiff);
        if ~success
            disp('')
            disp(["SBSTR model stage 2 issue: ", message])
            disp('Press any key to continue or CTRL-C to exit')
            pause
        end
        
        % concatenate the stages
        t = [t1; t2];
        nA = [dep1(:,1); dep2(:,1)];
        nB = [dep1(:,2); dep2(:,2)];
        nD = [dep1(:,3); dep2(:,3)];
        nU = [dep1(:,4); dep2(:,4)];
        T = [dep1(:,5); dep2(:,5)];
        V = [dep1(:,6); dep2(:,6)];
    end

    % SBSTR derivatives function
    function ddt = sbstr_derivatives(~, dep)
        % extract the dependent variables
        nA = dep(1);
        nB = dep(2);
        T = dep(5);
        V = dep(6);

        % calculate the additional unknowns
        if g_stage == 1
            Vdot_in = V_B/g_t1;
        else
            Vdot_in = 0.0;
        end
        k_1 = k0_1*exp(-E_1/Re/T);
        k_2 = k0_2*exp(-E_2/Re/T);
        CA = nA/V;
        CB = nB/V;
        r_1 = k_1*CA*CB;
        r_2 = k_2*CB^2;
        nDotB_in = Vdot_in*CB_in;

        % evaluate the derivatives
        dnAdt = -V*r_1;
        dnBdt = nDotB_in - V*(r_1 + 2*r_2);
        dnDdt = V*r_1;
        dnUdt = V*r_2;
        dTdt = (-Vdot_in*Cp*(T - T_in) - V*r_1*dH_1 - V*r_2*dH_2...
            + P*Vdot_in*Re/Rw)/(V*Cp);
        dVdt = Vdot_in;

        % return the derivatives as a column vector
        ddt = [dnAdt; dnBdt; dnDdt; dnUdt; dTdt; dVdt];
    end

    % items of interest function
    function items_of_interest()
        % allocate storage for the results
        nParams = length(t1_range);
        f_B = ones(1, nParams)*nan;
        S_DU = ones(1, nParams)*nan;
        Y_DB = ones(1, nParams)*nan;

        % loop through the t1 values
        for i = 1:nParams
            % solve the design equations
            [~, ~, nB, nD, nU, ~, ~] = sbstr_model_variables(t1_range(i));

            % calculate and save the quantities of interest
            f_B(i) = 100*(V_B*CB_in - nB(end))/(V_B*CB_in);
            S_DU(i) = nD(end)/nU(end);
            Y_DB(i) = 100*nD(end)/(V_B*CB_in);
        end
        
        % tabulate, show and save the results
        resultsTable = table(t1_range.', f_B.', S_DU.', Y_DB.');
        resultsTable.Properties.VariableNames = ["Stage 1 Duration (h)",...
            "Conversion (%)", "Selectivity D/U", "Yield D-B (%)"];
        disp('')
        disp(resultsTable)
        disp('')
        writetable(resultsTable, 'ex_8_4_2_results.csv')
    end

    % perform the calculations
    items_of_interest()
end

        