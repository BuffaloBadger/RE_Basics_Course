function example_7_4_2()
% Calculations for Example 7.4.2 from REB, The Book
    % global constants available to all functions
    % given
    V = 10.0E3; % cm^3
    Vex = 1.4E3; % cm^3
    Uex = 138.; % cal /ft^2 /min /K
    Aex = 1200./929.; % ft^2
    Tex_0 = 40. + 273.15; % K
    Tex_in = 40. + 273.15; % K
    mDot_ex = 100.; % g /min
    rho = 1.0; % g /cm^3
    rho_ex = 1.0; % g /cm^3
    Cp = 1.0; % cal /g /K
    Cp_ex = 1.0; % cal /g /K
    CA_0 = 5.0E-3; % mol /cm^3
    CB_0 = 7.0E-3; % mol /cm^3
    dH_1 = -16.7E3; % cal /mol
    dH_2 = -14.3E3; % cal /mol
    k0_1 = 9.74E12; % cm^3 /mol /min
    E_1 = 20.1E3; % cal /mol
    k0_2 = 2.38E13; % /min
    E_2 = 25.3E3; % cal /mol
    t_f = 30.; % min
    fA_f = 0.45;
    % known
    R = 1.987; % cal /mol /K
    % calculated
    nA_0 = CA_0*V;
    nB_0 = CB_0*V;
    nA_f = nA_0*(1 - fA_f);

    % BSTR model function
    function [t, nA, nB, nX, nY, nZ, T, Tex] = bstr_model_variables(T_0)
        % set the initial values
        ind_0 = 0.0;
        dep_0 = [nA_0, nB_0, 0.0, 0.0, 0.0, T_0, Tex_0];

        % define the stopping criterion
        stopVar = 0;
        stopVal = t_f;
        
        % solve the IVODEs
        odes_are_stiff = false;
        [t, dep, success, message] = solve_ivodes(ind_0, dep_0, stopVar...
            , stopVal, @bstr_derivatives, odes_are_stiff);

        % check that a solution was found
        if ~success
            disp('')
            disp(["BSTR function issue: ", message])
            disp('')
        end

        % extract the individual dependent variables
        nA = dep(:,1);
        nB = dep(:,2);
        nX = dep(:,3);
        nY = dep(:,4);
        nZ = dep(:,5);
        T = dep(:,6);
        Tex = dep(:,7);
    end

    % BSTR derivatives function
    function ddt = bstr_derivatives(~, dep)
        % extract the dependent variables for this integration step
        nA = dep(1);
        nB = dep(2);
        T = dep(6);
        Tex = dep(7);

        % calculate the additional unknowns
        CA = nA/V;
        CB = nB/V;
        k_1 = k0_1*exp(-E_1/R/T);
        r_1 = k_1*CA*CB;
        k_2 = k0_2*exp(-E_2/R/T);
        r_2 = k_2*CA;
        Qdot = Uex*Aex*(Tex - T);

        % evaluate the derivatives
        dnAdt = -(r_1 + r_2)*V;
        dnBdt = -r_1*V;
        dnXdt = r_1*V;
        dnYdt = r_1*V;
        dnZdt = r_2*V;
        dTdt = (Qdot-(r_1*dH_1 + r_2*dH_2)*V)/rho/V/Cp;
        dTexdt = (-Qdot - mDot_ex*Cp_ex*(Tex - Tex_in))/rho_ex/Vex/Cp_ex;

        % return the derivatives as a column vector
        ddt = [dnAdt; dnBdt; dnXdt; dnYdt; dnZdt; dTdt; dTexdt];
    end

    % coupled unknown residual function
    function epsilon = coupled_unknown_residual(T0_guess)
        % solve the BSTR design equations
        [~, nA, ~, ~, ~, ~, ~, ~] = bstr_model_variables(T0_guess);

        % extract the calculated final molar amount of A
        nA_f_calc = nA(end);

        % evaluate the residual to be returned
        epsilon = nA_f_calc - nA_f;
    end

    % deliverables function
    function deliverables()
        % set initial guess for T0
        T0_guess = Tex_in + 20;

        % solve the coupled unknown equation for T0
        [soln, success, message] = solve_ates(@coupled_unknown_residual...
            , T0_guess);

        % check for solver issues
        if ~success
            disp('')
            disp(["Issue when solving for the coupled unknown: ", message])
            disp('')
        end

        % extract the result
        T_0 = soln(1);

        % solve the reactor design equations
        [t, nA, ~, nX, ~, nZ, T, Tex] = bstr_model_variables(T_0);

        % check the final values
        disp('')
        disp('Check:')
        disp(["final time (should be 30): ", t(end)]);
        disp(["final conversion (should be 45%): ", 100*(nA_0 - nA(end))/nA_0])

        % calculate the other quantities of interest
        T_0 = T_0 - 273.15;
        T_f = T(end) - 273.15;
        Tex_out = Tex(end) - 273.15;
        sel_X_Z = nX(end)/nZ(end);

        % tabulate the results
        item = ["T0"; "Tf"; "Tex out"; "Selectivity X/Z"];
        value = [T_0; T_f; Tex_out; sel_X_Z];
        units = ["°C"; "°C"; "°C"; "mol X per mol Z"];
        resultsTable = table(item, value, units);

        % display the results
        disp(' ')
        disp(resultsTable)
        disp('')

        % save the results
        writetable(resultsTable,'results.csv')
    end

    % perform the calculations
    deliverables()
end