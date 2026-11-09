function ex_9_6_2()
% Calculations for Example 9.6.2 from REB, The Book
    % global constants available to all functions
    % given
    Vpfr = 40E3; % cm^3
    CA_in = 1E-3; % mol /cm^3
    CB_in = 1.2E-3; % mol /cm^3
    Vdot_in = 75E3; % cm^3 /min
    k0_1 = 8.72E8; % cm^3 /mol /min
    E_1 = 7200; % cal /mol
    dH_1 = -10700; % cal /mol
    Cp = 1.0; % cal /g /K
    rho = 1.0; % g /cm^3
    f_A = 0.95;
    % known
    R = 1.987; % cal /mol /K
    % calculated
    nA_in = Vdot_in*CA_in;
    nB_in = Vdot_in*CB_in;
    nA_out = nA_in*(1-f_A);

    % PFR reactor function
    function [V, nA, nB, nY, nZ, T] = pfr_model_variables(T_in)
        % set the initial values
        ind_0 = 0.0;
        dep_0 = [nA_in, nB_in, 0.0, 0.0, T_in];

        % define the stopping criterion
         stopVar = 0;
         stopVal = Vpfr;
        
        % solve the IVODEs
        odes_are_stiff = false;
        [V, dep, success, message] = solve_ivodes(ind_0, dep_0,  stopVar...
            ,  stopVal, @pfr_derivatives, odes_are_stiff);

        % check for solver issues
        if ~(success)
             disp('')
             disp(["PFR model function issue: ", message])
             disp('Press any key to continue or CTRL-C to exit')
            pause
        end

        % extract the individual dependent variables to return
        nA = dep(:,1);
        nB = dep(:,2);
        nY = dep(:,3);
        nZ = dep(:,4);
        T = dep(:,5);
    end

    % PFR derivatives function
    function ddz = pfr_derivatives(~, dep)
        % extract the dependent variables
        nA = dep(1);
        nB = dep(2);
        T = dep(5);

        % calculate the additional unknowns
        k_1 = k0_1*exp(-E_1/R/T);
        CA = nA/Vdot_in;
        CB = nB/Vdot_in;
        r_1 = k_1*CA*CB;

        % evaluate the derivatives
        dnAdV = -r_1;
        dnBdV = -r_1;
        dnYdV = r_1;
        dnZdV = r_1;
        dTdV = -r_1*dH_1/(Vdot_in*rho*Cp);

        % return the derivatives as a column vector
        ddz = [dnAdV; dnBdV; dnYdV; dnZdV; dTdV];
    end

    % coupled unknown residual function
    function epsilon = coupled_unknown_residual(guess)
        % extract the guess
        T_in = guess(1);

        % solve the reactor design equations
        [~, nA, ~, ~, ~, ~] = pfr_model_variables(T_in);

        % extract the calculated final value of nA
        nA_f = nA(end);

        % evaluate  and return the residual
        epsilon = nA_out - nA_f;
    end

    % items of interest function
    function items_of_interest()
        % set an initial guess for T_in
        initial_guess = 25 + 273.15;

        % solve the implict equation for [missing constant]
        [soln, success, message] = solve_ates(@coupled_unknown_residual...
            , initial_guess);

        % check for solver issues
        if ~(success)
            disp('')
            disp(["Issue solving for the coupled unknown: ", message])
            disp('Press any key to continue or CTRL-C to exit')
            pause
        end

        % extract the result
        T_in = soln(1);

        % solve the reactor design equations
        [~, ~, ~, ~, ~, T] = pfr_model_variables(T_in);

        % calculate the other quantities of interest
        T_f = T(end);

        % tabulate, show and save the results
        item = ["Inlet T"; "Outlet T"];
        value = [T_in - 273.15; T_f - 273.15];
        units = ["°C"; "°C"];
        resultsTable = table(item, value, units);
         disp(' ')
         disp(resultsTable)
         disp('')
         writetable(resultsTable, 'ex_9_6_2_results.csv');
    end

    % perform the calculations
    items_of_interest()
end
        