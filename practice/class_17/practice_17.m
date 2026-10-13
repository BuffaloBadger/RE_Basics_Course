function practice_17()
% Calculations for the Class 17 Practice Assignment from REB, The Course
    %global constants available to all functions
    %given
    k01 = 1.75E9; %L^0.5 mol^-1.5 s^-1
    E1 = 35700; %cal/mol
    dH1 = -12900; %cal/mol
    CpA = 37.0; %cal/mol/K
    CpB = 34.0; %cal/mol/K
    CpY = 36.0; %cal/mol/K
    CpZ = 38.0; %cal/mol/K
    P0 = 1.0; %atm
    T0 = 400 + 273.15; %K
    yA0 = 0.5;
    yB0 = 0.5;
    fA = 0.5;
    V = 5.0; %L
    %known
    Rpv = 0.08206; %L-atm/mol/K
    R = 1.987; %cal/mol/K
    %calculated
    n0 = P0*V/Rpv/T0;
    nA0 = yA0*n0;
    nB0 = yB0*n0;
    nAf = nA0*(1-fA);

    %BSTR reactor function
    function [t, nA, nB, nY, nZ, T, P] = bstr_model_variables()
        %define the initial values
        ind_0 = 0;
        dep_0 = [nA0, nB0, 0, 0, T0, P0];

        %define the stopping criterion
        stopVar = 1;
        stopVal = nAf;

        %solve the BSTR design equations
        odes_are_stiff = true;
        [t, dep, success, message] = solve_ivodes(ind_0, dep_0, stopVar...
            , stopVal, @bstr_derivatives, odes_are_stiff);
        
        %check for solver issues
        if ~success
            disp('')
            disp(["BSTR model function issue: ", message])
            disp('Press any key to continue or CTRL-C to exit')
            pause;
        end
        
        %return the individual bstr model variables
        nA = dep(:,1);
        nB = dep(:,2);
        nY = dep(:,3);
        nZ = dep(:,4);
        T = dep(:,5);
        P = dep(:,6);
    end

    %BSTR derivatives function
    function ddt = bstr_derivatives(~, dep)
        %extract the dependent variables
        nA = dep(1);
        nB = dep(2);
        nY = dep(3);
        nZ = dep(4);
        T = dep(5);

        %calculate the additional unknowns
        k1 = k01*exp(-E1/(R*T));
        CA = nA/V;
        CB = nB/V;
        r1 = k1*CA*sqrt(CB);

        %form the mass matrix
        mm = zeros(6,6);

        %add the coefficients for the mole balances
        mm(1,1) = 1;
        mm(2,2) = 1;
        mm(3,3) = 1;
        mm(4,4) = 1;

        %add the coefficients for the energy balance
        mm(5,5) = nA*CpA + nB*CpB + nY*CpY + nZ*CpZ;
        mm(5,6) = -V;

        %add the coefficients for the differential gas law
        mm(6,1) = Rpv*T;
        mm(6,2) = Rpv*T;
        mm(6,3) = Rpv*T;
        mm(6,4) = Rpv*T;
        mm(6,5) = Rpv*(nA + nB + nY + nZ);
        mm(6,6) = -V;

        %form the right side vector
        rhs = [-r1*V; -r1*V; r1*V; r1*V; -r1*V*dH1; 0];

        %evaluate and return the derivatives
        ddt = mm\rhs;
    end

    % items of interest function
    function items_of_interest()
        %solve the BSTR design equations
        [t, ~, ~, ~, ~, T, P] = bstr_model_variables();

        %tabulate, show, and save the results
        item = ["Time"; "Temperature"];
        value = [t(end)/60; T(end) - 273.15];
        units = ["min"; "°C"];
        resultsTable = table(item, value, units);
        disp('')
        disp(resultsTable)
        disp('')
        writetable(resultsTable, 'results.csv');

        %check on pressure
        P_ratio_from_mole_balance = P(end)/P0;
        P_ratio_from_temperatures = T(end)/T0;

        disp('')
        disp(["Pressure ratio from mole balance: ", P_ratio_from_mole_balance])
        disp(["Pressure ratio from temperatures: ", P_ratio_from_temperatures])
    end

    % perform the calculations
    items_of_interest()
end