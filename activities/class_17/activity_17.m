function activity_17()
% Calculations for the Class 17 Learning Activity in REB, The Course
    % global constants available to all functions
    % given
    k01 = 265; % L /mol /min
    E1 = 73000; % J /mol
    dH298 = -165000; % J /mol
    P0 = 1000; % Torr
    T0 = 1225; % K
    V = 1; % basis
    Tf = [1235,1325]; % K
    % known
    Re = 8.314; % J /mol /K
    Rw = 62.36367; % L Torr /mol \K
    % calculated
    nA0 = P0*V/(Rw*T0);

    % given equations
    Cp_A = @(T) 28 + 0.05*T;
    Cp_Y = @(T) 26 + 0.01*T;
    Cp_Z = @(T) 30 + 0.005*T;

    % BSTR reactor function
    function [t, nA, nY, nZ, T, P] = bstr_model_variables(Tf)
        % define the initial values
        ind_0 = 0;
        dep_0 = [nA0; 0; 0; T0; P0];

        % define the stopping criterion
        stopVar = 4;
        stopVal = Tf;

        % solve the design equations
        odes_are_stiff = false;
        [t, dep, success, message] = solve_ivodes(ind_0, dep_0, stopVar...
            , stopVal, @bstr_derivatives, odes_are_stiff);
        
        % check for solver issues
        if ~success
            disp('')
            disp(["BSTR model function issue: ", message])
            disp('')
        end

        % return the individual dependent variables
        nA = dep(:,1);
        nY = dep(:,2);
        nZ = dep(:,3);
        T = dep(:,4);
        P = dep(:,5);
    end

    % BSTR derivatives function
    function ddt = bstr_derivatives(~, dep)
        % extract the dependent variables
        nA = dep(1);
        nY = dep(2);
        nZ = dep(3);
        T = dep(4);

        % calculate the additional unknowns
        dH_1 = dH298 -2*integral(Cp_A,298,T) + 2*integral(Cp_Y,298,T)...
            + integral(Cp_Z,298,T);
        r_1 = k01*exp(-E1/Re/T)*(nA/V)^2;

        % Create mass matrix, setting all elements to zero
        mass_matrix = zeros(5,5);

        % Add 1 on the diagonal for the first 3 rows
        mass_matrix(1,1) = 1.0;
        mass_matrix(2,2) = 1.0;
        mass_matrix(3,3) = 1.0;

        % Add the elements for the energy balance
        mass_matrix(4,4) = nA*Cp_A(T) + nY*Cp_Y(T) + nZ*Cp_Z(T);
        mass_matrix(4,5) = -V*Re/Rw;

        % Add the elements for the ideal gas law equation
        mass_matrix(5,1) = Rw*T;
        mass_matrix(5,2) = Rw*T;
        mass_matrix(5,3) = Rw*T;
        mass_matrix(5,4) = Rw*(nA + nY + nZ);
        mass_matrix(5,5) = -V;

        % Create right side vector
        rhs1 = -2*V*r_1;
        rhs2 = 2*V*r_1;
        rhs3 = V*r_1;
        rhs4 = -V*r_1*dH_1;
        rhs5 = 0.0;
        rhs = [rhs1; rhs2; rhs3; rhs4; rhs5];

        % Evaluate the derivatives to be returned
        ddt = mass_matrix\rhs;
    end

    % items of interest function
    function items_of_interest()
        % allocate storage for the deliverables
        tf = [nan, nan];
        fA = [nan, nan];

        % loop through the Tf values
        for n=1:2
            % solve the bstr design equations
            [t, nA, ~, ~, ~, ~] = bstr_model_variables(Tf(n));

            % save the final time and the conversion
            tf(n) = t(end);
            fA(n) = 100*(nA0 - nA(end))/nA0;
        end
        
        % tabulate, show, and save the results
        item = ["Final T (K)"; "Time (min)"; "Conversion (%)"];
        value = [Tf; tf; fA];
        resultsTable = table(item, value);
        disp('')
        disp(resultsTable)
        disp('')
        writetable(resultsTable,'results.csv');
    end

    % perform the calculations
    items_of_interest();
end
    