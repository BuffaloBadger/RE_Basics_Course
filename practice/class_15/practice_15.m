function practice_15()
% Calculations for the Class 15 Practice Assignment from REB, The Course
    % global constants available to all functions
    % given
    V = 50; % gal
    T0 = 60 + 459.67; % °R
    CA0 = 0.02; % lbmol /gal
    Tex = 212 + 459.67; % °R
    fA = 0.9;
    tf = 50; % min
    Cp = 1; % BTU /lb
    rho = 8.5; % lb /gal
    dH1 = - 5150; % BTU /lbmol
    k01 = 4630; % gal /lbmol /min
    E1 = 6580; % BTU /lbmol
    U = 37.3/60; % BTU /min /ft^2 /°R
    % known
    R = 1.986; % BTU /lbmol /°R
    % calculated
    nA0 = CA0 * V;
    nAf = nA0 * (1 - fA);

    % global variable for the heat transfer area
    g_A = nan;

    % BSTR reactor function
    function [t, nA, nZ, T] = bstr_model_variables(Area)
        % make the area available to the derivatives function
        g_A = Area;

        % set the initial values
        ind_0 = 0;
        dep_0 = [nA0, 0, T0];

        % set the stopping criterion
        stopVar = 0;
        stopVal = tf;

        % solve the BSTR design equations
        odes_are_stiff = false;
        [t, dep, success, message] = solve_ivodes(ind_0, dep_0, stopVar...
            , stopVal, @bstr_model_derivatives, odes_are_stiff);

        % check for solver issues
        if ~success
            disp('')
            disp(["BSTR model function issue: ", message])
            disp('Pre')
        end
        
        % return the individual bstr model variables
        nA = dep(:,1);
        nZ = dep(:,2);
        T = dep(:,3);
    end

    % BSTR derivatives function
    function ddt = bstr_model_derivatives(~, dep)
        % extract the dependent variables
        nA = dep(1);
        T = dep(3);

        % calculate the additional unknowns
        CA = nA / V;
        k = k01 * exp(- E1 / (R * T));
        r = k * CA^2;
        Qdot = U * g_A * (Tex - T);

        % calculate the derivatives
        dnA_dt = - r * V;
        dnZ_dt = r * V;
        dT_dt = (Qdot - r*V*dH1) / (rho*V*Cp);

        % return the derivatives as a column vector
        ddt = [dnA_dt; dnZ_dt; dT_dt];
    end

    % coupled unknown residual function
    function epsilon = coupled_unknown_residual(Area)
        % calculate the bstr model variables
        [~, nA, ~, ~] = bstr_model_variables(Area);

        % calculate and return the residual
        epsilon = nA(end) - nAf;
    end

    % deliverables function
    function deliverables()
        % guess the heat transfer area
        A_guess = 1.0;

        % calculate the heat transfer area
        [soln, success, message] = solve_ates(@coupled_unknown_residual...
            , A_guess);
        Area = soln(1);

        % check for solver issues
        if ~success
            disp('')
            disp(["Coupled unknown residual function issue: ", message])
            disp('Press any key to continue or CTRL-C to exit')
            pause
        end
        
        % solve the BSTR design equations
        [t, nA, ~, T] = bstr_model_variables(Area);

        % tabulate, show, and save the results
        item = ["Area"; "Reaction Time"; "Final Conversion"...
            ; "Final Temperature"];
        value = [Area; t(end); 100*(nA0 - nA(end)) / nA0; T(end)  - 459.67];
        units = ["ft^2"; "min"; "%"; "°F"];
        resultsTable = table(item, value, units);
        disp('')
        disp(resultsTable)
        disp('')
        writetable(resultsTable,'results.csv');

        % plot the temperature profile for discussion
        figure
        plot(t, T - 459.67, 'k', 'LineWidth', 2)
        xlabel('Time (min)', 'FontSize', 14)
        ylabel('Temperature (°F)', 'FontSize', 14)
        set(gca, 'FontSize', 14)
        saveas(gcf,'T_vs_t.png')
    end

    % perform the calculations
    deliverables();
end
    