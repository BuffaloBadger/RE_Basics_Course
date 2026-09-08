function activity_10()
% Calculations for the Class 10 Learning Activity in REB, The Course
    % global constants available to all functions
    % given
    V = 50; % gal
    CA_in = 0.014; % lbmolA /gal
    CB_in = 0.020; % lbmolB /gal
    T_in = 70 + 459.67; % °R
    Cp = 65.0 * 0.1337; % BTU /gal /°R
    Tex = 220 + 459.67; % °R
    U = 60; % BTU /ft^2 /°R /h
    A = 13; % ft^2
    dH1 = 45000.; % BTU /lbmol
    dH2 = 39500.; % BTU /lbmol
    k01 = 1.6E18; % gal /lbmol /h
    E1 = 46000.; % BTU /lbmol
    k02 = 4.5E18; % gal /lbmol /h
    E2 = 54000.; % BTU /lbmol
    % known
    R = 1.987; % BTU/lbmol
    % calculated

    % allocate global storage for the current value of Vdot_in
    g_Vdot_in = nan;

    % CSTR model function
    function [nA, nB, nD, nU, nZ, T] = cstr_model_variables(Vdot_in, init_guess)
        % make Vdot_in available to the residuals function
        g_Vdot_in = Vdot_in;

        % solve the design equations
        [soln, success, message] = solve_ates(@cstr_residuals, init_guess);

        % check for solver issues
        if ~success
            disp(' ')
            disp(['    CSTR model function issue: ', message])
            disp(' ')
        end
        
        % extract and return the results
        nA = soln(1);
        nB = soln(2);
        nD = soln(3);
        nU = soln(4);
        nZ = soln(5);
        T = soln(6);
    end

    % CSTR residuals function
    function epsilon = cstr_residuals(guess)
        % extract the individual guesses
        nA = guess(1);
        nB = guess(2);
        nD = guess(3);
        nU = guess(4);
        nZ = guess(5);
        T = guess(6);

        % calculate the additional unknowns
        nA_in = CA_in * g_Vdot_in;
        nB_in = CB_in * g_Vdot_in;
        k1 = k01 * exp(-E1 / (R * T));
        k2 = k02 * exp(-E2 / (R * T));
        Vdot = g_Vdot_in;
        CA = nA / Vdot;
        CB = nB / Vdot;
        CD = nD / Vdot;
        r1 = k1 * CA * CB;
        r2 = k2 * CD * CB;
        Qdot = U * A * (Tex - T);
        
        % evaluate and return the CSTR design equation residuals
        epsilon_1 = nA_in - nA - r1 * V;
        epsilon_2 = nB_in - nB - r1 * V - r2 * V;
        epsilon_3 = -nD + r1 * V - r2 * V;
        epsilon_4 = -nU + r2 * V;
        epsilon_5 = -nZ + r1*V + r2*V;
        epsilon_6 = g_Vdot_in*Cp*(T - T_in) - Qdot + r1*V*dH1 + r2*V*dH2;

        % return the residuals as an array
        epsilon = [epsilon_1; epsilon_2; epsilon_3; epsilon_4; epsilon_5; epsilon_6];
    end

    % deliverables function
    function deliverables()
        % define the range of Vdot_in values to evaluate
        VFR_range = linspace(35, 50, 100);

        % allocate storage for the corresponding values of nD
        nD_range = ones(1,100) * nan;

        % define an initial guess for the first Vdot_in in the range
        %guess = [CA_in*VFR_range(1); CB_in*VFR_range(1); 0.0; 0.0; 0.0...
        %; T_in + 5.0];   did not converge
        guess = [CA_in*VFR_range(1); CB_in*VFR_range(1); 0.0; 0.0; 0.0...
        ; T_in + 100.0];

        % loop through the range of Vdot_in values
        for i = 1:100
            % solve the cstr design equations
            [nA, nB, nD, nU, nZ, T] = cstr_model_variables(VFR_range(i), guess);

            % save nD
            nD_range(i)= nD;

            % save the solution to use as the next guess
            guess = [nA; nB; nD; nU; nZ; T];
        end

        % plot the outlet molar flow of D versus the inlet volumetric flow rate
        figure
        plot(VFR_range, nD_range, 'k', 'LineWidth', 2)
        set(gca, 'FontSize', 14);
        xlabel('Inlet Volumetric Flow Rate (gal/h)', 'FontSize', 14)
        ylabel('Outlet Molar Flow of D (lbmol/h)', 'FontSize', 14)
        saveas(gcf,"act_10_nD_vs_VdotIn.pdf")

        % find the maximum flow rate of D
        [~, i_max] = max(nD_range);

        % identify the corresponding inlet volumetric flow rate
        Vdot_opt = VFR_range(i_max);

        % solve the cstr design equations at the optimum Vdot_in
        [nA, ~, ~, ~, ~, T] = cstr_model_variables(Vdot_opt, guess);

        % calculate the conversion
        nA_in = CA_in*Vdot_opt;
        fA = 100*(nA_in - nA)/nA_in;

        % tabulate, show, and save the results
        item = ["Optimum Vdot_in"; "Outlet T"; "Conversion"];
        value = [Vdot_opt; T - 459.67; fA];
        units = ["gal/h"; "°F"; "%"];
        resultsTable = table(item, value, units);
        disp(' ')
        disp(resultsTable)
        disp(' ')
        writetable(resultsTable,'act_10_results.csv')
    end

    % perform the calculations
    deliverables()
end