function act_20()
% Calculations for the Class 20 Learning Activity from REB, The Book
    % global constants available to all functions
    % given
    P = 1.0; % atm
    V0 = 550.; % gal
    T0 = 100 + 459.7; % °R
    CA0 = 0.015; % lbmol/gal
    VB = 2000.; % gal
    Tin = 100 + 459.7; % °R
    CBin = 0.004; % lbmol/gal
    fB_min = 0.9;
    tTurn = 0.5; % h
    k01 = 4.193e+11; % gal/lbmol/h
    E1 = 32400; % BTU/lbmol
    dH1 = -16190; % BTU/lbmol
    k02 = 6.087e+14; % gal/lbmol/h
    E2 = 36900; % BTU/lbmol
    dH2 = -14030; % BTU/lbmol
    Cp = 7.2; % BTU/gal/°R
    % known
    R = 1.987; % BTU/lbmol/°R
    Rpv = 5.4584; % gal-atm/lbmol/°R
    % calculated
    nA0 = CA0*V0;

    % global variable for the current volumetric flow rate
    g_Vdot = nan;

    % SBSTR reactor function
    function [t, nA, nB, nD, nU, T, V] = sbstr_model_variables(Vdot)
        % set the volumetric flow rate for stage 1
        g_Vdot = Vdot;

        % set the initial values for the first stage
        ind_0 = 0;
        dep_0 = [nA0, 0, 0, 0, T0, V0];

        % set the stopping criterion
        stopVar = 6;
        stopVal = V0 + VB;

        % solve the design equations for stage 1
        odes_are_stiff = false;
        [t, dep, success, message] = solve_ivodes(ind_0, dep_0, stopVar...
            , stopVal, @sbstr_derivatives, odes_are_stiff);
        
        % extract the individual dependent variables
        nA = dep(:,1);
        nB = dep(:,2);
        nD = dep(:,3);
        nU = dep(:,4);
        T = dep(:,5);
        V = dep(:,6);
        
        % check for solver issues
        if ~success
            disp('')
            disp(["SBSTR Model stage 1 issue: ", message])
            disp('Press any key to continue or CTRL-C to exit')
            pause
        end
        
        % calculate the conversion of B during stage 1
        fB1 = 100*(VB*CBin - nB(end))/(VB*CBin);

        % if the conversion is below 90% continue with stage 2
        if fB1 < 90
            
            % set the volumetric flow rate for stage 2
            g_Vdot = 0;

            % set the initial values equal to the final values from stage 1
            ind_0 = t(end);
            dep_0 = [nA(end), nB(end), nD(end), nU(end), T(end), V(end)];

            % set the stopping criterion for stage 2
            stopVar = 2;
            stopVal = VB*CBin*(1 - fB_min);
            
            % solve the design equations for stage 2
            odes_are_stiff = false;
            [t2, dep, success, message] = solve_ivodes(ind_0, dep_0, stopVar...
                , stopVal, @sbstr_derivatives, odes_are_stiff);
            
            % combine the results for the two stages
            t = [t; t2];
            nA = [nA; dep(:,1)];
            nB = [nB; dep(:,2)];
            nD = [nD; dep(:,3)];
            nU = [nU; dep(:,4)];
            T = [T; dep(:,5)];
            V = [V; dep(:,6)];
            
            % check for solver issues
            if ~success
                disp('')
                disp(["SBSTR Model stage 2 issue: ", message])
                disp('Press any key to continue or CTRL-C to exit')
                pause
            end
        end
    end

    % SBSTR derivatives function
    function ddt = sbstr_derivatives(~, dep)
        % Extract the dependent variables
        nA = dep(1);
        nB = dep(2);
        T = dep(5);
        V = dep(6);

        % calculate the additional unknowns
        k1 = k01 * exp(-E1/(R*T));
        k2 = k02 * exp(-E2/(R*T));
        CA = nA / V;
        CB = nB / V;
        r1 = k1 * CA * CB;
        r2 = k2 * CB^2;
        nB_in = g_Vdot * CBin;
        
        % Calculate the derivatives
        dnAdt = -V*r1;
        dnBdt = nB_in - V*(r1 + 2*r2);
        dnDdt = V*r1;
        dnUdt = V*r2;
        dTdt = (-g_Vdot*Cp*(T-Tin)-V*r1*dH1 -V*r2*dH2 + P*g_Vdot*R/Rpv)/V/Cp;
        dVdt = g_Vdot;
        
        % return the derivatives as a column vector
        ddt = [dnAdt; dnBdt; dnDdt; dnUdt; dTdt; dVdt];
    end

    % items of interest function
    function items_of_interest()
        % set a range of values for the volumetric flow rate
        %Vdot_range = linspace(20, 40, 100)
        Vdot_range = linspace(5, 20, 100); % narrow the range

        % allocate storage for the corresponding net rates
        rNet = ones(100,1)*nan;

        % for each value of the volumetric flow rate
        for i = 1:100
            % solve the design equations
            [t, ~, ~, nD, ~, ~, ~] = sbstr_model_variables(Vdot_range(i));

            % calculate and save the corresponding net rate
            rNet(i) = nD(end)/(t(end) + tTurn);
        end
        
        % find the index of the maximum net rate and corresponding Vdot
        [max_rNet, iMax] = max(rNet);
        Vdot_opt = Vdot_range(iMax);
        
        % tabulate, show and save the results
        item = ["Optimum Feed Rate"; "Net Rate"];
        value = [Vdot_opt; max_rNet];
        units = ["gal/h"; "lbmol/h"];
        resultsTable = table(item, value, units);
        disp('')
        disp(resultsTable)
        disp('')
        writetable(resultsTable, 'act_20_results.csv')

        % plot the net rate vs. the feed rate
        figure
        plot(Vdot_range, rNet, 'LineWidth', 2)
        xlabel('Feed Rate (gal/h)', 'FontSize', 18)
        ylabel('Net Rate (lbmol/h)','FontSize', 18)
        set(gca,'FontSize', 18)
        saveas(gcf,'act_20_rNet_vs_VFR.png')
    end

    % perform the calculations
    items_of_interest()
end
