function activity_16()
% Calculations for the Class 16 Learning Activity from REB, The Course"""
    % global constants available to all functions
    % given
    V = 10; % L
    U1 = 76; % cal /h /cm^2 /K
    A1 = 600; % cm^2
    U2 = 91; % cal /h /cm^2 /K
    A2 = 375; % cm^2
    Tex_1 = 115 + 273.15; % K
    Tex_2 = 30 + 273.15; % K
    CA0 = 5; % mol / L
    T0 = 25 + 273.15; % K
    Cp = 1000; % cal /L /K
    dH1 = 20000; % cal / mol
    k01 = 2.85E15; % 1/h
    E1 = 25000; % cal / mol
    t1 = 1.0; % h
    Tf = 35 + 273.15; % K
    t_turn = 0.5; % h
    % known
    R = 1.987; % cal / mol / K
    % calculated
    nA0 = CA0 * V;

    % global variable indicating the stage in the protocol
    g_stage = -1;

    % BSTR reactor function
    function [t, nA, nZ, T] = bstr_model_variables()
        % set the protocol stage to 1
        g_stage = 1;

        % set the initial values for stage 1 of the protocal
        ind_0 = 0;
        dep_0 = [nA0, 0, T0];

        % set the stopping criterion for stage 1 of the protocol
        stopVar = 0;
        stopVal = t1;

        % solve the design equations for stage 1 of the protocol
        odes_are_stiff = false;
        [t, dep, success, message] = solve_ivodes(ind_0, dep_0, stopVar...
            , stopVal, @bstr_derivatives, odes_are_stiff);
        nA = dep(:,1);
        nZ = dep(:,2);
        T = dep(:,3);
        
        % check for solver issues
        if ~success
            disp('')
            disp(["BSTR function stage 1 solver issue: ",message])
            disp('')
        end
        
        % set the protocol stage to 2
        g_stage = 2;

        % set the initial values for stage 2 of the protocal
        ind_0 = t(end);
        dep_0 = [nA(end); nZ(end); T(end)];

        % set the stopping criterion for stage 2 of the protocol
        stopVar = 3;
        stopVal = Tf;

        % solve the design equations for stage 2 of the protocol
        odes_are_stiff = false;
        [t2, dep2, success, message] = solve_ivodes(ind_0, dep_0, stopVar...
            , stopVal, @bstr_derivatives, odes_are_stiff);
        nA2 = dep2(:,1);
        nZ2 = dep2(:,2);
        T2 = dep2(:,3);
        
        % check for solver issues
        if ~success
            disp('')
            disp(["BSTR function stage 2 solver issue: ", message])
            disp('')
        end

        % combine the results from stage 1 and stage 2
        t = [t; t2];
        nA = [nA; nA2];
        nZ = [nZ; nZ2];
        T = [T; T2];
    end

    % BSTR derivatives function
    function ddt = bstr_derivatives(~, dep)
        % extract the dependent variables
        nA = dep(1);
        T = dep(3);

        % calculate the additional unknowns
        k1 = k01 * exp(-E1 / (R * T));
        CA = nA / V;
        r1 = k1 * CA;
        if g_stage == 1
            U = U1;
            A = A1;
            Tex = Tex_1;
        else
            U = U2;
            A = A2;
            Tex = Tex_2;
        end
        Qdot = U * A * (Tex - T);

        % evaluate the derivatives
        dnA_dt = -r1 * V;
        dnZ_dt = r1*V;
        dT_dt = (Qdot - r1 * V * dH1) / (Cp * V);

        % return the derivatives as a column vector
        ddt = [dnA_dt; dnZ_dt; dT_dt];
    end

    % deliverables function
    function deliverables()
        % solve the BSTR design equations
        [t, nA, nZ, T] = bstr_model_variables();

        % calculate the net rate of production of Z
        rNet_Z = nZ(end) / (t(end) + t_turn);

        % tabulate, show, and save the results
        item = ["Heating Time"; "Process Time"; "Conversion"; "Net Rate"];
        value = [t1; t(end); 100 * (nA0 - nA(end)) / nA0; rNet_Z];
        units = ["h"; "h"; "%"; "mol Z /h"];
        resultsTable = table(item, value, units);
        disp('')
        disp(resultsTable)
        disp('')
        writetable(resultsTable, "results.csv")

        % plot the conversion and temperature profiles
        conv = 100 * (nA0 - nA) / nA0;
        figure
        plot(t, conv, 'k', 'LineWidth', 2)
        xline(t1, 'r:')
        xlabel('Time (h)','FontSize', 14)
        ylabel('Conversion (%)','FontSize', 14)
        set(gca, 'FontSize', 14)
        saveas(gcf,'conversion_profile.png')

        figure
        plot(t, T - 273.15, 'k', 'LineWidth', 2)
        xline(t1, 'r:')
        xlabel('Time (h)','FontSize', 14)
        ylabel('Temperature (°C)','FontSize', 14)
        set(gca, 'FontSize', 14)
        saveas(gcf,'temperature_profile.png')
    end

    % perform the calculations
    deliverables();
end
    