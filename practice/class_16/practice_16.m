function practice_16()
% Calculations for the Class 16 Practice Assignment from REB, The Course
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
    Tf = 35 + 273.15; % K
    t_turn = [0.25, 0.5, 0.75]; % h
    % known
    R = 1.987; % cal / mol / K
    % calculated
    nA0 = CA0 * V;

    % global variable indication the stage in the protocol
    g_stage = 1;

    % BSTR reactor function
    function [t, nA, nZ, T] = bstr_model_variables(t1)
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
            disp(["BSTR function stage 1 solver issue: ", message])
            disp('')
        end
        
        % set the protocol stage to 2
        g_stage = 2;

        % set the initial values for stage 2 of the protocal
        ind_0 = t(end);
        dep_0 = [nA(end), nZ(end), T(end)];

        % set the stopping criterion for stage 2 of the protocol
        stopVar = 3;
        stopVal = Tf;

        % solve the design equations for stage 2 of the protocol
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
        % allocate storage for the maximum net rate and corresponding conversion
        rNet_max = [nan, nan, nan];
        conv_opt = [nan, nan, nan];

        % allocate storage for the plot data
        conv = ones(100,3)*nan;
        rNet_Z = ones(100,3)*nan;

        % generate the plot data for the different values of t_turn
        for i = 1:3
            % set a range of heating times
            t1_range = linspace(0.05, 1.5, 100);

            % loop through the different values of t1
            for j=1:100
                % solve the BSTR design equations
                [t, nA, nZ, ~] = bstr_model_variables(t1_range(j));

                % calculate the net rate of production of Z
                rNet_Z(j,i) = nZ(end)/ (t(end) + t_turn(i));
            
                % calculate the conversion
                conv(j,i) = 100 * (nA0 - nA(end)) / nA0;
            end

            % save the results
            [rNet_max(i), iMax] = max(rNet_Z(:,i));
            conv_opt(i) = conv(iMax);
        end

        % generate the graph
        figure
        plot(conv(:,1), rNet_Z(:,1), 'b'...
            , conv(:,2), rNet_Z(:,2), 'r'...
            , conv(:,3), rNet_Z(:,3), 'g', 'LineWidth',2)
        set(gca, 'FontSize', 14);
        xlabel('Conversion (%)','FontSize', 14)
        ylabel('Net Rate (mol Z/h)','FontSize', 14)
        legend({"0.25 h", "0.5 h", "0.75 h"})
        saveas(gcf,"rate_vs_conversion.png")

        % tabulate show and save the maximum net rates
        turnaround_time = t_turn';
        max_net_rate = rNet_max';
        conversion = conv_opt';
        resultsTable = table(turnaround_time, max_net_rate, conversion);
        disp('')
        disp(resultsTable)
        disp('')
        writetable(resultsTable,'results.csv')
        end

    % perform the calculations
    deliverables();
end