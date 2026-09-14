function activity_13()
% Calculations for the Class 13 Learning Activity in REB, The Course
    % global constants available to all function
    % given
    V = 4.430e3; % cm3
    CS_in = 0.04; % g/cm3
    CX_in = 0; % g/cm3
    Vmax = 0.014; % /min
    Km = 0.001; % g/cm3
    Vdot_in = 61; % cm3/min
    CS_0 = 2.97e-2; % g/cm3
    CX_0 = 4.68e-3; % g/cm3
    % calculated
    Vdot = Vdot_in;

    % CSTR model function
    function [t, mS, mX] = cstr_model_variables()
        % set the initial values
        ind_0 = 0;
        dep_0 = [CS_0, CX_0] *Vdot;

        % set the stopping criterion
        stopVar = 0;
        stopVal = 40000;

        % solve the cstr design equations
        odes_are_stiff = false;
        [ind, dep, success, message] = solve_ivodes(ind_0, dep_0, stopVar...
            , stopVal, @cstr_derivatives, odes_are_stiff);
        
        % check for solver issues
        if ~success
            disp('')
            disp(["CSTR model function issue: ",message])
            disp('')
        end
        
        % extract the cstr model variables
        t = ind;
        mS = dep(:,1);
        mX = dep(:,2);
    end

    % CSTR derivatives function
    function ddt = cstr_derivatives(~, dep)
        % extract the individual dependent variables
        mS = dep(1);
        mX = dep(2);

        % calculate the additional unknowns
        CS = mS/Vdot;
        CX = mX/Vdot;
        r1 = Vmax*CS*CX/(Km + CS);
        mS_in = CS_in*Vdot_in;
        mX_in = CX_in*Vdot_in;

        % evaluate the derivatives
        dmSdt = Vdot_in/V*(mS_in - mS -2.2*V*r1);
        dmXdt = Vdot/V*(mX_in - mX + V*r1);

        % return the derivatives as a column vector
        ddt = [dmSdt; dmXdt];
    end

    % deliverables function
    function deliverables()
        % solve the cstr design equations
        [t, ~, mX] = cstr_model_variables();

        % generate, show and save a plot of CX vs time
        CX = mX/Vdot;
        figure
        plot(t, CX, 'LineWidth', 2)
        xlabel("Time (min)", 'FontSize', 14)
        ylabel("Outlet Cell Mass Concentration (g/cc)", 'FontSize', 14)
        set(gca, 'FontSize', 14);
        saveas(gcf,"CX_vs_t.png")
    end

    % execution command
    deliverables();
end