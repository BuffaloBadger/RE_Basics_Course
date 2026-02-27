import numpy as np
import statsmodels.api as sm
from scipy.stats.distributions import t
from scipy.optimize import curve_fit
from scipy.optimize import root
from scipy.integrate import solve_ivp
from scipy.integrate import solve_bvp

def Arrhenius_parameters(k,T,R):
    # revised 12/11/25

    # define x and y in the linearized Arrhenius expression
    x = 1/T
    y = np.log(k)

    # add a constant to x because the model includes an intercept
    x = sm.add_constant(x)

    # fit a linear model to the data
    res = sm.OLS(y,x).fit()

    # extract the slope and intercept and their 95% confidence intervals
    beta = res.params
    beta_ci = res.conf_int(alpha=0.05)

    # extract the coefficient of determination
    r_squared = res.rsquared

    # calculate the Arrhenius parameters and their 95% confidence intervals
    k0 = np.exp(beta[0])
    k0_ci = np.exp(beta_ci[0,:])
    E = -R*beta[1]
    E_ci = np.zeros(2)
    E_ci[0] = -R*beta_ci[1,1]
    E_ci[1] = -R*beta_ci[1,0]

    # return the results
    return k0, k0_ci, E, E_ci, r_squared

def solve_ivodes(ind_0, dep_0, stop_var, stop_val, derivs_fcn, odes_are_stiff,
                 rel_tol = 1.0E-3, abs_tol = 1.0E-6):
    # revised 9/22/25
    
    # define an event for when the final value of a dependent variable is known
    def event(ind, dep):
        return dep[stop_var - 1] - stop_val
    event.terminal = True
    
    if stop_var == 0:
        # do not use the event
        ind = (ind_0, stop_val)
        if odes_are_stiff:
            soln = solve_ivp(derivs_fcn, ind, dep_0, method='LSODA'
                             , rtol = rel_tol, atol = abs_tol)
        else:
            soln = solve_ivp(derivs_fcn, ind, dep_0, method='RK45'
                             , rtol = rel_tol, atol = abs_tol)
        solved = soln.success
        message = soln.message
    else:
        # use the event
        count = 0
        solved = False
        ind_f = ind_0 + 1.0
        while (count < 10 and not solved):
            ind = (ind_0, ind_f)
            count +=1
            if odes_are_stiff:
                soln = solve_ivp(derivs_fcn, ind, dep_0, method='LSODA'
                                 , events=event, rtol = rel_tol, atol = abs_tol)
            else:
                soln = solve_ivp(derivs_fcn, ind, dep_0, method='RK45'
                                 , events=event, rtol = rel_tol, atol = abs_tol)
            if soln.t[-1] == ind_f: # ind_f was not large enough
                ind_f = (ind_0 + 1.0)*10**count
                message = 'The ivode stopping criterion was not reached.'
            elif soln.t[-1] < 0.1*ind_f: # ind_f was too large
                ind_f = (ind_0 + 1.0)/10**count
                message = 'The ivode integration results could be inaccurate.'
            else:
                success = soln.success
                message = soln.message
    return soln.t, soln.y, success, message

def fit_to_SR_data(beta_guess, x, y_meas, pred_resp_fcn, use_rel_error):
    # revised 9/22/25

    # set the weights
    if use_rel_error:
        weight = y_meas
    else:
        weight = np.ones(len(y_meas))
    
    # estimate the parameters
    beta, beta_cov, info, mesg, ier = curve_fit(pred_resp_fcn, x, y_meas,
            p0=beta_guess, sigma=weight, method = 'trf', full_output=True)

    # calculate r_squared
    y_pred = info['fvec'] + y_meas
    y_mean = np.mean(y_meas)
    ss_res = np.sum(np.square(y_meas - y_pred))
    ss_tot = np.sum(np.square(y_meas - y_mean))
    r_squared = 1 - ss_res/ss_tot

    # calculate 95% confidence interval
    beta_ci = np.zeros((len(beta),2))
    alpha = 0.05
    dof = max(0, len(y_meas) - len(beta))
    t_val = t.ppf(1.0 - alpha/2., dof)
    for i, p in enumerate(beta):
        beta_ci[i,0] = beta[i] - beta_cov[i,i]**0.5*t_val
        beta_ci[i,1] = beta[i] + beta_cov[i,i]**0.5*t_val
    
    # return the results
    return beta, beta_ci, r_squared

def solve_ates(residuals_fcn, guess):
    # revised 9/22/25
    
    # solve the ATEs
    soln = root(residuals_fcn,guess)

    # extract and return the results
    return soln.x, soln.success, soln.message

def solve_bvodes(xB1, xB2, dep_guess, derivs_fcn, resids_fcn):
    # revised 9/22/25
    
    # set up the initial mesh
    ind = np.linspace(xB1, xB2, 100)

    # solve the bvodes
    soln = solve_bvp(derivs_fcn, resids_fcn, ind, dep_guess)

    # return the results
    return soln.x, soln.y