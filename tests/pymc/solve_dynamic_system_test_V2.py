# %%
# Dynamic Flux Balance Analysis (dFBA) in COBRApy
# The following notebook shows a simple, but slow example of implementing dFBA using COBRApy and scipy.integrate.solve_ivp.
# This notebook shows a static optimization approach (SOA) implementation and should not be considered production ready.
# The model considers only basic Michaelis-Menten limited growth on glucose.
import arviz as az
import numpy as np
import pandas as pd
import pymc as pm
import pytensor
import pytensor.tensor as pt
from tqdm import tqdm

from numba import njit
from pytensor.compile.ops import as_op
from scipy.integrate import solve_ivp

import matplotlib.pyplot as plt
# %matplotlib inline

# %%
# %load_ext watermark
az.style.use("arviz-darkgrid")
rng = np.random.default_rng(1234)
plt.ioff()  # turn off interactive plotting
output_messages = []  # store printed output

# %%
import os
import pickle
import json
from datetime import datetime

def save_dFBA_results(ts, y0, parameters, sol, data, trace, summary=None):
    """
    Save important data from dFBA simulation and MCMC sampling to files.
    
    Parameters:
    -----------
    ts : array
        Time points for evaluation
    y0 : list
        Initial conditions [biomass0, glucose0]
    parameters : list
        Model parameters [V, Km]
    sol : OdeResult
        Solution from solve_ivp
    data : DataFrame
        DataFrame with simulation results
    trace : InferenceData or MultiTrace
        MCMC sampling results
    summary : DataFrame, optional
        Summary of MCMC trace
    """
    # Create a timestamp for the output directory
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    output_dir = f"dFBA_results_{timestamp}"
    
    # Create the output directory if it doesn't exist
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    print(f"Saving results to directory: {output_dir}")
    
    # 1. Save simulation parameters
    params_dict = {
        # "timepoints": ts.tolist(),
        "initial_conditions": {
            "biomass": y0[0],
            "glucose": y0[1]
        },
        "parameters": {
            "V": parameters[0],
            "Km": parameters[1]
        }
    }
    
    with open(os.path.join(output_dir, "simulation_parameters.json"), "w") as f:
        json.dump(params_dict, f, indent=4)
    
    # 2. Save the ODE solution
    sol_dict = {
        "t": sol.t.tolist(),
        "y": [sol.y[0].tolist(), sol.y[1].tolist()],
        "status": sol.status,
        "success": sol.success,
        "message": sol.message
    }

    # Handle t_events which might be a list of arrays
    if hasattr(sol, 't_events') and sol.t_events is not None:
        sol_dict["t_events"] = [event.tolist() if isinstance(event, np.ndarray) and len(event) > 0 
                               else [] for event in sol.t_events]
    
    # Handle y_events which might be a list of arrays
    if hasattr(sol, 'y_events') and sol.y_events is not None:
        sol_dict["y_events"] = [event.tolist() if isinstance(event, np.ndarray) and len(event) > 0 
                               else [] for event in sol.y_events]

    with open(os.path.join(output_dir, "ode_solution.json"), "w") as f:
        json.dump(sol_dict, f, indent=4)
    
    # 3. Save the simulation data
    data.to_csv(os.path.join(output_dir, "simulation_data.csv"), index=False)
    
    # 4. Save the MCMC trace
    try:
        trace_file = os.path.join(output_dir, "mcmc_trace.pickle")
        with open(trace_file, "wb") as f:
            pickle.dump(trace, f)
        print(f"MCMC trace saved to {trace_file}")
    except Exception as e:
        print(f"Error saving trace with pickle: {e}")
        try:
            # Try to save using ArviZ if available
            import arviz as az
            az.to_netcdf(trace, os.path.join(output_dir, "mcmc_trace.nc"))
            print(f"MCMC trace saved using ArviZ")
        except Exception as e2:
            print(f"Error saving trace with ArviZ: {e2}")
            print(f"Warning: Could not save the MCMC trace")
    
    # 5. Save the MCMC summary
    if summary is not None:
        try:
            summary.to_csv(os.path.join(output_dir, "mcmc_summary.csv"))
            print(f"MCMC summary saved to {output_dir}/mcmc_summary.csv")
        except Exception as e:
            print(f"Error saving summary: {e}")
        
    print(f"Results saved successfully to {output_dir}")
    return output_dir

# %%
# Create or load a cobrapy model. Here, we use the ‘textbook’ e-coli core model.
import cobra
from cobra.io import load_model
cobra_model = load_model('textbook')

# %%
# Set up the dynamic system
# Dynamic flux balance analysis couples a dynamic system in external cellular concentrations to a pseudo-steady state metabolic model.
# In this notebook, we define the function add_dynamic_bounds(model, y) to convert the external metabolite concentrations into bounds on the boundary fluxes in the metabolic model.
def add_dynamic_bounds(cobra_model, y, parameters):
    """Use external concentrations to bound the uptake flux of glucose."""
    biomass, glucose = y  # expand the boundary species
    # Assumption: Michaelis-Menten kinetics for glucose uptake
    # V = maximum reaction rate
    # Km = Michaelis constant, which represents the substrate concentration at which the reaction rate is half of V.
    V, Km = parameters
    
    glucose_max_import = -V * glucose / (Km + glucose)
    cobra_model.reactions.EX_glc__D_e.lower_bound = glucose_max_import

    # Pre-optimize other constraints for faster solving
    cobra_model.solver.update()


def dynamic_system(t, y, parameters):
    """Calculate the time derivative of external species."""

    biomass, glucose = y  # expand the boundary species

    # Calculate the specific exchanges fluxes at the given external concentrations.
    with cobra_model:
        add_dynamic_bounds(cobra_model, y, parameters)

        cobra.util.add_lp_feasibility(cobra_model)
        feasibility = cobra.util.fix_objective_as_constraint(cobra_model)
        lex_constraints = cobra.util.add_lexicographic_constraints(
            cobra_model, ['Biomass_Ecoli_core', 'EX_glc__D_e'], ['max', 'max'])

    # Since the calculated fluxes are specific rates, we multiply them by the
    # biomass concentration to get the bulk exchange rates.
    fluxes = lex_constraints.values
    fluxes *= biomass

    # This implementation is **not** efficient, so I display the current
    # simulation time using a progress bar.
    if dynamic_system.pbar is not None:
        dynamic_system.pbar.update(1)
        dynamic_system.pbar.set_description('t = {:.3f}'.format(t))

    return fluxes

dynamic_system.pbar = None


def infeasible_event(t, y, parameters):
    """
    Determine solution feasibility.

    Avoiding infeasible solutions is handled by solve_ivp's built-in event detection.
    This function re-solves the LP to determine whether or not the solution is feasible
    (and if not, how far it is from feasibility). When the sign of this function changes
    from -epsilon to positive, we know the solution is no longer feasible.

    """

    with cobra_model:

        add_dynamic_bounds(cobra_model, y, parameters)

        cobra.util.add_lp_feasibility(cobra_model)
        feasibility = cobra.util.fix_objective_as_constraint(cobra_model)

    return feasibility - infeasible_event.epsilon

infeasible_event.epsilon = 1E-6
infeasible_event.direction = 1
infeasible_event.terminal = True

# %%
# Set up the simulation
def solve_dynamic_system(parameters):

    # with tqdm() as pbar:
        # dynamic_system.pbar = pbar

        V, Km = parameters

        sol = solve_ivp(
            fun=dynamic_system,
            events=[infeasible_event],
            t_span=(ts.min(), ts.max()),
            y0=y0,
            t_eval=ts,
            rtol=1e-6,
            atol=1e-8,
            method='BDF',
            args=(parameters,)
        )

        return sol

# %%
# Gradient-Free Sampler Options
# Like other Numpy or Scipy-based functions, the scipy.integrate.solve_ivp function cannot be used directly in a PyMC model
# because PyMC needs to know the variable input and output types to compile.
# Therefore, we use a Pytensor wrapper to give the variable types to PyMC.
# Then the function can be used in PyMC in conjunction with gradient-free samplers.

# %%
# Convert Python Function to a Pytensor Operator using @as_op decorator
# We tell PyMC the input variable types and the output variable types using the @as_op decorator.
# solve_ivp returns Numpy arrays, but we tell PyMC that they are Pytensor double float tensors for this purpose.

# Remember, when using @as_op with PyMC, the function becomes a "black box" to PyMC, meaning it can't compute gradients through it.
# This limits you to using gradient-free samplers.

# decorator with input and output types a Pytensor double float tensors
@as_op(itypes=[pt.dvector, pt.dvector], otypes=[pt.dmatrix])
def pytensor_model(parameters, y0):
    # with tqdm() as pbar:
        # dynamic_system.pbar = pbar

        # Extract parameters
        V, Km = parameters

        # Solve ODE system
        sol = solve_dynamic_system(parameters)

        # Important: Ensure we return the same dimensions as our data
        # Create a result array of the same shape as the data
        result = np.zeros_like(data[["biomass", "glucose"]].values)

        # Fill in the values we have (might be fewer than expected due to early termination)
        rows_to_fill = min(len(sol.t), result.shape[0])
        result[:rows_to_fill, 0] = sol.y[0, :rows_to_fill]
        result[:rows_to_fill, 1] = sol.y[1, :rows_to_fill]

        return result

# %%
# PyMC Model
# Now, we can specify the PyMC model using the ODE solver.
# For priors, we will use the manually set values: parameters = [10, 5] for V and Km,
# and y0 = [0.1, 10] for initial biomass and glucose concentrations.
# These are empirically derived weakly informative priors that are constrained to be positive.
# We will use a normal likelihood on untransformed data (i.e., not log transformed) to best fit the peaks of the data.

def pm_model(data, parameters, y0):

    with pm.Model() as model:
        # Priors for kinetic parameters and initial conditions
        V = pm.TruncatedNormal("V", mu=parameters[0], sigma=1.0, lower=0, initval=parameters[0])
        Km = pm.TruncatedNormal("Km", mu=parameters[1], sigma=0.5, lower=0, initval=parameters[1])
        biomassT0 = pm.TruncatedNormal("biomassT0", mu=y0[0], sigma=1, lower=0, initval=y0[0])
        glucoseT0 = pm.TruncatedNormal("glucoseT0", mu=y0[1], sigma=1, lower=0, initval=y0[1])
        sigma_biomass = pm.HalfNormal("sigma_biomass", 0.1)
        sigma_glucose = pm.HalfNormal("sigma_glucose", 1.0)

        # Ode solution function
        ode_solution = pytensor_model(
            pm.math.stack([V, Km]),
            pm.math.stack([biomassT0, glucoseT0])
        )

        # Extract the predicted biomass and glucose values
        biomass_pred = ode_solution[:, 0]
        glucose_pred = ode_solution[:, 1]

        # Likelihood
        pm.Normal("biomass_obs", mu=biomass_pred, sigma=sigma_biomass, observed=data["biomass"].values)
        pm.Normal("glucose_obs", mu=glucose_pred, sigma=sigma_glucose, observed=data["glucose"].values)

    return model

# %%
# pm.model_to_graphviz(model=model)

# %%
# Gradient-Free Sampler Options
# Having good gradient free samplers can open up the models that can be fit within PyMC. There are five options for gradient-free samplers in PyMC that are applicable to this problem (?):
# Slice - the default gradient-free sampler
# DEMetropolisZ - a differential evolution Metropolis sampler that uses the past to inform sampling jumps
# DEMetropolis - a differential evolution Metropolis sampler
# Metropolis - the vanilla Metropolis sampler
# SMC - Sequential Monte Carlo

# A few notes on running these inferences.
# For each sampler, the number of tuning steps and draws have been reduced to run the inference in a reasonable amount of time (on the order of minutes).
# This is not a sufficient number of draws to get a good inferences, in some cases, but it works for demonstration purposes.
# In addition, multicore processing was not working for the Pytensor op function on all machines, so inference is performed on one core.

def pm_sample(model, tune=1000, draws=1000, sampler="DEMetropolisZ"):
    with model:
        trace = pm.sample(
             tune=tune,
             draws=draws,
             step=pm.DEMetropolisZ(),
        )

    return trace

# %%
# Run the dynamic FBA simulation
ts = np.linspace(0, 15, 100)  # Desired integration resolution and interval
y0 = [0.1, 10]  # Initial conditions for biomass and glucose
parameters = [10, 5]  # Parameters for the Michaelis-Menten kinetics

sol = solve_dynamic_system(parameters)

# Because the culture runs out of glucose, the simulation terminates early.
# The exact time of this ‘cell death’ is recorded in sol.t_events.

# Plot timelines of biomass and glucose
# ax = plt.subplot(111)
# ax.plot(sol.t, sol.y.T[:, 0], color='b')
# ax2 = plt.twinx(ax)
# ax2.plot(sol.t, sol.y.T[:, 1], color='r')

# ax.set_ylabel('Biomass', color='b')
# ax2.set_ylabel('Glucose', color='r')
# plt.show()

# Store the simulation results in a pandas DataFrame
data = pd.DataFrame(dict(time=sol.t, biomass=sol.y[0], glucose=sol.y[1]))

# %%
# Run the full analysis
def run_full_analysis():
    # Get the existing global variables
    global ts, y0, parameters, sol, data

    # Run the PyMC model
    model = pm_model(data, parameters, y0)

    # Run with minimal settings for testing
    print("Running sampling with minimal settings for testing...")
    trace = pm_sample(model, tune=10, draws=10, sampler="DEMetropolisZ")
    
    # Print summary
    summary = az.summary(trace)
    print("\nSampling summary:")
    print(summary)

    # Save all the important data
    try:
        save_dFBA_results(ts, y0, parameters, sol, data, trace, summary)
        print("All data saved successfully")
    except Exception as e:
        print(f"Error during data saving: {e}")
    
    return trace, data, ts

# %%
# Inference!
if __name__ == '__main__':
    trace, data, ts = run_full_analysis()
