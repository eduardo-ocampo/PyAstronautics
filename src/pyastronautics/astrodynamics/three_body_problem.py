""" 
   Copyright 2024 Eduardo Ocampo
   https://github.com/eduardo-ocampo/PyAstronautics
"""

import math
import pickle
import numpy as np
from typing import Union
from numpy.linalg import norm
from scipy.integrate import solve_ivp

class CR3BP(object):
    """
    Non-Dimensional Circular Restricted Three-Body Problem as defined in the
    jacobi coordinate frame and shifted into the rotating frame.

    This Python class is derived from and setup as two_body_problem:TwoBodyModel

    Attributes
    ----------
    mu : float
        Mass ratio of the primary bodies, set to Earth-Moon System by default.
    position : list[float]
        The 3D position vector of the object.
    velocity : list[float] 
        The 3D velocity vector of the object.
    position_norm : float
        The magnitude of the position vector.
    velocity_norm : float
        The magnitude of the velocity vector.
    initial_state_vector : list[float]
        The concatenated initial state vector combining
        position and velocity, used for numerical analysis.
    abs_tol : float
        Absolute tolerance value for numerical analysis, default is 1e-10.
    rel_tol : float
        Relative tolerance value for numerical analysis, default is 1e-10.
    num_sol_pickle_file : str
        The filename for the pickle file to store numerical solution data.
   
    """
    def __init__(self, position: list[float], velocity: list[float]):
        """
        Initialize the Non-Dimensional CR3BP instance with initial
        position and velocity vectors.

        Parameters
        ----------
        position : list[float]
            The 3D position vector
        velocity : list[float]
            The 3D velocity vector

        Raises
        ------
        TypeError
            If position or velocity is not a list.
        """
        # Set Non-Dimensional CR3BP Primary Bodies Mass Ratio
        # ---------------------------------------   
        # Default Set to Earth-Moon System
        self.mu = 0.012150515586657583

        # Set Position and Velocity Vector
        # ---------------------------------------   
        if not isinstance(position, list):
            raise TypeError("position must be a list.")
        if not isinstance(velocity, list):
            raise TypeError("velocity must be a list.")
        
        self.position = position
        self.velocity = velocity

        self.position_norm = norm(position)
        self.velocity_norm = norm(velocity)

        # Numerical Analysis Setup
        # ---------------------------------------
        self.initial_state_vector = position + velocity
        # Default tolerance values
        self.abs_tol = 1e-10
        self.rel_tol = 1e-10

        # Set File Names
        # ---------------------------------------   
        self.num_sol_pickle_file = "cr3bp_solution.pickle"  

    def non_dim_differential_equations(self, t: float, state: Union[list, np.ndarray]) -> np.ndarray:
        """
        Define the non-dimensional differential equations for the Circular Restricted Three-Body
        Problem using their Equations of Motions.

        This method computes the derivatives of the state vector `state`, which includes both position and 
        velocity components. The equations describe the motion of a third body under its mutual 
        gravitational influence of two bodies. Assuming the third body has zero mass. 

        The vector `state` is expected to be structured as follows:
        [x, y, z, vx, vy, vz], where:
            - x, y, z: position coordinates of the body
            - vx, vy, vz: velocity components of the body

        Parameters
        ----------
        t : float
            The current time in the simulation.

        state : Union[list, np.ndarray]
            The state vector containing the current position and velocity of the body. 
            This can be a list or a NumPy array.

        Returns
        -------
        np.ndarray
            The derivatives of the state vector (ndarray), consisting of position and velocity derivatives.
        """

        if isinstance(state, list):
            # Convert list to NumPy array
            state = np.array(state)
        elif not isinstance(state, np.ndarray):
            raise ValueError("state must be a list or a NumPy array.")
            
        x,y,z, vx,vy,vz = state
        
        # Compute Differential Equation Constants: Position to Primary Bodies
        r1 = math.sqrt((x+self.mu)**2 + y**2 + z**2)
        r2 = math.sqrt((x-1+self.mu)**2 + y**2 + z**2)

        # Differential Equations: ddot is a second derivative
        x_ddot =  2*vy + x - (1-self.mu)*(x+self.mu)/r1**3 - self.mu*(self.mu+x-1)/r2**3
        y_ddot = -2*vx + y - y*(1-self.mu)/r1**3 - self.mu*y/r2**3
        z_ddot =  -z*(1-self.mu)/r1**3 - self.mu*z/r2**3

        # Return d/dt vector of
        # [x, y, z, vx, vy, vx]
        return np.concatenate(([vx,vy,vz],[x_ddot,y_ddot,z_ddot]))

    def solve_non_dim_trajectory(self, save_analysis:bool = False) -> None:
        """
        Solve the trajectory of the Non-Dimensional Circular Restricted Three-Body Problem using the
        initial value problem (IVP).

        This method uses the `scipy.integrate.solve_ivp()` function to numerically integrate the differential equations 
        governing the motion of the bodies, given the initial conditions specified in `self.initial_state_vector`.

        Before calling this method, ensure that `self.time` is defined as a sequence of non-dimensional time points
        over which the simulation will be evaluated. If `self.time` is not defined, a ValueError will be raised.

        The initial value problem must be set up in the following format:
        [x, y, z, vx, vy, vz], where:
            - x, y, z: initial position coordinates
            - vx, vy, vz: initial velocity components

        The results of the integration are stored in `self.num_sol`, and the position and velocity 
        results are extracted into `self.numerical_position` and `self.numerical_velocity`, respectively.

        The final state which corresponds to `self.time[-1]` is stored as `self.final_state`.

        Parameters
        ----------
        save_analysis : bool, optional
            If True, analysis results will be saved for later use. Defaults to False.

        Raises
        ------
        ValueError
            If `self.time` is not defined or is empty.

        Returns
        -------
        None
        """

        # Check if self.time is defined
        if not hasattr(self, 'time'):
            raise ValueError("Attribute 'time' must be defined before calling solve_trajectory.")

        ivp = self.initial_state_vector

        self.num_sol = solve_ivp(self.non_dim_differential_equations,
                                [self.time[0],self.time[-1]],
                                 ivp,
                                 t_eval=self.time,
                                 rtol=self.rel_tol,
                                 atol=self.abs_tol)

        # Check if solver reached interval end or a termination event occurred 
        if not self.num_sol.success:
            print(f"Solver termination status: {self.num_sol.status}")
        else:
            print(f"Solver Success: {self.num_sol.success}")

        # Extract Position and Velocity Results
        self.numerical_position = self.num_sol.y[:3,:].T
        self.numerical_velocity = self.num_sol.y[3:,:].T

        self.final_state = self.num_sol.y[:,-1].T

        # Allow user to save numerical analysis
        if save_analysis:
            with open(self.num_sol_pickle_file, 'wb') as handle:
                pickle.dump(self, handle)
