""" 
   Copyright 2024 Eduardo Ocampo
   https://github.com/eduardo-ocampo/PyAstronautics
"""

import pickle
import numpy as np
from typing import Union
from numpy.linalg import norm
from scipy.integrate import solve_ivp

from .base_model import TwoBodyOrbitalModel

class TwoBodyModel(TwoBodyOrbitalModel):
    """
    Represents a two-body system with position and velocity vectors.

    Attributes
    ----------
    orbit_elements : OrbitElements 
        An instance that holds the orbital parameters.
    mu : float
        The gravitational parameter (standard gravitational constant) of the central body,
        set to Earth's gravitational constant (3.986004418E+05 km^3/sec^2).
    position : list[float]
        The 3D position vector of the object in kilometers (km).
    velocity : list[float] 
        The 3D velocity vector of the object in kilometers per second (km/s).
    position_norm : float
        The magnitude of the position vector in kilometers (km).
    velocity_norm : float
        The magnitude of the velocity vector in kilometers per second (km/s).
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
        Initialize the TwoBodyModel instance with initial position and velocity vectors.

        Parameters
        ----------
        position : list[float]
            The 3D position vector in km
        velocity : list[float]
            The 3D velocity vector in km/sec

        Raises
        ------
        TypeError
            If position or velocity is not a list.
        
        """

        # Set Gravitational Constant
        # ---------------------------------------   
        # Default Set to Earth
        self.mu = 3.986004418E+05 # km^3/sec^2

        # Set Position and Velocity Vector
        # ---------------------------------------   
        # Assumed input is in (km, sec)    
        if not isinstance(position, list):
            raise TypeError("position must be a list.")
        if not isinstance(velocity, list):
            raise TypeError("velocity must be a list.")
        
        self.position = position
        self.velocity = velocity

        self.position_norm = norm(position)
        self.velocity_norm = norm(velocity)

        # Initialize Orbit Element Object
        # ---------------------------------------
        # Call the parent class's __init__ method
        super().__init__(position, velocity)

        # Numerical Analysis Setup
        # ---------------------------------------
        self.initial_state_vector = position + velocity
        # Default tolerance values
        self.abs_tol = 1e-10
        self.rel_tol = 1e-10

        # Set File Names
        # ---------------------------------------   
        self.num_sol_pickle_file = "twoBody_trajectory.pickle"        

    def differential_equations(self, t: float, state: Union[list, np.ndarray]) -> np.ndarray:
        """
        Define the differential equations for the Two-Body Problem using their
        Equations of Motions.

        This method computes the derivatives of the state vector `state`, which includes both position and 
        velocity components. The equations describe the motion of two bodies under their mutual 
        gravitational influence.

        The state vector `state` is expected to be structured as follows:
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

        pos = state[0:3]
        vel = state[3:]

        # Compute Differential Equation Constants
        constant = -self.mu/(norm(pos)**3)

        # Differential Equations
        accel = np.dot(constant,pos)

        # Return d/dt vector of
        # [x, y, z, vx, vy, vx]
        return np.concatenate((vel,accel))

    def solve_trajectory(self, save_analysis:bool = False) -> None:
        """
        Solve the trajectory of a Two-Body system using the initial value problem (IVP).

        This method uses the `scipy.integrate.solve_ivp()` function to numerically integrate the differential equations 
        governing the motion of the bodies, given the initial conditions specified in `self.initial_state_vector`.

        Before calling this method, ensure that `self.time` is defined as a sequence of time points in seconds
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

        self.num_sol = solve_ivp(self.differential_equations,
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
        self.numerical_velocity = self.num_sol.y[3::,:].T

        self.final_state = self.num_sol.y[:,-1].T

        # Allow user to save numerical analysis
        if save_analysis:
            with open(self.num_sol_pickle_file, 'wb') as handle:
                pickle.dump(self, handle)

