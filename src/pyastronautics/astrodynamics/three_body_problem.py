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
from scipy.optimize import newton

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

    @staticmethod
    def get_jacobi_velocity(x: float, y: float, jacobi: float, mass_ratio: float) -> float:
        """
        Compute the maximum velocity magnitude at a given position (x, y) in space 
        for a specific Jacobi constant. This is the velocity corresponding to the 
        given Jacobi constant at the point (x, y) in a planar restricted three-body problem 
        or similar system.

        Parameters
        ----------
        x : float
            The x-coordinate in the planar non-dimensional system.
        y : float
            The y-coordinate in the planar non-dimensional system.
        jacobi : float
            The Jacobi constant at the point (x, y).
        mass_ratio : float
            The mass ratio between the two bodies in the system.

        Returns
        -------
        float
            The maximum velocity magnitude (in non-dimensional units) at the specified 
            position corresponding to the given Jacobi constant.

        Examples
        --------
        velocity = get_jacobi_velocity(0.5, 0.5, 2.5, 0.1)

        """
        
        r1 = np.sqrt((x + mass_ratio)**2 + y**2)
        r2 = np.sqrt((x - 1 + mass_ratio)**2 + y**2)

        vel_mag = np.sqrt(((x**2 + y**2) + 2*(1 - mass_ratio)/r1 + 2*mass_ratio/r2) - jacobi)

        return vel_mag
    
    @staticmethod
    def calculate_jacobi(x: float, y: float, vx: float, vy: float, mass_ratio: float) -> float:
        """
        Calculate the Jacobi constant at a given position (x, y) and velocity (vx, vy)
        in a planar restricted three-body problem or similar system, based on the given mass ratio.

        Parameters
        ----------
        x : float
            The x-coordinate in the planar non-dimensional system.
        y : float
            The y-coordinate in the planar non-dimensional system.
        vx : float
            The x-component of velocity.
        vy : float
            The y-component of velocity.
        mass_ratio : float
            The mass ratio between the two bodies in the system.

        Returns
        -------
        float
            The calculated Jacobi constant at the specified position and velocity.

        Notes
        -----
        This function assumes a planar system where only the x and y components of position and velocity
        are considered (no motion in the z-direction).

        Examples
        --------
        jacobi = calculate_jacobi(0.5, 0.5, 0.1, 0.1, 0.1)
        """
        
        # Calculate distances to the two bodies
        r1 = np.sqrt((x + mass_ratio)**2 + y**2)
        r2 = np.sqrt((x - 1 + mass_ratio)**2 + y**2)
        
        # Calculate the Jacobi constant
        jacobi = (x**2 + y**2) + 2 * (1 - mass_ratio) / r1 + 2 * mass_ratio / r2 - (vx**2 + vy**2)
        
        return jacobi

    @staticmethod
    def forbidden_region(jacobi_max: float, mass_ratio: float, 
                         x_range: list = [-1.5, 1.5],
                         y_range: list = [-1.5, 1.5],
                         linspace_num: int = 200) -> tuple:
        """
        Compute the forbidden region in the x-y plane based on the Jacobi constant and 
        mass ratio. The forbidden region corresponds to the area where the Jacobi constant 
        exceeds the provided `jacobi_max`. Knowing the velocity term in the Jacobi equation
        is certain to always be positive this method relies only on position. The function
        assumes a non-dimensional system and generates a meshgrid in the given x-y range
        to identify the forbidden region.

        Parameters
        ----------
        jacobi_max : float
            The Jacobi constant threshold. Points with Jacobi constant greater than this 
            value are considered in the forbidden region.
        mass_ratio : float
            The mass ratio between the two bodies in the system (e.g., in the restricted 
            three-body problem).
        x_range : list of float, optional
            The range of x-values for the meshgrid (default is [-1.5, 1.5]).
        y_range : list of float, optional
            The range of y-values for the meshgrid (default is [-1.5, 1.5]).
        linspace_num : int, optional
            The number of points to use for the linspace along each axis (default is 200).

        Returns
        -------
        tuple
            A tuple containing:
                - x_vals : numpy.ndarray
                    The x-values of the meshgrid.
                - y_vals : numpy.ndarray
                    The y-values of the meshgrid.
                - jc : numpy.ndarray
                    The Jacobi constant values on the meshgrid, with values greater than 
                    `jacobi_max` set to NaN to indicate the forbidden region.

        Examples
        --------
        x_vals, y_vals, jc = forbidden_region(3.0, 0.1)

        """
        
        x_min, x_max = x_range
        y_min, y_max = y_range
        
        x_vals = np.linspace(x_min, x_max, linspace_num)
        y_vals = np.linspace(y_min, y_max, linspace_num)

        [X, Y] = np.meshgrid(x_vals, y_vals)

        r1 = np.sqrt((X + mass_ratio)**2 + Y**2)
        r2 = np.sqrt((X - 1 + mass_ratio)**2 + Y**2)

        jc = X**2 + Y**2 + 2*(1 - mass_ratio)/r1 + 2*mass_ratio/r2
        jc[jc > jacobi_max] = np.nan

        return x_vals, y_vals, jc


class planar_lagrange_points(object):
    """
    A class to calculate the Planar Lagrange points (L1, L2, L3, L4, L5) in a Circular Restricted 
    Three-Body Problem, based on a given mass ratio between the two bodies.

    The Lagrange points are solutions to the equations of motion in the restricted 
    three-body problem, where the two massive bodies (e.g., Sun and Earth) are in orbit, 
    and the third object (e.g., satellite) can remain in a "stable" position relative to the two bodies.

    Attributes
    ----------
    mu : float
        The mass ratio between the two primary bodies in the system. This is the ratio of the mass 
        of the secondary body to the total mass of the two-body system.
    l1x : float
        The x-coordinate of the L1 Lagrange point.
    l2x : float
        The x-coordinate of the L2 Lagrange point.
    l3x : float
        The x-coordinate of the L3 Lagrange point.
    l1y, l2y, l3y : float
        The y-coordinates of the collinear Lagrange points (L1, L2, L3) are initialized to zero. 
        These points lie along the x-axis in a restricted 3-body problem.
    l4x, l4y, l4 : float
        The x and y coordinates, and the distance from the barycenter, of the L4 triangular 
        Lagrange point.
    l5x, l5y, l5 : float
        The x and y coordinates, and the distance from the barycenter, of the L5 triangular 
        Lagrange point (which is symmetric to L4).

    """
    
    def __init__(self, mass_ratio: float):
        """
        Initializes the Lagrange point object with the given mass ratio of the two primary bodies.
        
        Parameters
        ----------
        mass_ratio : float
            The mass ratio between the two primary bodies in the system. The mass ratio (mu) 
            is defined as the ratio of the secondary mass to the total system mass.
        """
        self.mu = mass_ratio

        self.l1y, self.l2y, self.l3y = 3*[0]

    def get_points(self):
        """
        Calculates both the collinear and triangular Lagrange points (L1, L2, L3, L4, L5).
        
        This method calls the `colinear_points()` and `triangular_points()` methods to 
        compute the positions of the Lagrange points in the system.
        """
        self.colinear_points()
        self.triangular_points()
    
    def colinear_points(self):
        """
        Calculates the collinear Lagrange points (L1, L2, L3) by approximating the positions 
        and then solving the root equation numerically using Newton's method.
        
        The collinear points lie along the line connecting the two primary bodies, 
        and their positions are determined by solving for the points where the net force 
        from both bodies is zero.
        
        Returns
        -------
        list of float
            The x-coordinates of the collinear Lagrange points [L1x, L2x, L3x].
        """
        # First, calculate the initial guesses for the collinear points
        x_guess = self.colinear_approximation()
        
        # Solve the root equation for each collinear point using Newton's method
        lagrange_xsol = newton(func=self.root_equation, x0=x_guess, tol=1e-16, rtol=1e-16) 
        
        # Assign the solutions to the collinear Lagrange points
        self.l1x, self.l2x, self.l3x = lagrange_xsol

        return [self.l1x, self.l2x, self.l3x]

    def colinear_approximation(self) -> list:
        """
        Provides an approximation for the positions of the collinear Lagrange points (L1, L2, L3).
        
        These approximations are based on an analytical formula that provides a good starting 
        guess for solving the root equation.
        
        Returns
        -------
        list of float
            Initial guesses for the positions of the collinear Lagrange points [L1, L2, L3].
        """
        # Approximations based on the mass ratio and system geometry
        alpha_guess = ((self.mu / 3) * (1 - self.mu))**(1/3)
        l1_approx = 1 - self.mu - alpha_guess

        beta_guess = ((self.mu / 3) * (1 - self.mu))**(1/3)
        l2_approx = beta_guess - self.mu + 1

        gamma_guess = (-7 * self.mu / 12) + 1
        l3_approx = -self.mu - gamma_guess

        return [l1_approx, l2_approx, l3_approx]

    def root_equation(self, x: float) -> float:
        """
        Defines the equation to be solved for the collinear Lagrange points (L1, L2, L3).
        
        The equation is derived from the balance of forces at the collinear points using the force potential equation.
        
        Parameters
        ----------
        x : float
            The x-coordinate to be evaluated in the root equation.
        
        Returns
        -------
        float
            The value of the function at the point `x`.
        """
        # Lagrange x-coordinate equation based on force potential
        func = x - self.mu*(self.mu + x - 1) / np.abs(self.mu + x - 1)**3 + (self.mu - 1)*(self.mu + x) / np.abs(self.mu + x)**3
        return func

    def triangular_points(self):
        """
        Calculates the triangular Lagrange points (L4 and L5) using the mass ratio of the system.
        
        L4 and L5 form an equilateral triangle with the two primary bodies. The calculation 
        assumes that the two massive bodies are located along the x-axis, and L4 and L5 
        are positioned symmetrically with respect to the line joining the two bodies.
        """

        self.l4x  = (1/2 - self.mu)
        self.l4y = (math.sqrt(3)/2)
        # Distance from barycenter
        self.l4 = math.sqrt(self.l4x**2 + self.l4y**2)
        
        # L5 is symmetric to L4
        self.l5x, self.l5y, self.l5 = self.l4x, -1*self.l4y, self.l4 