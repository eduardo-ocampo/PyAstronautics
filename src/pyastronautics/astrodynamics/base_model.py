""" 
   Copyright 2024 Eduardo Ocampo
   https://github.com/eduardo-ocampo/PyAstronautics
"""

import math
import numpy as np
from typing import Union
from numpy.linalg import norm

class Parameter:
    """
    Represents a parameter with a value and its unit.

    Attributes
    ----------
    value : Union[list[float], np.ndarray, float, int]
        The parameter value, which can be a list, a NumPy array,
        a single float, or an integer.
    unit : str
        The unit of the parameter.
    description : str
        A short description fo the parameter. Default to None
    """
    def __init__(self, value: Union[list[float], np.ndarray, float, int],
                 unit: str, description: str = None):
        if not isinstance(value, (list, np.ndarray, float, int)):
            raise TypeError("value must be a list, NumPy array, float, or integer.")
        self.value = value
        self.unit = unit
        self.description = description

    def __repr__(self) -> str:
        """Return a string representation of the parameter."""
        return f"Parameter(value={self.value}, unit='{self.unit}', description='{self.description})"


class OrbitElements:
    """
    Stores various orbital parameters.

    Attributes
    ----------
    parameters : dict[str, Parameter]
        A dictionary storing parameters by their name (e.g., position, velocity).
    """
    def __init__(self):
        self.parameters: dict[str, Parameter] = {}

    def add_parameter(self, name: str, value: Union[list[float], np.ndarray, float, int],
                      unit: str, description: str = None):
        """
        Add a new parameter to the orbit elements repository.

        Parameters
        ----------
        name : str 
            The name of the parameter (e.g., "position").
        value : Union[list[float], np.ndarray, float, int]
            The value of the parameter.
        unit : str
            The unit of the parameter.
        description : str
            A short description of the parameter.
        """
        self.parameters[name] = Parameter(value, unit, description)

    def __getattr__(self, name: str) -> Parameter:
        """Allow access to parameters via dot notation.

        Parameters
        ----------
        name : str
            The name of the parameter to retrieve.

        Returns
        -------
        Parameter

        Raises:
            AttributeError: If the parameter does not exist.
        """
        if name in self.parameters:
            return self.parameters[name]
        raise AttributeError(f"'OrbitElements' object has no attribute '{name}'")

    def __repr__(self) -> str:
        """Return a string representation of the parameters."""
        return f"OrbitElements({', '.join(self.parameters.keys())})"

class TwoBodyOrbitalModel():

    def __init__(self, position: list[float], velocity: list[float]):
        """

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

        # Check Position and Velocity
        # ---------------------------------------   
        # Assumed input is in (km, sec)    
        if not isinstance(position, list):
            raise TypeError("position must be a list.")
        if not isinstance(velocity, list):
            raise TypeError("velocity must be a list.")

        # Initialize Orbit Element Dictionary
        # ---------------------------------------
        self.orbit_elements = OrbitElements()
        
        # Add position and velocity to orbit elements
        self.orbit_elements.add_parameter("position", position, "km")
        self.orbit_elements.add_parameter("velocity", velocity, "km/s")

    def calc_orbit_elements(self):
        """
        Using the initial position and velocity vectors, calculates and sets Conservation
        Parameters, Keplerian Orbit Elements, Initial Anomaly Parameters, and the
        Orbital Period for the Two-Body Problem. In both scalar and vector form when
        made possible.

        Below is a detailed list of parameters added to the existing OrbitElement Dictionary

        Conservation Parameters
        - Specific Energy (E)
        - Specific Angular Momentum (h)
        
        Keplerian Orbit Elements
        - Eccentricity (e)
        - Inclination (i)
        - Longitude of Ascending Node Ω (Omega)
        - Argument of Perigee ω (omega)
        - Semi-major Axis (a)
        - Initial True Anomaly (fi)

        Initial Anomalies
        - Initial True Anomaly (fi)
        - Initial Eccentric Anomaly (Ei)
        - Initial Mean Anomaly (Mi)

        Orbital Period
        - Mean Motion (n)
        - Period (period)

        Additional Orbital Parameters
        - Semi-latus Rectum (p)
        - Ascending Node Vector (N)

        Returns
        -------
        None
        """

        r = self.position_norm
        v = self.velocity_norm

        # Specific Energy
        # ---------------------------------------------------------------------
        energy = v**2/2 - self.mu/r # Vis-Viva Equation
        self.orbit_elements.add_parameter("E", energy, "km^2/sec^2",
                                          description="Specific Energy")

        # Specific Angular Momentum Vector
        # ---------------------------------------------------------------------
        ang_momentum = np.cross(self.position,self.velocity)
        self.orbit_elements.add_parameter("h_vector", 
                                          ang_momentum, "km^2/sec",
                                          description="Specific Angular Momentum Vector")
        self.orbit_elements.add_parameter("h",
                                          norm(ang_momentum), "km^2/sec",
                                          description="Magnitude of The Specific Angular Momentum Vector")

        # Eccentricity Vector
        # ---------------------------------------------------------------------
        eccentricity = (1/self.mu)*(np.cross(self.velocity,ang_momentum))-self.position/r
        self.orbit_elements.add_parameter("e_vector",
                                          eccentricity, "",
                                          description="Eccentricity Vector")
        self.orbit_elements.add_parameter("e",
                                          norm(eccentricity), "",
                                          description="Magnitude of The Eccentricity Vector")

        # Inclination
        # ---------------------------------------------------------------------
        incl = np.arccos(np.dot([0,0,1],ang_momentum)/norm(ang_momentum))
        self.orbit_elements.add_parameter("i",
                                          np.degrees(incl), "degrees",
                                          description="Inclination")

        # Ascending Node Vector 
        # ---------------------------------------------------------------------
        node_vec = np.cross([0,0,1],ang_momentum)
        self.orbit_elements.add_parameter("N",
                                          node_vec, "km^2/sec",
                                          description="Ascending Node Vector")

        # Longitude of Ascending Node (Ω)
        # ---------------------------------------------------------------------
        long_ascend_node = np.arccos(np.dot([1,0,0],node_vec)/norm(node_vec))

        if node_vec[1] >= 0.0:
            lan =  np.degrees(norm(long_ascend_node))
        elif node_vec[1] < 0.0:
            lan = np.degrees(2*np.pi - norm(long_ascend_node))
        self.orbit_elements.add_parameter("Omega",
                                          lan, "degrees",
                                          description="Longitude of Ascending Node (Ω)")

        # Argument of Perigee (ω)
        # ---------------------------------------------------------------------
        arg_peri = np.arccos(np.dot(node_vec,eccentricity) / 
                            (norm(eccentricity)*norm(node_vec)))
        # Check eccentricity z component
        if eccentricity[-1] < 0:
            arg_peri = 2*np.pi - arg_peri

        self.orbit_elements.add_parameter("omega",
                                          np.degrees(arg_peri), "degrees",
                                          description="Argument of Perigee (ω)")

        # Semi-latus Rectum
        # ---------------------------------------------------------------------
        p = norm(ang_momentum)*norm(ang_momentum)/self.mu
        self.orbit_elements.add_parameter("p", p, "km",
                                          description="Semi-latus Rectum")

        # Semi-major Axis
        # ---------------------------------------------------------------------
        # a = p/self.mu*(1-norm(eccentricity)*norm(eccentricity))
        a = -self.mu/(2*energy)
        self.orbit_elements.add_parameter("a", a, "km",
                                          description="Semi-major Axis")

        # Mean Motion
        # ---------------------------------------------------------------------
        try:
            n = math.sqrt(self.mu/(a**3))
        except ValueError as e:
            n = np.nan
        self.orbit_elements.add_parameter("n", n, "rad/sec",
                                          description="Mean Motion")

        # Initial True Anomaly
        # ---------------------------------------------------------------------
        # f = np.arccos((p-r)/(r*norm(eccentricity)))
        f = np.arccos(np.dot(eccentricity,self.position) / 
                      (norm(eccentricity)*norm(self.position_norm)))
        # Check orientation of True Anomaly
        if np.dot(self.position,self.velocity) < 0:
            f = 2*np.pi - f
        self.orbit_elements.add_parameter("fi",
                                          np.degrees(f), "degrees",
                                          description="Initial True Anomaly")

        # Initial Eccentric Anomaly
        # ---------------------------------------------------------------------
        ecc_anomaly = np.arccos((norm(eccentricity) + 
                      np.cos(f))/(1+norm(eccentricity)*np.cos(f)))
        self.orbit_elements.add_parameter("Ei",
                                          np.degrees(ecc_anomaly), "degrees",
                                          description="Initial Eccentric Anomaly")

        # Mean Anomaly
        # ---------------------------------------------------------------------
        mean_anomaly = ecc_anomaly - norm(eccentricity)*np.sin(ecc_anomaly)
        self.orbit_elements.add_parameter("Mi",
                                          np.degrees(mean_anomaly), "degrees",
                                          description="Initial Mean Anomaly")

        # Orbital Period
        # ---------------------------------------------------------------------
        period = 2*np.pi / n
        self.orbit_elements.add_parameter("period", period, "secs",
                                          description="Period")
 