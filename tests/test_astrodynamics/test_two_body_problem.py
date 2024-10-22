""" 
   Copyright 2024 Eduardo Ocampo
   https://github.com/eduardo-ocampo/PyAstronautics
"""

import math
import pytest
import numpy as np

from pyastronautics.astrodynamics.two_body_problem import TwoBodyModel

class TestTwoBody:

    @pytest.fixture(autouse=True)
    def setup_method(self):
        """Set up the spacecraft model instance for testing."""

        # Earth 
        mu = 398600 # km^3/sec^2

        position = [5000, 100, 0] # km
        vy = math.sqrt(mu/position[0])+1
        velocity = [1, vy, 1] # km/s

        # Spacecraft
        self.sc = TwoBodyModel(position, velocity)
        # Earth
        self.sc.mu = mu

        # Create Orbital Parameters
        self.sc.calc_orbit_elements()

    def test_create_two_body(self):
        
        assert self.sc.orbit_elements.position.value == [5000, 100, 0] # km
        assert self.sc.orbit_elements.position.unit == "km"
        assert np.allclose(self.sc.orbit_elements.velocity.value, [1, 9.9286057, 1]) # km/s
        assert self.sc.orbit_elements.velocity.unit == "km/s"

    def test_create_two_body_with_invalid_position(self):
        with pytest.raises(TypeError, match="position must be a list."):
            TwoBodyModel(position=np.array([100.0, 200.0, 300.0]),
                        velocity=[10.0, 15.0, 20.0])
            
    def test_create_two_body_with_invalid_velocity(self):
        with pytest.raises(TypeError, match="velocity must be a list."):
            TwoBodyModel(position=[100.0, 200.0, 300.0],
                         velocity=np.array([10.0, 15.0, 20.0]))

    # Test Calculated Orbital Parameters
    # against known analytical solutions        
    def test_conservation_parameters(self):
        print(self.sc.orbit_elements.parameters)
        assert np.isclose(self.sc.orbit_elements.E.value, -29.41554340)
        assert np.isclose(self.sc.orbit_elements.h.value, 49794.795711)
        assert np.allclose(self.sc.orbit_elements.h_vector.value, [100, -5000, 49543.0285711])
        
    def test_keplerian_orbit_elements(self):
        assert np.allclose(self.sc.orbit_elements.e_vector.value, [0.24679464, -0.14403758, -0.01503476])
        assert np.isclose(self.sc.orbit_elements.e.value, 0.286147621)
        assert np.isclose(self.sc.orbit_elements.i.value, 5.7640578743)
        assert np.isclose(self.sc.orbit_elements.Omega.value, 1.14576283)
        assert np.isclose(self.sc.orbit_elements.omega.value, 328.4556365)
        assert np.isclose(self.sc.orbit_elements.a.value, 6775.337042869201)

    def test_initial_anomalies(self):
        assert np.isclose(self.sc.orbit_elements.fi.value, 31.544363499)
        assert np.isclose(self.sc.orbit_elements.Ei.value, 23.7661274)
        assert np.isclose(self.sc.orbit_elements.Mi.value, 17.15885116835)
    
    def test_additional_orbital_parameters(self):
        assert np.allclose(self.sc.orbit_elements.N.value, [5000, 100, 0])
        assert np.isclose(self.sc.orbit_elements.p.value, 6220.5693220)
        assert np.isclose(self.sc.orbit_elements.n.value, 0.00113206419)
        assert np.isclose(self.sc.orbit_elements.period.value, 5550.1834355392)

class TestTwoBodyModel:

    @pytest.fixture(autouse=True)
    def setup_method(self):
        """Set up the spacecraft model instance for testing."""

        position = [5000, 100, 0] # km
        velocity = [1, 10, 5] # km/s
        
        # Spacecraft
        self.sc = TwoBodyModel(position, velocity)
        # Earth
        self.sc.mu = 398600 # km^3/sec^2

        # Known orbital period
        period = 12969.97314383982 # sec

        self.sc.time = np.arange(0, 20*period, 15*60)

    def test_solve_trajectory_success(self):

        self.sc.solve_trajectory()
        
        # Check that the results are populated
        assert self.sc.numerical_position is not None
        assert self.sc.numerical_velocity is not None

        # Check the shape of the output arrays
        assert self.sc.numerical_position.shape == (len(self.sc.time), 3)
        assert self.sc.numerical_velocity.shape == (len(self.sc.time), 3)

    def test_solve_trajectory_final_state(self):

        self.sc.solve_trajectory()

        # Known final state
        expected_state = np.array([4.48191388e+03, -1.85756502e+03, -9.75552752e+02,
                                   4.15484717e+00, 9.41162692e+00,  4.67361221e+00])

        assert np.allclose(self.sc.final_state, expected_state)

    def test_solve_trajectory_empty_time(self):
 
        self.sc.time = []
        with pytest.raises(IndexError, match="list index out of range"):
            self.sc.solve_trajectory()

    def test_solve_trajectory_no_time(self):
        
        # Remove the time attribute
        del self.sc.time
        with pytest.raises(ValueError, match="Attribute 'time' must be defined before calling solve_trajectory."):
            self.sc.solve_trajectory()
