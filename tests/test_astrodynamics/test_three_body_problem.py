""" 
   Copyright 2024 Eduardo Ocampo
   https://github.com/eduardo-ocampo/PyAstronautics
"""

import pytest
import numpy as np

from pyastronautics.astrodynamics.three_body_problem import CR3BP

class TestCR3BP:

    @pytest.fixture(autouse=True)
    def setup_method(self):
        """Set up the spacecraft model instance for testing."""

        initial_position = [ 0.50, 0.50, 0.0]
        initial_velocity = [-0.05, 0.10, 0.0]
        
        # Spacecraft
        self.sc = CR3BP(initial_position,initial_velocity)
        # Earth-Moon Mass Ratio 
        self.sc.mu = 0.012150515586657583

        # Time
        time_num = 1000
        self.sc.time = np.linspace(0, 2*np.pi*4, time_num)

        # Numerical Analsysi Setup
        self.sc.rel_tol = 1e-12
        self.sc.abs_tol = 1e-13

    def test_solve_trajectory_success(self):

        self.sc.solve_non_dim_trajectory()
        
        # Check that the results are populated
        assert self.sc.numerical_position is not None
        assert self.sc.numerical_velocity is not None

        # Check the shape of the output arrays
        assert self.sc.numerical_position.shape == (len(self.sc.time), 3)
        assert self.sc.numerical_velocity.shape == (len(self.sc.time), 3)

    def test_solve_trajectory_final_state(self):

        self.sc.solve_non_dim_trajectory()

        # Known final state
        expected_state = np.array([-0.00131049,  0.53114043,  0.,
                                   -0.59227612, -0.62395045, 0.])

        assert np.allclose(self.sc.final_state, expected_state)

    def test_solve_trajectory_empty_time(self):
 
        self.sc.time = []
        with pytest.raises(IndexError, match="list index out of range"):
            self.sc.solve_non_dim_trajectory()

    def test_solve_trajectory_no_time(self):
        
        # Remove the time attribute
        del self.sc.time
        with pytest.raises(ValueError, match="Attribute 'time' must be defined before calling solve_trajectory."):
            self.sc.solve_non_dim_trajectory()

class TestForbiddenRegions:
    def test_forbidden_region_default(self):
        # Use some reasonable defaults
        jacobi_max = 3.0
        mass_ratio = 0.1
        x_range = [-1.5, 1.5]
        y_range = [-1.5, 1.5]
        linspace_num = 200

        x_vals, y_vals, jc = CR3BP.forbidden_region(jacobi_max, mass_ratio, x_range, y_range, linspace_num)

        # Check the shape of the output arrays
        assert x_vals.shape == (linspace_num,)
        assert y_vals.shape == (linspace_num,)
        assert jc.shape == (linspace_num, linspace_num)

        # Ensure no values exceed the jacobi_max (they should be NaN)
        assert np.all(np.isnan(jc[jc > jacobi_max]))

    def test_forbidden_region_with_custom_range(self):
        # Test with custom x_range and y_range
        jacobi_max = 3.0
        mass_ratio = 0.1
        x_range = [-2.0, 2.0]
        y_range = [-2.0, 2.0]
        linspace_num = 100

        x_vals, y_vals, jc = CR3BP.forbidden_region(jacobi_max, mass_ratio, x_range, y_range, linspace_num)

        # Check the shape of the output arrays
        assert x_vals.shape == (linspace_num,)
        assert y_vals.shape == (linspace_num,)
        assert jc.shape == (linspace_num, linspace_num)

    def test_forbidden_region_with_high_jacobi_max(self):
        # Test with a high jacobi_max to see if it affects NaN values
        jacobi_max = 450.0  # Set high value for jacobi_max
        mass_ratio = 0.1
        x_range = [-1.5, 1.5]
        y_range = [-1.5, 1.5]
        linspace_num = 200

        x_vals, y_vals, jc = CR3BP.forbidden_region(jacobi_max, mass_ratio, x_range, y_range, linspace_num)

        # Verify that the Jacobi constant is below the given threshold (no NaNs)
        assert np.all(~np.isnan(jc))  # There should be no NaNs if jacobi_max is high

    # Test for calculate_jacobi function
    def test_calculate_jacobi(self):
        # Known test case
        x = 0.5
        y = 0.5
        vx = 0.1
        vy = 0.1
        mass_ratio = 0.1
        
        expected_jacobi = (x**2 + y**2) + 2*(1 - mass_ratio) / np.sqrt((x + mass_ratio)**2 + y**2) + 2 * mass_ratio / np.sqrt((x - 1 + mass_ratio)**2 + y**2) - (vx**2 + vy**2)
        
        jacobi = CR3BP.calculate_jacobi(x, y, vx, vy, mass_ratio)

        # Check if the result is as expected
        assert np.isclose(jacobi, expected_jacobi), f"Expected {expected_jacobi}, but got {jacobi}"

    def test_calculate_jacobi_zero_velocity(self):
        # Test case with zero velocity (at rest)
        x = 0.5
        y = 0.5
        vx = 0.0
        vy = 0.0
        mass_ratio = 0.1
        
        # Calculate Jacobi constant manually for vx=0, vy=0
        expected_jacobi = (x**2 + y**2) + 2*(1 - mass_ratio) / np.sqrt((x + mass_ratio)**2 + y**2) + 2 * mass_ratio / np.sqrt((x - 1 + mass_ratio)**2 + y**2)
        
        jacobi = CR3BP.calculate_jacobi(x, y, vx, vy, mass_ratio)

        # Check if the result is as expected
        assert np.isclose(jacobi, expected_jacobi), f"Expected {expected_jacobi}, but got {jacobi}"

    def test_calculate_jacobi_edge_case(self):
        # Test case with edge position and velocity
        x = 0.0
        y = 0.0
        vx = 1.0
        vy = 1.0
        mass_ratio = 0.5
        
        # Calculate Jacobi constant manually for these values
        r1 = np.sqrt((x + mass_ratio)**2 + y**2)
        r2 = np.sqrt((x - 1 + mass_ratio)**2 + y**2)
        expected_jacobi = (x**2 + y**2) + 2 * (1 - mass_ratio) / r1 + 2 * mass_ratio / r2 - (vx**2 + vy**2)
        
        jacobi = CR3BP.calculate_jacobi(x, y, vx, vy, mass_ratio)

        # Check if the result is as expected
        assert np.isclose(jacobi, expected_jacobi), f"Expected {expected_jacobi}, but got {jacobi}"
