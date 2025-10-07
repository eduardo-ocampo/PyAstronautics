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

from pyastronautics.astrodynamics.three_body_problem import planar_lagrange_points
class TestPlanarLagrangePoints:
    def test_lagrange_points_init(self):
        """Test the initialization of the lagrange_points class."""
        mu = 0.1
        lagrange = planar_lagrange_points(mu)

        # Check if the mass ratio is correctly set
        assert lagrange.mu == mu, f"Expected mu = {mu}, but got {lagrange.mu}"

        # Check if the initial Lagrange points are correctly set to 0
        assert lagrange.l1y == lagrange.l2y == lagrange.l3y == 0

    def test_colinear_points(self):
        """Test the calculation of collinear Lagrange points (L1, L2, L3)."""
        mu = 0.1
        lagrange = planar_lagrange_points(mu)
        
        # Get collinear points
        lagrange.colinear_points()
        
        # Lagrange points should be within a reasonable range (not NaN or infinite)
        assert np.all(np.isfinite([lagrange.l1x, lagrange.l2x, lagrange.l3x])), "Collinear points contain invalid values"
        
        # Optionally, check for expected ranges (based on mass ratio)
        # The L1 point should be near the first primary body, etc.
        assert lagrange.l1x < 1, f"L1x should be less than 1, but got {lagrange.l1x}"
        assert lagrange.l2x > 1, f"L2x should be greater than 1, but got {lagrange.l2x}"

        # Check exact known values for mu = 0.1
        assert np.isclose(lagrange.l1x, 0.609035110, atol=1e-8), f"L1x is incorrect: {lagrange.l1x}"
        assert np.isclose(lagrange.l2x, 1.259699832, atol=1e-8), f"L2x is incorrect: {lagrange.l2x}"
        assert np.isclose(lagrange.l3x, -1.041608908, atol=1e-8), f"L3x is incorrect: {lagrange.l3x}"

    def test_triangular_points(self):
        """Test the calculation of triangular Lagrange points (L4 and L5)."""
        mu = 0.1
        lagrange = planar_lagrange_points(mu)
        
        # Get triangular points
        lagrange.triangular_points()
        
        # Check if L4 and L5 are correctly calculated (should be symmetric)
        assert np.isclose(lagrange.l4x, 0.5 - mu), f"L4x is incorrect: {lagrange.l4x}"
        assert np.isclose(lagrange.l4y, np.sqrt(3)/2), f"L4y is incorrect: {lagrange.l4y}"
        assert np.isclose(lagrange.l4, np.sqrt(lagrange.l4x**2 + lagrange.l4y**2)), f"L4 distance is incorrect: {lagrange.l4}"

        # Check that L5 is symmetric to L4
        assert np.isclose(lagrange.l5x, lagrange.l4x), f"L5x is incorrect: {lagrange.l5x}"
        assert np.isclose(lagrange.l5y, -lagrange.l4y), f"L5y is incorrect: {lagrange.l5y}"
        assert np.isclose(lagrange.l5, lagrange.l4), f"L5 distance is incorrect: {lagrange.l5}"

    def test_get_points(self):
        """Test the calculation of both collinear and triangular points using get_points."""
        mu = 0.1
        lagrange = planar_lagrange_points(mu)
        
        # Call get_points to compute both sets of points
        lagrange.get_points()
        
        # Check if all the points are computed correctly (should not be NaN or infinite)
        assert np.all(np.isfinite([lagrange.l1x, lagrange.l2x, lagrange.l3x, lagrange.l4x, lagrange.l5x])), "One or more Lagrange points contain invalid values"

        # Check exact known values for mu = 0.1
        assert np.isclose(lagrange.l1x, 0.609035110, atol=1e-8), f"L1x is incorrect: {lagrange.l1x}"
        assert np.isclose(lagrange.l2x, 1.259699832, atol=1e-8), f"L2x is incorrect: {lagrange.l2x}"
        assert np.isclose(lagrange.l3x, -1.041608908, atol=1e-8), f"L3x is incorrect: {lagrange.l3x}"

    def test_root_equation(self):
        """Test the root equation used to solve for collinear Lagrange points."""
        mu = 0.1
        lagrange = planar_lagrange_points(mu)
        
        # Known L1 point for mu = 0.1
        x_L1 = 0.609035110
        
        # The force function (root_equation) should be zero at the L1 point
        residual = lagrange.root_equation(x_L1)
    
        assert np.isclose(residual, 0, atol=1e-8), f"Expected near-zero force at L1 (x = {x_L1}), got: {residual}"
    
    def test_colinear_approximation(self):
        """Test the colinear approximation method for Lagrange points."""
        mu = 0.1
        lagrange = planar_lagrange_points(mu)
        
        # Get the approximation for L1, L2, L3
        approx = lagrange.colinear_approximation()
        
        # Check if the approximation provides reasonable results (based on known theoretical results)
        assert len(approx) == 3, "Collinear approximation should return 3 points."
        
        # Verify that the approximations are within reasonable bounds for the mass ratio
        assert approx[0] < 1, f"L1 approximation should be less than 1, but got {approx[0]}"
        assert approx[1] > 1, f"L2 approximation should be greater than 1, but got {approx[1]}"
        assert approx[2] < 0, f"L3 approximation should be less than 0, but got {approx[2]}"

        # Compare approximated values to known expected ones
        expected = [0.589276749, 1.210723250, -1.041666667]
        for i, val in enumerate(approx):
            assert np.isclose(val, expected[i], atol=1e-2), f"Approximation for L{i+1} is off: {val} (expected ~{expected[i]})"
