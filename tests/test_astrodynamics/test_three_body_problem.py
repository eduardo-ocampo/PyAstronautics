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
