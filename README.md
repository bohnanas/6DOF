# 6DOF
Big thanks to Ben Dickinson for showing how to implement this in Python!

Simulation of a 6-DOF rigid body aicraft

Verified EOM against NASA test cases 1, 2, and 3:
Case 1: Dragless sphere dropped from altitude (verifies gravitational and translational EOM)
Case 2: Dragless brick dropped from altitude with no damping (verifies rotational EOM)
Case 3: Dragless brick dropped from altitude with damping (verifies inertial coupling)

Added X-15 longitudinal dynamics, lateral to come

Added FlightGear visualization 