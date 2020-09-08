from rkCalcUVE import *

def getEmpties(os):
    return np.empty(os.U.shape), np.empty(os.V.shape), np.empty(os.E.shape)

# This function integrates the U, V and E states using a simple forward Euler scheme.
def forwardEuler(os, sp, t):
    # Get derivative of model states:
    f_U, f_V, f_E = rkCalcUVE(os, sp, t, sp.dt)
    # Update model state:
    os.U = os.U + sp.dt*f_U
    os.V = os.V + sp.dt*f_V
    os.E = os.E + sp.dt*f_E

# This function integrates the U, V and E states using the Runge-Kutta 4 integration
# scheme. RK4 is a fourth order scheme that provides better stability and allows
# longer time steps than a simple forward Euler scheme.
def rk4Integration(os, sp, t):
    h = sp.dt

    # Copy original state values:
    U0, V0, E0 = getEmpties(os)
    U0[:] = os.U[:]
    V0[:] = os.V[:]
    E0[:] = os.E[:]

    # Calculate k1:
    f_U, f_V, f_E = rkCalcUVE(os, sp, t, h)
    k1_U = h*f_U
    k1_V = h*f_V
    k1_E = h*f_E

    # Calculate k2:
    os.U = os.U + 0.5*k1_U
    os.V = os.V + 0.5*k1_V
    os.E = os.E + 0.5*k1_E
    f_U, f_V, f_E = rkCalcUVE(os, sp, t, h)
    k2_U = h * f_U
    k2_V = h * f_V
    k2_E = h * f_E

    # Calculate k3:
    os.U = os.U + 0.5*k2_U
    os.V = os.V + 0.5*k2_V
    os.E = os.E + 0.5*k2_E
    f_U, f_V, f_E = rkCalcUVE(os, sp, t, h)
    k3_U = h * f_U
    k3_V = h * f_V
    k3_E = h * f_E

    # Calculate k4:
    os.U = os.U + k3_U
    os.V = os.V + k3_V
    os.E = os.E + k3_E
    f_U, f_V, f_E = rkCalcUVE(os, sp, h, t)
    k4_U = h * f_U
    k4_V = h * f_V
    k4_E = h * f_E

    # Finally update model state:
    os.U = U0 + (1.0 / 6.0) * (k1_U + 2 * k2_U + 2 * k3_U + k4_U)
    os.V = V0 + (1.0 / 6.0) * (k1_V + 2 * k2_V + 2 * k3_V + k4_V)
    os.E = E0 + (1.0 / 6.0) * (k1_E + 2 * k2_E + 2 * k3_E + k4_E)
