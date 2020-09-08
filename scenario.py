import oceanStateSplit

class Scenario:

    # Floating: If False, boundaries will be constant as their initial values unless bounds are given.
    # If True, boundaries will float unless bounds are given.
    E_floating = False
    U_floating = False
    V_floating = False


    def initialize(self, sp):
        print("Scenario must override initialize method!")

    def getOs(self):
        print("Scenario must override getOs method!")

    # The initPassiveTracer method sets initial values for the passive tracer variable os.X.
    # If not overridden, os.X will contain all zeros.
    def initPassiveTracer(self):
        return

    # The setBounds method gets dimensions and depth matrix for the full domain as input, and
    # should return a dictionary containing boundaries for those variables for which boundary
    # conditions are set. Each dictionary entry is keyed by the variable name, and should be a
    # list of numpy arrays containing bounds for [left, bottom, right, upper] sides.
    def setBounds(self, imax, jmax, kmax, fullDepth, t, os):
        print("Scenario must override setBounds method!")


    # The setAtmo method gets dimensions for the full domain as input, and should return
    # matrices containing the U and V components of the wind field.
    # If the scenario needs no wind, it doesn't need to override this method.
    def setAtmo(self, imax, jmax, t, os):
        return None, None


    # The getFreshwater method should return:
    # - list of river positions (each an i,j tuple)
    # - list of river directions (each a string: east, west, top, bottom)
    # - list of flow properties (each a Q, T, S tuple)
    def getFreshwater(self, imax, jmax, t, os):
        return [], [], []
