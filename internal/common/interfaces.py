# abstract drug delivery system with any form
class AbstractDDS(object):
    def __init__(self):
        super(AbstractDDS, self).__init__()

    def common_ode(self, c, t):
        """ Setup the system of n equations"""
        # kp = permeability coefficient (um/hour)
        # kpn = permeability coefficient for the n^th compartment
        pass

    def single_particle_ode(self, c, t):
        """This single ODE is solved if only one compartment remains"""
        # kpn = permeability coefficient (um/hour)
        pass

    def two_particles_ode(self, c, t):
        """ This system is solved if two compartments remain"""
        # kpn = permeability coefficient (um/hour)
        pass