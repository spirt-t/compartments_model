from internal.common.interfaces import AbstractDDS


class Slab(AbstractDDS):


    def __init__(self):
        # self.time = time
    pass

    def common_ode(self, c, t):
        """ Setup the system of n equations"""
        # kp = permeability coefficient (um/hour)
        # kpn = permeability coefficient for the n^th compartment
        kp = N * D / R0
        kpn = 2 * D / ((N - n + 2) * R0 / N - v * t)
        kpn2 = 2 * D / (((N - n + 1) * R0 / N - v * t))
        if kpn2 >= 2 * D / (tolerance * R0 / N):
            kpn2 = 2 * D / (tolerance * R0 / N)
        dc = []
        dc.append(kp * N / R0 * (-c[0] + c[1]))
        for i in range(1, n - 2):
            dc.append(kp * N / R0 * (c[i - 1] - 2 * c[i] + c[i + 1]))
        dc.append(
            kp * N / R0 * (c[n - 3] - c[n - 2])
            + kpn * N / R0 * (c[n - 1] - c[n - 2])
        )
        dc.append(1 / (R0 - v * t - R0 * (n - 1) / N)
                  * (kpn * (c[n - 2] - c[n - 1]) - kpn2 * c[n - 1])
                  )
        return np.array(dc)

    def single_particle_ode(self, c, t):
        """This single ODE is solved if only one compartment remains"""
        # kpn = permeability coefficient (um/hour)
        kpn2 = 2 * D / (((N - n + 1) * R0 / N - v * t))
        if kpn2 >= 2 * D / (tolerance * R0 / N):
            kpn2 = 2 * D / (tolerance * R0 / N)
        dc = -kpn2 * c[n - 1] / (R0 - v * t)
        return np.array(dc)

    def two_particles_ode(self, c, t):
        """ This system is solved if two compartments remain"""
        # kpn = permeability coefficient (um/hour)
        dc = []
        kpn = 2 * D / (((N - n + 2) * R0 / N - v * t))
        kpn2 = 2 * D / (((N - n + 1) * R0 / N - v * t))
        if kpn2 >= 2 * D / (tolerance * R0 / N):
            kpn2 = 2 * D / (tolerance * R0 / N)
        dc.append(kpn * N / R0 * (-c[0] + c[1]))
        dc.append(1 / (R0 - v * t - R0 * (n - 1) / N) * (kpn * (-c[1] + c[0]) - kpn2 * c[1]))
        return np.array(dc)
