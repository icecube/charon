"""
Author  : C.A. Arguelles,  Q.R. Liu
"""
import numpy as np
import scipy.special as spe
import scipy.interpolate as interpolate
import matplotlib as mpl
import socket

if "cobalt" in socket.gethostname():
    mpl.use("Agg", warn=False)
import matplotlib.pyplot as plt
import physicsconstants as PC

pc = PC.PhysicsConstants()


def HelmFormFactor(q, A):
    """Calculates Helm's nuclcear form factor.

    Ref : Spin-independent elastic WIMP scattering and the DAMA annual modulation signal. M. Fairbairn.
    arXiv : 0808.0704
    Ec. between ec. (1) and (2).

    @type  q            :      float
    @param q            :      transfer momentum [eV]
    @type  A            :      float
    @param A            :      mass number

    @rtype              :      float
    @return             :      Helm's form factor
    """
    s = 1.0 * pc.fermi
    R = 1.2 * A ** (1 / 3) * pc.fermi
    r = np.sqrt(R ** 2 - 5.0 * s ** 2)

    return (
        3.0
        * np.exp(-((q * s) ** 2) / 2)
        * (np.sin(q * r) - q * r * np.cos(q * r))
        / (q * r) ** 3
    )


def Woods_SaxonFormFactor(q, A):
    """Calculates Wood-Saxong's nuclcear form factor.

    Ref : SUPERSYMMETRIC DARK MATTER, G. JUNGMAN.
    Ec. 7.31 p. 271.

    @type  q            :      float
    @param q            :      transfer momentum [eV]
    @type  A            :      float
    @param A            :      mass number


    @rtype              :      float
    @return             :      Woods-Saxon's form factor
    """
    s = 1.0 * pc.fermi
    R = 1.2 * A ** (1 / 3) * pc.fermi
    r = np.sqrt(R ** 2 - 5.0 * s ** 2)

    return (3.0 * spe.sph_jn(q * r) / (q * r)) * np.exp(-((q * r) ** 2))


def JungmanFormFactor(element_id, DM_mass):
    """Calculates Jungman's form factor according to ec. 9.24 p.301

    Ref : SUPERSYMMETRIC DARK MATTER, G. JUNGMAN.
    Ec. 9.24 p. 301.

    @type  i            :      int
    @param i            :      element array number.
    @type  DM_mass      :      float
    @param DM_mass      :      DM mass [eV].


    @rtype              :      float
    @return             :      Jungman's form factor
    """
    # USING Table 10, p. 301 from Jungman Rev.
    elements = {
        0: "H",
        1: "He",
        2: "C",
        3: "N",
        4: "O",
        5: "Ne",
        6: "Mg",
        7: "Si",
        8: "S",
        9: "Fe",
    }  # names and order
    m_e_GeV = [1.0, 18.2, 61.6, 75.2, 75.2, 75.2, 71.7, 71.7, 57.0, 29.3]  # [GeV]
    m_e = [x * pc.GeV for x in m_e_GeV]
    F_inf = [1.0, 0.986, 0.788, 0.613, 0.613, 0.613, 0.281, 0.281, 0.101, 0.00677]
    alpha = [1.0, 1.58, 2.69, 2.69, 2.69, 2.69, 2.97, 2.97, 3.1, 3.36]

    DM_mass_GeV = DM_mass / pc.GeV

    i = element_id
    if i == 0:  # Ref : SUPERSYMMETRIC DARK MATTER, G. JUNGMAN. around 9.24 F_H=1
        return 1
    else:
        return F_inf[i] + (1.0 - F_inf[i]) * np.exp(
            -((np.log(DM_mass_GeV) / np.log(m_e_GeV[i])) ** alpha[i])
        )


def GouldAuxF1(A, a, b, eta):
    """Calculates part of the DM annihilation as estimated by Gould. This is
    the first termn in {}-brackets on eq. 2.18.

    Ref : Cosmological density of WIMPs from solar and terrestrial annihilations. A. Gould.
    Astrophysical Journal, 388:338-344, 1992 April.

    @type  A      :      float
    @param A      :      Gould's parameter A.
    @type  a      :      float
    @param a      :      Gould's parameter a.
    @type  b      :      float
    @param b      :      Gould's parameter b.
    @type  eta    :      float
    @param eta    :      Gould's parameter eta.

    @rtype        :      float
    @return       :      part of Gould's annihilation formulae.
    """
    A_hat = A * np.sqrt(1.0 + a)
    eta_hat = eta / np.sqrt(1.0 + a)

    A_hat_plus = A_hat + eta_hat
    A_hat_minus = A_hat - eta_hat

    return (A_hat_plus * A_hat_minus - 0.5 - (1 + a) / (a - b)) * (
        spe.erf(A_hat_plus) - spe.erf(A_hat_minus)
    ) + (1.0 / np.sqrt(np.pi)) * (
        A_hat_minus * np.exp(-(A_hat_plus ** 2))
        - A_hat_plus * np.exp(-(A_hat_minus ** 2))
    )


def GouldAuxF2(A, a, b, eta):
    """Calculates part of the DM annihilation as estimated by Gould. This is
    the first termn in {}-brackets on eq. 2.18.

    Ref : Cosmological density of WIMPs from solar and terrestrial annihilations. A. Gould.
    Astrophysical Journal, 388:338-344, 1992 April.

    @type  A      :      float
    @param A      :      Gould's parameter A.
    @type  a      :      float
    @param a      :      Gould's parameter a.
    @type  b      :      float
    @param b      :      Gould's parameter b.
    @type  eta    :      float
    @param eta    :      Gould's parameter eta.

    @rtype        :      float
    @return       :      part of Gould's annihilation formulae.
    """
    A_hat = A * np.sqrt(1.0 + b)
    eta_hat = eta / np.sqrt(1.0 + b)

    A_hat_plus = A_hat + eta_hat
    A_hat_minus = A_hat - eta_hat

    return (
        2.0 * spe.erf(eta_hat) - spe.erf(A_hat_plus) + spe.erf(A_hat_minus)
    ) * np.exp(-(a - b) * A ** 2)


def DMSunCaptureRateGould(
    DM_mass,
    DM_cs,
    DM_rho=0.3,
    vel_rot=220.0,
    model=None,
    n=None,
    mass_num=None,
    mass_gr_per_mol=None,
    eps=None,
):
    """Calculates DM capture rate in the Sun using Gould's formula.

    Ref : Cosmological density of WIMPs from solar and terrestrial annihilations. A. Gould.
    Astrophysical Journal, 388:338-344, 1992 April.

    @type  DM_mass        :      float
    @param DM_mass        :      DM mass [GeV]
    @type  DM_cs          :      float
    @param DM_cs          :      DM cross section [cm^2]
    @type  DM_rho         :      float
    @param DM_rho         :      local DM density [GeV/cm^3]
    @type  vel_rot        :      float
    @param vel_rot        :      local DM velocity[km/s]
    @type  model          :      str
    @param model          :      model name 'bs05_agsop.dat', 'struct_b16_agss09.dat', 'struct_b16_gs98.dat'
    @type  n              :      int
    @param n              :      number of elements
    @type  mass_num       :      list
    @param mass_num       :      mass number of the n elements
    @type  mass_gr_per_mol:      list
    @param mass_gr_per_mol:      mass per mole of the n elements [g]
    @type  eps            :      list
    @param eps            :      abundance of the n elements


    @type                 :      float
    @return               :      DM capture rate at the Sun. [s^-1]
    """

    DM_mass = DM_mass * pc.GeV
    DM_cs = DM_cs * pc.cm ** 2
    if any(e != None for e in [n, mass_num, eps, mass_gr_per_mol]):
        pass

    elif model != None:
        path = "../charon/models/"
        data = np.genfromtxt(path + model + ".dat")
        m = data[:, 0]
        m_shell = np.append(m[0], (m[1:] - m[:-1]))

        if model in ["struct_b16_agss09", "struct_b16_gs98"]:
            elements = {
                0: "H1",
                1: "He4",
                2: "He3",
                3: "C12",
                4: "C13",
                5: "N14",
                6: "N15",
                7: "O16",
                8: "O17",
                9: "O18",
                10: "Ne",
                11: "Na",
                12: "Mg",
                13: "Al",
                14: "Si",
                15: "P",
                16: "S",
                17: "Cl",
                18: "Ar",
                19: "K",
                20: "Ca",
                21: "Sc",
                22: "Ti",
                23: "V",
                24: "Cr",
                25: "Mn",
                26: "Fe",
                27: "Co",
                28: "Ni",
            }
            mass_num = [
                1.0,
                4.0,
                3.0,
                12.0,
                13.0,
                14.0,
                15.0,
                16.0,
                17.0,
                18.0,
                20.0,
                23.0,
                24.0,
                27.0,
                28.0,
                31.0,
                32.0,
                35.0,
                40.0,
                39.0,
                40.0,
                45.0,
                48.0,
                51.0,
                52.0,
                55.0,
                56.0,
                59.0,
                58.0,
            ]
            mass_gr_per_mol = [
                1.0079,
                4.0026,
                3.0160,
                12.0000,
                13.0034,
                14.0031,
                15.0001,
                15.9950,
                16.9991,
                17.9992,
                20.1797,
                22.9898,
                24.3050,
                26.9815,
                28.0850,
                30.9738,
                32.0600,
                35.4500,
                39.9481,
                39.0983,
                40.0784,
                44.9559,
                47.8671,
                50.9415,
                51.9962,
                54.9380,
                55.8452,
                58.9332,
                58.6934,
            ]
            eps = [np.sum(data[:, i + 6] * m_shell) for i in range(len(elements))]

        elif model == "bs05_agsop":
            elements = {0: "H", 1: "He4", 2: "He3", 3: "C12", 4: "N14", 5: "O16"}
            mass_num = [1.0, 4.0, 3.0, 12.0, 14.0, 16.0]
            mass_gr_per_mol = [1.0079, 4.0026, 3.0160, 12.0000, 14.0031, 15.9950]
            eps = [np.sum(data[:, i + 6] * m) / m_sum for i in range(len(elements))]
        n = len(elements)
    else:
        # Jungman
        elements = {
            0: "H",
            1: "He",
            2: "C",
            3: "N",
            4: "O",
            5: "Ne",
            6: "Mg",
            7: "Si",
            8: "S",
            9: "Fe",
        }  # names
        mass_num = [
            1.0,
            4.0,
            12.0,
            14.0,
            16.0,
            20.0,
            24.0,
            28.0,
            32.0,
            56.0,
        ]  # A : mass number
        mass_gr_per_mol = [
            1.0079,
            4.0026,
            12.0107,
            14.0067,
            15.9994,
            20.1797,
            24.3051,
            28.0855,
            32.0655,
            55.8452,
        ]  # gr mol^-1
        eps = [
            0.772,
            0.209,
            3.87e-3,
            9.4e-4,
            8.55e-3,
            1.51e-3,
            7.39e-4,
            8.13e-4,
            4.65e-4,
            1.46e-3,
        ]  # relative aboundances in the Sun from Table 8 Jungman p.131
        n = len(elements)

    # input data
    mass_eV = [m * pc.gr / pc.Na for m in mass_gr_per_mol]
    atom_radius = [(1.2 * np.power(A, 1.0 / 3.0) * pc.fermi) for A in mass_num]

    energy_element = [3.0 / (2.0 * mass_eV[i] * atom_radius[i] ** 2) for i in range(n)]

    DM_rho = DM_rho * pc.GeV / pc.cm ** 3  # new profile pdg2018
    sun_mass = 1.9891 * 1.0e30 * pc.kg

    eta = 1.0
    # Gould's velocities definitions
    vel_rot_km = vel_rot
    vel_char = vel_rot_km / eta
    vel_rot = vel_rot_km * pc.km / pc.sec  # local DM rel speed.
    # Jungman velocities definitions
    vel_bar = vel_rot_km * pc.km / pc.sec
    vel_star = np.sqrt(3.0 / 2.0) * vel_bar / eta

    vel_surface = 795.0 * pc.km / pc.sec
    vel_center = 1354.0 * pc.km / pc.sec

    # special assumptions
    q_i = [m for m in mass_eV]

    # Gould's auxiliaries variables
    beta_plus = [(4.0 * DM_mass * m_i) / (DM_mass + m_i) ** 2 for m_i in mass_eV]
    beta_minus = [(4.0 * DM_mass * m_i) / (DM_mass - m_i) ** 2 for m_i in mass_eV]

    a = [DM_mass * vel_bar ** 2 / (2.0 * E) for E in energy_element]
    b = [a[i] * beta_plus[i] for i in range(n)]

    A_center = [np.sqrt(b_min * (vel_center / vel_bar) ** 2) for b_min in beta_minus]
    A_surface = [np.sqrt(b_min * (vel_surface / vel_bar) ** 2) for b_min in beta_minus]

    eta_up_hat = [eta / np.sqrt(1.0 + aa) for aa in a]
    eta_dn_hat = [eta / np.sqrt(1.0 + bb) for bb in b]

    reduce_mass = lambda m1, m2: m1 * m2 / (m1 + m2)

    p_mass = 938.272 * pc.MeV  # proton mass
    q_i = [
        (mass_num[i] * reduce_mass(DM_mass, mass_eV[i]) / reduce_mass(DM_mass, p_mass))
        / np.sqrt(beta_plus[i] * mass_eV[i] * DM_mass)
        for i in range(n)
    ]
    sigma_0 = [DM_cs for i in range(n)]

    T_0 = [
        (sigma_0[i] * DM_rho * sun_mass * vel_star / (2.0 * eta))
        * (q_i[i] ** 2 * eps[i] / a[i])
        for i in range(n)
    ]

    T_1 = [
        (
            2.0
            * np.exp(-a[i] * eta_up_hat[i] ** 2)
            * spe.erf(eta_up_hat[i])
            / np.sqrt(1.0 + a[i])
        )
        for i in range(n)
    ]
    T_2 = [
        -np.exp(-a[i] * eta_up_hat[i] ** 2)
        / ((A_center[i] ** 2 - A_surface[i] ** 2) * np.power(1.0 + a[i], 3.0 / 2.0))
        * (
            GouldAuxF1(A_center[i], a[i], b[i], eta)
            - GouldAuxF1(A_surface[i], a[i], b[i], eta)
        )
        for i in range(n)
    ]
    T_3 = [
        np.exp(-b[i] * eta_dn_hat[i] ** 2)
        / ((a[i] - b[i]) * (A_center[i] ** 2 - A_surface[i] ** 2) * np.sqrt(1.0 + b[i]))
        * (
            GouldAuxF2(A_center[i], a[i], b[i], eta)
            - GouldAuxF2(A_surface[i], a[i], b[i], eta)
        )
        for i in range(n)
    ]

    cap_rate_sun = [T_0[i] * (T_1[i] + T_2[i] + T_3[i]) for i in range(n)]

    return np.array(cap_rate_sun) * pc.sec


def DMSunAnnihilationRateGould(
    DM_mass,
    DM_cs,
    model=None,
    n=None,
    mass_num=None,
    mass_gr_per_mol=None,
    eps=None,
    DM_rho=0.3,
    vel_rot=220.0,
):
    """Calculate annihilation rate"""
    return (
        DMSunCaptureRateGould(
            DM_mass,
            DM_cs,
            model=model,
            n=n,
            mass_num=mass_num,
            mass_gr_per_mol=mass_gr_per_mol,
            eps=eps,
            DM_rho=DM_rho,
            vel_rot=vel_rot,
        )
        / 2.0
    )


def DMSunCaptureRateHalzen(DM_mass, DM_cs, DM_rho=0.3, vel_rot=270.0):
    """Calculates DM capture rate in the Sun using Halzen-Hooper formula for mDM > 30 GeV.

    Ref : The Indirect Search for Dark Matter with IceCube. F. Halzen. D. Hooper.
    arXiv : 0910.4513

    @type  DM_mass      :      float
    @param DM_mass      :      DM mass [GeV]
    @type  DM_cs        :      float
    @param DM_cs        :      DM cross section [cm^2]
    @type  DM_rho       :      float
    @param DM_rho       :      local DM density [GeV/cm^3]
    @type  vel_rot      :      float
    @param vel_rot      :      local DM velocity[km/s]

    @rtype              :      float
    @return             :      DM capture rate at the Sun. [s^-1]
    """
    DM_mass = DM_mass * pc.GeV
    DM_cs = DM_cs * pc.cm ** 2

    DM_rho = DM_rho * pc.GeV / pc.cm ** 3  # local DM density
    vel_rot = vel_rot * pc.km / pc.sec  # local DM rel speed.

    return (
        3.35
        * 1.0e20
        * pc.sec ** -1
        * (DM_rho / (0.3 * pc.GeV / pc.cm ** 3))
        * ((270.0 * pc.km / pc.sec) / vel_rot) ** 3
        * (100.0 * pc.GeV / DM_mass) ** 2
        * (DM_cs / (1.0e-6 * pc.picobarn))
        * pc.sec
    )


def DMSunAnnihilationRateHalzen(DM_mass, DM_cs, DM_rho=0.3, vel_rot=270.0):
    """Calculate annihilation rate"""
    return DMSunCaptureRateHalzen(DM_mass, DM_cs, DM_rho=DM_rho, vel_rot=vel_rot) / 2.0


def DMCaptureJungmanSI(DM_mass, DM_cs, location="Sun", DM_rho=0.3, vel_rot=270.0):
    """Calculates DM annihilation rate in the Sun using Jungman/Kamionkowski aproximated formula
    in the case of scalar interaction.

    Ref : From Gustav Wikstrom Thesis and paper by Jungman/Kamionkowski et al.
    Eq. 9.20  pp132
    @type  DM_mass      :      float
    @param DM_mass      :      DM mass [GeV]
    @type  DM_cs        :      float
    @param DM_cs        :      DM cross section [cm^2]
    @type  DM_rho       :      float
    @param DM_rho       :      local DM density [GeV/cm^3]
    @type  vel_rot      :      float
    @param vel_rot      :      local DM velocity[km/s]


    @rtype              :      float
    @return             :      DM annihiliation rate at the Sun [annh/eV^-1].
    """

    DM_mass = DM_mass * pc.GeV
    DM_cs = DM_cs * pc.cm ** 2
    DM_rho = DM_rho * pc.GeV / pc.cm ** 3  # local DM density
    vel_rot = vel_rot * pc.km / pc.sec  # local DM rel speed.

    p_mass = 938.272 * pc.MeV  # proton mass

    if location == "Sun":
        elements = {
            0: "H",
            1: "He",
            2: "C",
            3: "N",
            4: "O",
            5: "Ne",
            6: "Mg",
            7: "Si",
            8: "S",
            9: "Fe",
        }  # names
        n = len(elements)

        mass_num = [
            1.0,
            4.0,
            12.0,
            14.0,
            16.0,
            20.0,
            24.0,
            28.0,
            32.0,
            56.0,
        ]  # A : mass number
        mass_gr_per_mol = [
            1.0079,
            4.0026,
            12.0107,
            14.0067,
            15.9994,
            20.1797,
            24.3051,
            28.0855,
            32.0655,
            55.8452,
        ]  # gr mol^-1
        eps = [
            0.772,
            0.209,
            3.87e-3,
            9.4e-4,
            8.55e-3,
            1.51e-3,
            7.39e-4,
            8.13e-4,
            4.65e-4,
            1.46e-3,
        ]  # relative aboundances in the Sun from Table 8 Jungman p.130
        phi = [
            3.16,
            3.40,
            3.23,
            3.23,
            3.23,
            3.23,
            3.23,
            3.23,
            3.23,
            3.23,
        ]  # relative gravitational potential in the Sun from Table 8 Jungman p.299
        mass = 1.9891 * 1.0e30 * pc.kg  # Sun mass
        v_escape = 1156.0 * pc.km * pc.sec ** -1
        c = 4.8e24

    elif location == "Earth":
        elements = {
            0: "O",
            1: "Si",
            2: "Mg",
            3: "Fe",
            4: "Ca",
            5: "P",
            6: "Na",
            7: "S",
            8: "Ni",
        }  # names
        n = len(elements)
        mass_num = [16.0, 28.0, 24.0, 56.0, 40.0, 30.0, 23.0, 32.0, 59.0]
        mass_gr_per_mol = [
            15.9994,
            28.0855,
            24.305,
            55.8452,
            40.078,
            30.973761,
            22.98977,
            32.0655,
            58.6934,
        ]
        eps = [0.3, 0.15, 0.14, 0.3, 0.015, 0.011, 0.004, 0.05, 0.03]
        phi = [1.2, 1.2, 1.2, 1.6, 1.2, 1.2, 1.2, 1.6, 1.6]

        mass = 5.9722e24 * pc.kg
        v_escape = 13.2 * pc.km * pc.sec ** -1
        c = 4.8e15

    mass_eV = [m * pc.gr / pc.Na for m in mass_gr_per_mol]
    atom_radius = [(1.2 * np.power(A, 1.0 / 3.0) * pc.fermi) for A in mass_num]
    energy_element = [3.0 / (2.0 * mass_eV[i] * atom_radius[i] ** 2) for i in range(n)]

    # auxiliary function
    reduce_mass = lambda m1, m2: m1 * m2 / (m1 + m2)
    # from Jungman review. ec. 9.21 and 9.22 p.299.
    A_aux = lambda x: (3.0 / 2.0) * (x / (x - 1.0) ** 2) * (v_escape / vel_rot) ** 2
    B_aux = 1.5
    kinematical_supression_factor = lambda x: (
        (A_aux(x) ** B_aux) / (1.0 + A_aux(x) ** B_aux)
    ) ** (1.0 / B_aux)

    nuclear_term = [
        (
            (
                (
                    mass_num[i]
                    * reduce_mass(DM_mass, mass_eV[i])
                    / reduce_mass(DM_mass, p_mass)
                )
                ** 2
            )
            / (mass_eV[i] / pc.GeV)
        )
        * eps[i]
        * phi[i]
        * JungmanFormFactor(i, DM_mass)
        * kinematical_supression_factor(mass_eV[i] / DM_mass)
        for i in range(n)
    ]

    C_c = (
        (c * pc.sec ** -1)
        * (DM_rho / (0.3 * pc.GeV / pc.cm ** 3))
        * ((270 * pc.km * pc.sec ** -1) / vel_rot)
        * (1.0 * pc.GeV / DM_mass)
        * (DM_cs / (1.0e-40 * pc.cm ** 2))
        * sum(nuclear_term)
    )

    return C_c * pc.sec


def DMAnnihilationJungmanSI(DM_mass, DM_cs, DM_rho=0.3, vel_rot=270.0):
    """Calculate annihilation rate"""
    return DMCaptureJungmanSC(DM_mass, DM_cs, DM_rho=DM_rho, vel_rot=vel_rot) / 2.0


def DMCaptureJungmanSD(DM_mass, DM_cs, DM_rho=0.3, vel_rot=270.0):
    """Calculates DM annihilation rate in the Sun using Jungman/Kamionkowski aproximated formula
    in the case of axial-vector interaction.

    Ref : From Gustav Wikstrom Thesis and paper by Jungman/Kamionkowski et al.
    Ec. 9.19 on p. 131

    @type  DM_mass      :      float
    @param DM_mass      :      DM mass [GeV]
    @type  DM_cs        :      float
    @param DM_cs        :      DM cross section [cm^2]
    @type  DM_rho       :      float
    @param DM_rho       :      local DM density [GeV/cm^3]
    @type  vel_rot      :      float
    @param vel_rot      :      local DM velocity[km/s]

    @rtype              :      float
    @return             :      DM capture rate at the Sun [s^-1].
    """
    DM_mass = DM_mass * pc.GeV
    DM_cs = DM_cs * pc.cm ** 2
    DM_rho = DM_rho * pc.GeV / pc.cm ** 3  # local DM density
    vel_rot = vel_rot * pc.km / pc.sec  # local DM rel speed.

    mass_gr_per_mol = [1.0079]  # H   gr mol^-1
    mass_eV = [m * pc.gr / pc.Na for m in mass_gr_per_mol]

    # sun
    v_escape = 1156.0 * pc.km * pc.sec ** -1

    # from Jungman review. ec. 9.21 and 9.22 p.299.
    A_aux = lambda x: (3.0 / 2.0) * (x / (x - 1.0) ** 2) * (v_escape / vel_rot) ** 2
    B_aux = 1.5
    kinematical_supression_factor = lambda x: (
        (A_aux(x) ** B_aux) / (1.0 + A_aux(x) ** B_aux)
    ) ** (1.0 / B_aux)

    C_c = (
        (1.3e25 * pc.sec ** -1)
        * (DM_rho / (0.3 * pc.GeV / pc.cm ** 3))
        * ((270 * pc.km * pc.sec ** -1) / vel_rot)
        * (1.0 * pc.GeV / DM_mass)
        * (DM_cs / (1.0e-40 * pc.cm ** 2))
        * kinematical_supression_factor(DM_mass / mass_eV)
    )

    return C_c * pc.sec


def DMAnnihilationJungmanSD(DM_mass, DM_cs, DM_rho=0.3, vel_rot=270.0):
    """Calculate annihilation rate"""
    return DMCaptureJungmanAX(DM_mass, DM_cs, DM_rho=DM_rho, vel_rot=vel_rot) / 2.0
