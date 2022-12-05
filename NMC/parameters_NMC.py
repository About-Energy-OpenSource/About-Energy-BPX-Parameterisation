# 12.5 Ah NMC111|Graphite Pouch Cell - PyBaMM DFN Model Parameters for BPX
#
# About:Energy Limited, UK
#
# https://www.aboutenergy.io
#


from pybamm import exp, tanh, constants, Parameter, ParameterValues, sqrt
from numpy import array, nan


def electrolyte_diffusivity_Nyman2008(c_e, T):
    """
    Diffusivity of LiPF6 in EC:EMC (3:7) as a function of ion concentration. The data
    comes from [1]

    References
    ----------
    .. [1] A. Nyman, M. Behm, and G. Lindbergh, "Electrochemical characterisation and
    modelling of the mass transport phenomena in LiPF6-EC-EMC electrolyte,"
    Electrochim. Acta, vol. 53, no. 22, pp. 6356–6365, 2008.

    Parameters
    ----------
    c_e: :class:`pybamm.Symbol`
        Dimensional electrolyte concentration
    T: :class:`pybamm.Symbol`
        Dimensional temperature

    Returns
    -------
    :class:`pybamm.Symbol`
        Electrolyte diffusivity
    """

    D_c_e = 8.794e-11 * (c_e / 1000) ** 2 - 3.972e-10 * (c_e / 1000) + 4.862e-10

    # Nyman et al. (2008) does not provide temperature dependence
    # Add an arrhenius fit
    Eact = 17100
    Tref = Parameter("Reference temperature [K]")
    arrhenius = exp(Eact / constants.R * (1 / Tref - 1 / T))

    return D_c_e * arrhenius


def electrolyte_conductivity_Nyman2008(c_e, T):
    """
    Conductivity of LiPF6 in EC:EMC (3:7) as a function of ion concentration. The data
    comes from [1].

    References
    ----------
    .. [1] A. Nyman, M. Behm, and G. Lindbergh, "Electrochemical characterisation and
    modelling of the mass transport phenomena in LiPF6-EC-EMC electrolyte,"
    Electrochim. Acta, vol. 53, no. 22, pp. 6356–6365, 2008.

    Parameters
    ----------
    c_e: :class:`pybamm.Symbol`
        Dimensional electrolyte concentration
    T: :class:`pybamm.Symbol`
        Dimensional temperature

    Returns
    -------
    :class:`pybamm.Symbol`
        Electrolyte conductivity
    """

    sigma_e = (
        0.1297 * (c_e / 1000) ** 3 - 2.51 * (c_e / 1000) ** 1.5 + 3.329 * (c_e / 1000)
    )

    # Nyman et al. (2008) does not provide temperature dependence
    # Add an arrhenius fit
    Eact = 17100
    Tref = Parameter("Reference temperature [K]")
    arrhenius = exp(Eact / constants.R * (1 / Tref - 1 / T))

    return sigma_e * arrhenius


def ocp_n(sto):
    """
    Graphite open circuit potential as a function of stochiometry.

    Parameters
    ----------
    sto: :class:`pybamm.Symbol`
        Electrode stochiometry

    Returns
    -------
    :class:`pybamm.Symbol`
        Open circuit potential
    """

    p = array(
        [
            9.47057878e-01,  1.59418743e+02, -3.50928033e+04,  1.64230269e-01,
            -4.55509094e+01,  3.24116012e-02,  3.69968491e-02, -1.96718868e+01,
            1.68334476e-01,  1.91517003e+04,  3.19648312e+00,  1.85139824e+00,
            5.42448511e+04, -3.19009848e+00,  2.01660395e+00
        ]
    )

    u_eq = (
        # fit a linear combination of exponentials
        + p[0] * exp(-p[1] * sto)
        + p[2]
        + p[3] * tanh(p[4] * (sto - p[5]))
        + p[6] * tanh(p[7] * (sto - p[8]))
        + p[9] * tanh(p[10] * (sto - p[11]))
        + p[12] * tanh(p[13] * (sto - p[14]))          
    )
    return u_eq



def ocp_p(sto):
    """
    NMC open circuit potential as a function of stochiometry.

    Parameters
    ----------
    sto: :class:`pybamm.Symbol`
        Electrode stochiometry

    Returns
    -------
    :class:`pybamm.Symbol`
        Open circuit potential
    """

    p = array(
        [   
            -3.04420906, 10.04892207, -0.65637536, -4.02134095,  0.80063948,
            4.24678547, 12.17805062,  7.57659337, -0.3757068 , 59.33067782,
            0.99784492
        ]
    )

    u_eq = (
        + p[0] * sto
        + p[1]
        + p[2] * tanh(p[3] * (sto - p[4]))
        + p[5] * tanh(p[6] * (sto - p[7]))
        + p[8] * tanh(p[9] * (sto - p[10]))
    )

    return u_eq


def negative_electrode_exchange_current_density(c_e, c_s_surf, c_s_max, T):
    """
    Exchange-current density for Butler-Volmer reactions.
    Parameters
    ----------
    c_e : :class:`pybamm.Symbol`
        Electrolyte concentration [mol.m-3]
    c_s_surf : :class:`pybamm.Symbol`
        Particle concentration [mol.m-3]
    c_s_max : :class:`pybamm.Symbol`
        Maximum particle concentration [mol.m-3]
    T : :class:`pybamm.Symbol`
        Temperature [K]
    Returns
    -------
    :class:`pybamm.Symbol`
        Exchange-current density [A.m-2]
    """

    j0_50 = 0.250805067133009
    c_e_typ = Parameter("Initial concentration in electrolyte [mol.m-3]")

    k_ref = 2*j0_50 / c_s_max / sqrt(c_e_typ) # (A/m2)(mol/m3)**1.5 - includes ref concentrations

    Eact = 55000
    Tref = Parameter("Reference temperature [K]")
    arrhenius = exp(Eact / constants.R * (1 / Tref - 1 / T))

    return (
        k_ref * arrhenius * c_e ** 0.5 * c_s_surf ** 0.5 * (c_s_max - c_s_surf) ** 0.5
    )

def positive_electrode_exchange_current_density(c_e, c_s_surf, c_s_max, T):
    """
    Exchange-current density for Butler-Volmer reactions.
    Parameters
    ----------
    c_e : :class:`pybamm.Symbol`
        Electrolyte concentration [mol.m-3]
    c_s_surf : :class:`pybamm.Symbol`
        Particle concentration [mol.m-3]
    c_s_max : :class:`pybamm.Symbol`
        Maximum particle concentration [mol.m-3]
    T : :class:`pybamm.Symbol`
        Temperature [K]
    Returns
    -------
    :class:`pybamm.Symbol`
        Exchange-current density [A.m-2]
    """

    j0_50 = 1.11178673498325
    c_e_typ = Parameter("Initial concentration in electrolyte [mol.m-3]")

    k_ref = 2*j0_50 / c_s_max / sqrt(c_e_typ) # (A/m2)(mol/m3)**1.5 - includes ref concentrations
    
    Eact = 35000
    Tref = Parameter("Reference temperature [K]")
    arrhenius = exp(Eact / constants.R * (1 / Tref - 1 / T))

    return (
        k_ref * arrhenius * c_e ** 0.5 * c_s_surf ** 0.5 * (c_s_max - c_s_surf) ** 0.5
    )

def graphite_entropic_change_ORegan2021(sto, c_s_max):
    """
    Graphite entropic change in open circuit potential (OCP) at a temperature of
    298.15K as a function of the stochiometry. The fit is taken from [1].
    References
    ----------
    .. [1] Kieran O’Regan, Ferran Brosa Planella, W. Dhammika Widanage, and Emma
    Kendrick. "Thermal-electrochemical parametrisation of a lithium-ion battery:
    mapping Li concentration and temperature dependencies." Journal of the
    Electrochemical Society, submitted (2021).
    Parameters
    ----------
    sto: :class:`pybamm.Symbol`
       Electrode stochiometry
    Returns
    -------
    :class:`pybamm.Symbol`
       Entropic change [V.K-1]
    """

    a0 = -0.1112
    a1 = -0.09002 * 0  # fixed fit (see discussion O'Regan et al 2021)
    a2 = 0.3561
    b1 = 0.4955
    b2 = 0.08309
    c0 = 0.02914
    c1 = 0.1122
    c2 = 0.004616
    d1 = 63.9

    dUdT = (
        a0 * sto
        + c0
        + a2 * exp(-((sto - b2) ** 2) / c2)
        + a1 * (tanh(d1 * (sto - (b1 - c1))) - tanh(d1 * (sto - (b1 + c1))))
    ) / 1000  # fit in mV / K

    return dUdT

def graphite_diffusivity(sto, T):
    """
    Graphite diffusivity.

    Parameters
    ----------
    sto: :class:`pybamm.Symbol`
        Electrode stochiometry
    T: :class:`pybamm.Symbol`
        Dimensional temperature

    Returns
    -------
    :class:`pybamm.Symbol`
        Solid diffusivity
    """
    D_ref = 2.728e-14
    Eact = 30000
    arrhenius = exp(Eact / constants.R * (1 / 298.15 - 1 / T))

    return D_ref * arrhenius


def nmc_diffusivity(sto, T):
    """
    NMC diffusivity.

    Parameters
    ----------
    sto: :class:`pybamm.Symbol`
      Electrode stochiometry
    T: :class:`pybamm.Symbol`
       Dimensional temperature

    Returns
    -------
    :class:`pybamm.Symbol`
       Solid diffusivity
    """
    D_ref = 3.2e-14

    Eact = 15000
    arrhenius = exp(Eact / constants.R * (1 / 298.15 - 1 / T))

    return D_ref * arrhenius


parameter_values = ParameterValues(
    {   
        '1 + dlnf/dlnc': 1.0,
        'Ambient temperature [K]': 298.15,
        'Cation transference number': 0.2594,
        'Cell cooling surface area [m2]': 0.0378568,
        'Cell volume [m3]': 0.00012774,
        'Current function [A]': 12.5,
        'Electrode height [m]': 0.191,
        'Electrode width [m]': 0.088,
        "Electrolyte conductivity [S.m-1]": electrolyte_conductivity_Nyman2008,
        "Electrolyte diffusivity [m2.s-1]": electrolyte_diffusivity_Nyman2008,
        'Initial concentration in electrolyte [mol.m-3]': 1000.0,
        'Initial concentration in negative electrode [mol.m-3]': 22496.0,
        'Initial concentration in positive electrode [mol.m-3]': 19600.0,
        'Initial temperature [K]': 298.15,
        'Lower voltage cut-off [V]': 2.7,
        'Maximum concentration in negative electrode [mol.m-3]': 29730.0,
        'Maximum concentration in positive electrode [mol.m-3]': 46200.0,
        'Negative current collector conductivity [S.m-1]': 58700000.0,
        'Negative current collector density [kg.m-3]': 8933.0,
        'Negative current collector specific heat capacity [J.kg-1.K-1]': 383.0,
        'Negative current collector thermal conductivity [W.m-1.K-1]': 401.0,
        'Negative current collector thickness [m]': 8e-06,
        'Negative electrode active material volume fraction': 0.68601,
        'Negative electrode Bruggeman coefficient (electrode)': 1.5,
        'Negative electrode Bruggeman coefficient (electrolyte)': 1.5,
        'Negative electrode conductivity [S.m-1]': 0.39,
        'Negative electrode density [kg.m-3]': 1454.2,
        'Negative electrode diffusivity [m2.s-1]': graphite_diffusivity,
        'Negative electrode electrons in reaction': 1.0,
        'Negative electrode exchange-current density [A.m-2]': negative_electrode_exchange_current_density,
        'Negative electrode OCP [V]': ocp_n,
        'Negative electrode OCP entropic change [V.K-1]': graphite_entropic_change_ORegan2021,
        'Negative electrode porosity': 0.253991,
        'Negative electrode specific heat capacity [J.kg-1.K-1]': 846.0,
        'Negative electrode thermal conductivity [W.m-1.K-1]': 3.793,
        'Negative electrode thickness [m]': 5.62e-05,
        'Negative particle radius [m]': 4.12e-06,
        'Nominal cell capacity [A.h]': 12.5,
        'Number of cells connected in series to make a battery': 1.0,
        'Number of electrodes connected in parallel to make a cell': 34.0,
        'Positive current collector conductivity [S.m-1]': 36900000.0,
        'Positive current collector density [kg.m-3]': 2700.0,
        'Positive current collector specific heat capacity [J.kg-1.K-1]': 876.0,
        'Positive current collector thermal conductivity [W.m-1.K-1]': 237.0,
        'Positive current collector thickness [m]': 1.5e-05,
        'Positive electrode active material volume fraction': 0.66251,
        'Positive electrode Bruggeman coefficient (electrode)': 1.5,
        'Positive electrode Bruggeman coefficient (electrolyte)': 1.5,
        'Positive electrode conductivity [S.m-1]': 1.464,
        'Positive electrode density [kg.m-3]': 2383.6,
        'Positive electrode diffusivity [m2.s-1]': nmc_diffusivity,
        'Positive electrode electrons in reaction': 1.0,
        'Positive electrode exchange-current density [A.m-2]': positive_electrode_exchange_current_density,
        'Positive electrode OCP [V]': ocp_p,
        'Positive electrode OCP entropic change [V.K-1]': -1e-4,
        'Positive electrode porosity': 0.277493,
        'Positive electrode specific heat capacity [J.kg-1.K-1]': 902.0,
        'Positive electrode thermal conductivity [W.m-1.K-1]': 0.807,
        'Positive electrode thickness [m]': 5.23e-05,
        'Positive particle radius [m]': 4.6e-06,
        'Reference temperature [K]': 298.15,
        'Separator Bruggeman coefficient (electrolyte)': 1.5,
        'Separator density [kg.m-3]': 1548.0,
        'Separator porosity': 0.47,
        'Separator specific heat capacity [J.kg-1.K-1]': 1128.0,
        'Separator thermal conductivity [W.m-1.K-1]': 0.3344,
        'Separator thickness [m]': 2e-05,
        'Total heat transfer coefficient [W.m-2.K-1]': 10.0,
        'Typical current [A]': 12.5,
        'Typical electrolyte concentration [mol.m-3]': 1000.0,
        'Upper voltage cut-off [V]': 4.2
    }
)