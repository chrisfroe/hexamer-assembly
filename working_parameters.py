import numpy as np
import readdy

def configure_lipid_system(**kwargs):
    """System of capsid proteins forming hexagons."""
    if "unit_system" not in kwargs:
        kwargs.update({"unit_system": None})
    if "box_size" not in kwargs:
        kwargs.update({"box_size": [25, 25, 25]})
    system = readdy.ReactionDiffusionSystem(**kwargs)

    system.topologies.add_type("CA")
    system.add_topology_species("core", .1)
    system.add_topology_species("site", .1)

    system.topologies.configure_harmonic_bond("core", "site", force_constant=60, length=1.)
    system.topologies.configure_harmonic_bond("core", "core", force_constant=60, length=2.)
    system.topologies.configure_harmonic_bond("site", "site", force_constant=60, length=0.1)
    system.topologies.configure_harmonic_angle("site", "core", "site", force_constant=10., equilibrium_angle=120./180. * np.pi)
    system.topologies.configure_harmonic_angle("site", "core", "core", force_constant=10., equilibrium_angle=120./180. * np.pi)
    system.topologies.configure_harmonic_angle("core", "core", "core", force_constant=10., equilibrium_angle=120./180. * np.pi)

    system.potentials.add_harmonic_repulsion("core", "core", force_constant=80., interaction_distance=2.)
    system.potentials.add_weak_interaction_piecewise_harmonic("Tail", "Tail", force_constant=800.,                                                              desired_distance=1.122462048309373,                                                              depth=0.91 * system.kbt,                                                              cutoff=2.722462048309373)
    return system
