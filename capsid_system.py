import numpy as np
import readdy


def configure_capsid_system(**kwargs):
    """System of capsid proteins forming hexagons"""
    if "unit_system" not in kwargs:
        kwargs.update({"unit_system": None})
    if "box_size" not in kwargs:
        kwargs.update({"box_size": [25, 25, 25]})
    system = readdy.ReactionDiffusionSystem(**kwargs)

    system.topologies.add_type("CA")
    system.add_topology_species("core", 0.1)
    system.add_topology_species("site", 0.1)

    system.add_species("dummy", 0.1)  # instantly decays
    system.reactions.add("dummydecay: dummy ->", rate=1e12)

    system.topologies.configure_harmonic_bond("core", "site", force_constant=100, length=1.)
    system.topologies.configure_harmonic_bond("core", "core", force_constant=100, length=2.)
    system.topologies.configure_harmonic_bond("site", "site", force_constant=100, length=0.1)
    system.topologies.configure_harmonic_angle("site", "core", "site", force_constant=100.,
                                               equilibrium_angle=120. / 180. * np.pi)
    system.topologies.configure_harmonic_angle("site", "core", "core", force_constant=100.,
                                               equilibrium_angle=120. / 180. * np.pi)
    system.topologies.configure_harmonic_angle("core", "core", "core", force_constant=150.,
                                               equilibrium_angle=120. / 180. * np.pi)
    dihedral_angle = 0.
    system.topologies.configure_cosine_dihedral("core", "core", "core", "core", 100., 1., dihedral_angle)
    system.topologies.configure_cosine_dihedral("site", "core", "core", "core", 100., 1., dihedral_angle)
    system.topologies.configure_cosine_dihedral("site", "core", "core", "site", 100., 1., dihedral_angle)

    system.potentials.add_harmonic_repulsion("core", "core", force_constant=80., interaction_distance=2.)

    system.topologies.add_spatial_reaction("attach: CA(site)+CA(site)->CA(site--site) [self=true]", rate=10.,
                                           radius=0.5)

    def clean_sites_rate_function(topology):
        edges = topology.get_graph().get_edges()
        vertices = topology.get_graph().get_vertices()

        if len(vertices) > 3:
            for e in edges:
                v1_ref, v2_ref = e[0], e[1]
                v1 = v1_ref.get()
                v2 = v2_ref.get()
                v1_type = topology.particle_type_of_vertex(v1)
                v2_type = topology.particle_type_of_vertex(v2)
                if v1_type == "site" and v2_type == "site":
                    return 1e12
        else:
            return 0.
        return 0.

    def clean_sites_reaction_function(topology):

        recipe = readdy.StructuralReactionRecipe(topology)
        vertices = topology.get_graph().get_vertices()

        def search_configuration():
            # dfs for finding configuration core-site-site-core
            for v1 in vertices:
                if topology.particle_type_of_vertex(v1) == "core":
                    for v2_ref in v1.neighbors():
                        v2 = v2_ref.get()
                        if topology.particle_type_of_vertex(v2) == "site":
                            for v3_ref in v2.neighbors():
                                v3 = v3_ref.get()
                                if v3.particle_index != v1.particle_index:
                                    if topology.particle_type_of_vertex(v3) == "site":
                                        for v4_ref in v3.neighbors():
                                            v4 = v4_ref.get()
                                            if v4.particle_index != v2.particle_index:
                                                if topology.particle_type_of_vertex(v4) == "core":
                                                    return v1.particle_index, v2.particle_index, v3.particle_index, v4.particle_index

        core1_p_idx, site1_p_idx, site2_p_idx, core2_p_idx = search_configuration()

        # find corresponding vertex indices from particle indices
        core1_v_idx = None
        site1_v_idx = None
        site2_v_idx = None
        core2_v_idx = None
        for i, v in enumerate(vertices):
            if v.particle_index == core1_p_idx and core1_v_idx is None:
                core1_v_idx = i
            elif v.particle_index == site1_p_idx and site1_v_idx is None:
                site1_v_idx = i
            elif v.particle_index == site2_p_idx and site2_v_idx is None:
                site2_v_idx = i
            elif v.particle_index == core2_p_idx and core2_v_idx is None:
                core2_v_idx = i
            else:
                pass

        if (core1_v_idx is not None) and (core2_v_idx is not None) and (site1_v_idx is not None) and (
                site2_v_idx is not None):
            recipe.add_edge(core1_v_idx, core2_v_idx)
            recipe.separate_vertex(site1_v_idx)
            recipe.separate_vertex(site2_v_idx)
            recipe.change_particle_type(site1_v_idx, "dummy")
            recipe.change_particle_type(site2_v_idx, "dummy")
        else:
            raise RuntimeError("core-site-site-core wasn't found")

        return recipe

    system.topologies.add_structural_reaction(topology_type="CA", reaction_function=clean_sites_reaction_function,
                                              rate_function=clean_sites_rate_function, raise_if_invalid=True,
                                              expect_connected=False)

    return system
