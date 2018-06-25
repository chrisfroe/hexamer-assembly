import numpy as np
import readdy
import os

import working_parameters

bs = 25.
system = working_parameters.configure_capsid_system(box_size=[bs, bs, bs], unit_system=None)

simulation = system.simulation(kernel="SingleCPU")

top = simulation.add_topology("CA", ["site", "core", "site"],
                                      np.array([head, tail1, tail2]))
top.get_graph().add_edge(0, 1)
top.get_graph().add_edge(1, 2)

# data_dir = "/home/chris/workspace/data/"
data_dir = "/srv/public/workspace/data/viral-capsid-hexagons/"
simulation.output_file = os.path.join(data_dir, "capsids.h5")
if os.path.exists(simulation.output_file):
    os.remove(simulation.output_file)
simulation.observe.topologies(10)
simulation.record_trajectory(10)
simulation.progress_output_stride = 10
simulation.show_progress = True

simulation.run(4000, .001)

