import readdy
from readdy.api.utils import load_trajectory_to_npy
import readdyviewer
import matplotlib.pyplot as plt

trajfile = '/srv/public/workspace/data/viral-capsid-hexagons/capsids.h5'

stride = 1000

n_particles_per_frame, positions, types, ids = load_trajectory_to_npy(trajfile, stride=stride)

t = readdy.Trajectory(trajfile)

config = readdyviewer.Configuration()

config.colors[t.particle_types['Head']] = readdyviewer.Color(*(plt.get_cmap("tab10")(0)[:3]))
config.radii[t.particle_types['Head']] = .5
config.colors[t.particle_types['Tail']] = readdyviewer.Color(*(plt.get_cmap("tab10")(1)[:3]))
config.radii[t.particle_types['Tail']] = .5

config.clearcolor = readdyviewer.Color(240, 240, 240)
config.smoothing = 1
config.cutoff = 2.
config.bond_radius = .3

_, topology_records = t.read_observable_topologies()

edges = []
for tops in topology_records[::stride]:
    current_edges = []
    for top in tops:
        for e1, e2 in top.edges:
            ix1 = top.particles[e1]
            ix2 = top.particles[e2]
            current_edges.append((ix1, ix2))
    edges.append(current_edges)

readdyviewer.watch_npy(positions, types, ids, n_particles_per_frame, config, edges)
