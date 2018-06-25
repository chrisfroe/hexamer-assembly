import readdy

def get_final_topologies_from_file(filepath):
    traj = readdy.Trajectory(filepath)
    final_topologies = traj.read_observable_topologies()[1][-1]
    return final_topologies

def get_final_particles_from_file(filepath):
    final_particles = readdy.Trajectory(filepath)[-1]
    return final_particles

