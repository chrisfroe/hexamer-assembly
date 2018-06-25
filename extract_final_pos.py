#!/usr/bin/env python
import numpy as np
import argparse
import readdy
from readdy.api.utils import load_trajectory_to_npy

parser = argparse.ArgumentParser()
parser.add_argument('traj_file_path', type=str)
parser.add_argument('position_file_path', type=str)

if __name__ == '__main__':
    args = parser.parse_args()
    traj_file_path = args.traj_file_path
    position_file_path = args.position_file_path

    stride = 100

    n_particles_per_frame, positions, types, ids = load_trajectory_to_npy(traj_file_path, stride=stride)

    traj = readdy.Trajectory(traj_file_path)

    # @todo generic topology
    _, topology_records = traj.read_observable_topologies()

    final_tops = topology_records[-1]
    final_pos = positions[-1]

    lipid_pos = np.zeros((len(final_tops), 3, 3))
    for ix, top in enumerate(final_tops):
        head_index = top.particles[0]
        tail1_index = top.particles[1]
        tail2_index = top.particles[2]
        lipid_pos[ix, 0] = final_pos[head_index]
        lipid_pos[ix, 1] = final_pos[tail1_index]
        lipid_pos[ix, 2] = final_pos[tail2_index]

    np.save(position_file_path, lipid_pos)
