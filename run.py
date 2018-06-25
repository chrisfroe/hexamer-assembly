import numpy as np
import os
import argparse

import capsid_system

parser = argparse.ArgumentParser()
parser.add_argument('traj_file_path', type=str)


def generate_3d_rotation():
    """Generate a 3D random rotation matrix.
    Returns:
        np.matrix: A 3D rotation matrix.

    See https://github.com/qobilidop/randrot/blob/master/randrot/__init__.py
    """
    x1, x2, x3 = np.random.rand(3)
    R = np.matrix([[np.cos(2 * np.pi * x1), np.sin(2 * np.pi * x1), 0],
                   [-np.sin(2 * np.pi * x1), np.cos(2 * np.pi * x1), 0],
                   [0, 0, 1]])
    v = np.matrix([[np.cos(2 * np.pi * x2) * np.sqrt(x3)],
                   [np.sin(2 * np.pi * x2) * np.sqrt(x3)],
                   [np.sqrt(1 - x3)]])
    H = np.eye(3) - 2 * v * v.T
    M = -H * R
    return M


if __name__ == '__main__':
    args = parser.parse_args()
    traj_file_path = args.traj_file_path

    boxsize = 25.
    system = capsid_system.configure_capsid_system(box_size=[boxsize, boxsize, boxsize], unit_system=None)

    simulation = system.simulation(kernel="SingleCPU")

    # generate initial positions for particles
    number_ca = 100
    for i in range(number_ca):
        core = np.array([0., 0., 0.])
        site1 = np.array([0., 0., 1.])
        site2 = np.array([np.sin(np.pi * 60. / 180.), 0., - 1. * np.cos(np.pi * 60. / 180.)])
        rot = generate_3d_rotation()

        site1 = np.dot(rot, site1)
        site2 = np.dot(rot, site2)

        site1 = np.squeeze(np.asarray(site1))
        site2 = np.squeeze(np.asarray(site2))

        origin = np.array([-boxsize / 2., -boxsize / 2., -boxsize / 2.]) + 2.
        extent = np.array([boxsize, boxsize, boxsize]) - 4.

        translation = np.random.uniform(size=3) * extent + origin
        core += translation
        site1 += translation
        site2 += translation

        top = simulation.add_topology("CA", ["site", "core", "site"], np.array([site1, core, site2]))
        top.get_graph().add_edge(0, 1)
        top.get_graph().add_edge(1, 2)

    if os.path.exists(traj_file_path):
        os.remove(traj_file_path)
    simulation.output_file = traj_file_path
    simulation.observe.topologies(100)
    simulation.record_trajectory(100)
    simulation.progress_output_stride = 100
    simulation.show_progress = True

    simulation.run(500000, .01)
