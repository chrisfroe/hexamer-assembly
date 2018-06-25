#!/usr/bin/env python
import argparse
import os
import readdy

parser = argparse.ArgumentParser()
parser.add_argument('traj_file_path', type=str)
parser.add_argument('out_name', type=str)

if __name__ == '__main__':
    print("run convert_traj ..")

    args = parser.parse_args()
    traj_file_path = args.traj_file_path
    out_name = args.out_name
    
    traj = readdy.Trajectory(traj_file_path)
    particle_radii = {"Head": 0.5, "Tail": 0.5}
    traj.convert_to_xyz(xyz_filename=out_name, particle_radii=particle_radii)

    print("done convert_traj ..")

