#read a ply file point cloud and output as a simple XYZ file

import sys
import numpy as np
import plyfile
import time

def main():
    # if len(sys.argv) != 3:
    #     print("usage: ply2xyz.py <input.ply> <output.xyz>")
    #     sys.exit(1)

    #input = 'carbon_plastic_items_quality.ply'
    input = 'MultiShiny-Carbon-Original.ply'
    #name of the system data and time
    output = time.strftime("%Y_%m_%d-%H-%M-%S") + '.xyz'

    plydata = plyfile.PlyData.read(input)

    xyz = np.vstack([plydata['vertex'][x] for x in ['x', 'y', 'z']]).T

    np.savetxt(output, xyz, fmt='%.6f')

if __name__ == '__main__':
    main()
    
