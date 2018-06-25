```bash
outfile=/path/to/data/capsids.h5
python run.py $outfile
./convert_traj.py $outfile traj.xyz
vmd -e traj.xyz.tcl
```
