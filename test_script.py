#!/usr/bin/env python3
import numpy as np
import os, sys, re
import mantid.simpleapi as s_api
currdir = os.path.join(os.path.dirname(__file__))
sys.path.append(currdir)
# Loads the reduction scripts
from reduction_utils import iliad

# Tests a saved livestream dataset from MnO on MERLIN
datapath = os.path.join(currdir, 'datarepo', 'datafiles')
s_api.config.appendDataSearchDir(datapath)
s_api.config.appendDataSearchDir(os.path.join(currdir, 'InstrumentFiles', 'merlin'))
ws = s_api.Load('mno_live_1.nxs', OutputWorkspace='mno_live_1')

# Parameters for MnO
inst = 'merlin'                  # Instrument
maskfile = 'mask_25_1.xml'       # Mask file name
ei = 30                          # Incident energy of run in meV
vanrun = 71617                   # Run number of vanadium calibration run
lattice = [4.446, 4.446, 4.446]  # Lattice parameters in Angstrom
lat_ang = [90., 90., 90.]        # Lattice angles in degrees
uvec = [1, 1, 0]                 # Vector of incident beam when psi=0
vvec = [0, 0, 1]                 # A vector perpendicular to incident beam when psi=0
psi0 = 132.5                     # Logged rotation angle when psi=0
blockname = 'Rot'                # Name of the logged rotation "block"
ang_bin_size = 5.                # Angle bin size in degrees (really should be 0.5 or 1.0)
ang_tol = 0.001                  # Tolerance on angle bin size
ignorezeros = False              # Whether to keep or ignore zeros in data

# Current workflow - bin in angles, run reduction (iliad) at each angle then use MergeMD
logval = ws.getRun().getLogData(blockname).filtered_value
ang_min, ang_max = (min(logval), max(logval))
ang_range = ang_max - ang_min
ang_nbin = int(ang_range / ang_bin_size)
ang_ws = []
for irot, ang_start in enumerate([x * ang_bin_size + ang_min for x in range(ang_nbin)]):
    ang_end = min(ang_start + ang_bin_size - ang_tol, ang_max)
    ws_ang = s_api.FilterByLogValue(ws, blockname, ang_start, ang_end)
    iliad(runno=['ws_ang'], wbvan=vanrun, ei=ei, FixEi=True, inst=inst, hard_mask_file=maskfile,
          powder=False, cs_single=True) # This produces a workspace called 'single_md'
    s_api.SetGoniometer('single_md', Axis0=f'{blockname},0,1,0,1', Axis1=f'{psi0},0,1,0,-1')
    s_api.SetUB('single_md', a=lattice[0], b=lattice[1], c=lattice[2], 
                alpha=lat_ang[0], beta=lat_ang[1], gamma=lat_ang[2],
                u=f'{uvec[0]},{uvec[1]},{uvec[2]}', v=f'{vvec[0]},{vvec[1]},{vvec[2]}')
    s_api.ConvertToMD('single_md', QDimensions='Q3D', dEAnalysisMode='Direct', Q3DFrames='HKL',
                      QConversionScales='HKL', PreprocDetectorsWS='-',
                      OutputWorkspace=f'ang_md_{irot}', IgnoreZeroSignals=ignorezeros)
    ang_ws.append(f'ang_md_{irot}')
# Now merge the workspaces
wsout = s_api.MergeMD(','.join(ang_ws), OutputWorkspace=f'ws_out_{ei}meV_1to1_md')
for aws in ang_ws:
    s_api.DeleteWorkspace(aws)

"""
# Future workflow without angular binning (Feature #1)
iliad(runno=['ws_ang'], wbvan=vanrun, ei=ei, FixEi=True, inst=inst, hard_mask_file=maskfile,
      powder=False, cs_single=True)
s_api.SetGoniometer('single_md', Axis0=f'{blockname},0,1,0,1', Axis1=f'{psi0},0,1,0,-1')
s_api.SetUB('single_md', a=lattice[0], b=lattice[1], c=lattice[2], alpha=lat_ang[0], beta=lat_ang[1], gamma=lat_ang[2],
            u=f'{uvec[0]},{uvec[1]},{uvec[2]}', v=f'{vvec[0]},{vvec[1]},{vvec[2]}')
mn, mx = s_api.ConvertToMDMinMaxGlobal('single_md', QDimensions='Q3D', dEAnalysisMode='Direct', Q3DFrames='HKL')
ws_out = s_api.ConvertToMD('single_md', QDimensions='Q3D', dEAnalysisMode='Direct', Q3DFrames='HKL',
                           QConversionScales='HKL', MinValues=mn, MaxValues=mx, PreprocDetectorsWS='-',
                           IgnoreZeroSignals=ignorezeros)
"""

# And histogram it to plot
mn = [wsout.getDimension(i).getMinimum() for i in range (4)]
mx = [wsout.getDimension(i).getMaximum() for i in range (4)]
wsbin = s_api.BinMD(wsout, AlignedDim0=f'[H,0,0],{mn[0]},{mx[0]},10', AlignedDim1=f'[0,K,0],{mn[1]},{mx[1]},10',
    AlignedDim2=f'[0,0,L],{mn[2]},{mx[2]},10', AlignedDim3=f'DeltaE,{mn[3]},{mx[3]},10')
wsbin = s_api.CompactMD(wsbin)
s_api.SaveMD(wsbin, os.path.join(currdir, 'ws_out_1to1_md.nxs'))

# wsbin can be plotted with SliceViewer.
