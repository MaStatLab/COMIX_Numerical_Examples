import os
from glob import glob
import flowkit as fk

# Convert fcs to csv:

base_dir = os.getcwd()

fcs_dir = os.path.join(base_dir, 'fcs')
xform_dir = os.path.join(base_dir, 'csv_xform')

fcs_paths = glob(os.path.join(fcs_dir, '*.fcs'))

xform = fk.transforms.AsinhTransform(
    'my_xform',
    param_t = 12000,
    param_m = 4.0,
    param_a = 0.7
)

for f in fcs_paths:
    sample = fk.Sample(f)
    sample.apply_transform(xform)
    
    new_name = os.path.basename(f).replace('fcs', 'csv')

    sample.export_csv(source = 'xform', filename = new_name, directory = xform_dir)

# Extract Batch Control:
from matplotlib.patches import Rectangle
from matplotlib.path import Path
import pandas as pd
import numpy as np

fs = glob('./csv_xform/*')

for f in fs:
    df = pd.read_csv(f)
    x = df['Pr141Di']
    y = df['Y89Di']
    r1 = Rectangle((0.1, 0.5), 0.4, 0.3, linewidth = 2, edgecolor = 'red', facecolor = 'none')
    p1 = Path(r1.get_verts())
    bc = df[p1.contains_points(np.c_[x, y])]
    fn = f.split('/')[-1]
    
    bc.to_csv(os.path.join('csv_xform', fn.replace('.csv', '_bc.csv')))
