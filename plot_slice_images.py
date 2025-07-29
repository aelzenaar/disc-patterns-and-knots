import pandas as pd
import numpy as np
import holoviews as hv
from holoviews import opts
from holoviews.operation.datashader import rasterize,regrid
hv.extension('matplotlib')
hv.output(fig='png', dpi=600)
fs = {'labels': 8, 'ticks': 6}


overlay = hv.Points([(1,0)]).opts(color="white",marker="o",s=50)*hv.Text(.4,-.3,'Whitehead\ncusp', fontsize=8).opts(color="white")
data = pd.read_csv('wielenberg_slice.csv',names=['Re','Im','minimal tr.len']).sort_values('minimal tr.len', ascending=True).drop_duplicates(subset=['Re','Im'],keep='first').set_index(['Re', 'Im'])
array = data.to_xarray()
ds = hv.Dataset(array,['Re', 'Im'], 'minimal tr.len')
img = ds.to(hv.Image, ['Re', 'Im'], 'minimal tr.len')
img = regrid(img, upsample=True, interpolation='bilinear',y_sampling=.01, x_sampling=.01, aggregator='min').opts(cmap='gray', aspect=1, fig_inches=4, colorbar=True,xlim=(-2,2),ylim=(-2,2), fontsize=fs)
hv.save(img*overlay, f'wielenberg_slice_horiz.png', fmt='png')


data = pd.read_csv('schottky_slice_easy.csv',names=['Re','Im','minimal tr.len']).drop_duplicates(subset=['Re','Im'],keep='last').set_index(['Re', 'Im'])
array = data.to_xarray()
ds = hv.Dataset(array,['Re', 'Im'], 'minimal tr.len')
img = ds.to(hv.Image, ['Re', 'Im'], 'minimal tr.len')
img = regrid(img, upsample=True, interpolation='bilinear',y_sampling=.01, x_sampling=.01, aggregator='min').opts(cmap='gray', aspect=1, fig_inches=4, colorbar=True,xlim=(-2,2),ylim=(-2,2), fontsize=fs)
hv.save(img, f'schottky_slice_easy.png', fmt='png')

overlay = hv.Points([(1,0)]).opts(color="white",marker="o",s=50)*hv.Text(1,-.15,'cusp', fontsize=8).opts(color="white")*hv.Points([(0,0)]).opts(color="white",marker="o",s=50)*hv.Text(0,.15,'orbifold', fontsize=8).opts(color="white")
data = pd.read_csv('schottky_slice_hard.csv',names=['Re','Im','minimal tr.len']).drop_duplicates(subset=['Re','Im'],keep='last').set_index(['Re', 'Im'])
array = data.to_xarray()
ds = hv.Dataset(array,['Re', 'Im'], 'minimal tr.len')
img = ds.to(hv.Image, ['Re', 'Im'], 'minimal tr.len')
img = regrid(img, upsample=True, interpolation='bilinear',y_sampling=.01, x_sampling=.01, aggregator='min').opts(cmap='gray', aspect=1, fig_inches=4, colorbar=True,xlim=(-2,2),ylim=(-2,2), fontsize=fs)
hv.save(img*(overlay), f'schottky_slice_hard.png', fmt='png')

overlay = hv.Points([(1,0)]).opts(color="white",marker="o",s=50)*hv.Text(1,.15,'cusp', fontsize=8).opts(color="white")*hv.Points([(0,0)]).opts(color="white",marker="o",s=50)*hv.Text(0,.15,'knot', fontsize=8).opts(color="white")
data = pd.read_csv('whitehead_cusp.csv',names=['Re','Im','minimal tr.len']).drop_duplicates(subset=['Re','Im'],keep='last').set_index(['Re', 'Im'])
array = data.to_xarray()
ds = hv.Dataset(array,['Re', 'Im'], 'minimal tr.len')
img = ds.to(hv.Image, ['Re', 'Im'], 'minimal tr.len')
img = regrid(img, upsample=True, interpolation='bilinear',y_sampling=.01, x_sampling=.01, aggregator='min').opts(cmap='gray', aspect=1, fig_inches=4, colorbar=True,xlim=(-1,3),ylim=(-2,2), fontsize=fs)
hv.save(img*(overlay), f'whitehead_cusp.png', fmt='png')

overlay = hv.Points([(1,0)]).opts(color="white",marker="o",s=50)*hv.Text(1,.15,'cusp', fontsize=8).opts(color="white")*hv.Points([(0,0)]).opts(color="white",marker="o",s=50)*hv.Text(0,.30,'Fuchsian\ngroup', fontsize=8).opts(color="white")
data = pd.read_csv('solomon_cusp.csv',names=['Re','Im','minimal tr.len']).drop_duplicates(subset=['Re','Im'],keep='last').set_index(['Re', 'Im'])
array = data.to_xarray()
ds = hv.Dataset(array,['Re', 'Im'], 'minimal tr.len')
img = ds.to(hv.Image, ['Re', 'Im'], 'minimal tr.len')
img = regrid(img, upsample=True, interpolation='bilinear',y_sampling=.01, x_sampling=.01, aggregator='min').opts(cmap='gray', aspect=1, fig_inches=4, colorbar=True,xlim=(-1,3),ylim=(-2,2), fontsize=fs)
hv.save(img*(overlay), f'solomon_cusp.png', fmt='png')

data = pd.read_csv('pendulum.csv',names=['Re','Im','minimal tr.len']).drop_duplicates(subset=['Re','Im'],keep='last').set_index(['Re', 'Im'])
array = data.to_xarray()
ds = hv.Dataset(array,['Re', 'Im'], 'minimal tr.len')
img = ds.to(hv.Image, ['Re', 'Im'], 'minimal tr.len')
img = regrid(img, upsample=True, interpolation='bilinear',y_sampling=.01, x_sampling=.01, aggregator='min').opts(cmap='gray', aspect=1, fig_inches=4, colorbar=True,xlim=(0,2),ylim=(0,2), fontsize=fs)
hv.save(img, f'pendulum.png', fmt='png')

overlay = hv.Points([(1,0)]).opts(color="white",marker="o",s=50)*hv.Text(1,-.15,'cusp', fontsize=8).opts(color="white")*hv.Points([(0,0)]).opts(color="white",marker="o",s=50)*hv.Text(0,.15,'orbifold', fontsize=8).opts(color="white")
data = pd.read_csv('schottky_slice_hard_rand.csv',names=['Re','Im','minimal tr.len']).drop_duplicates(subset=['Re','Im'],keep='last').set_index(['Re', 'Im'])
array = data.to_xarray()
ds = hv.Dataset(array,['Re', 'Im'], 'minimal tr.len')
img = ds.to(hv.Image, ['Re', 'Im'], 'minimal tr.len')
img = regrid(img, upsample=True, interpolation='bilinear',y_sampling=.01, x_sampling=.01, aggregator='min').opts(cmap='gray', aspect=1, fig_inches=4, colorbar=True,xlim=(-2,2),ylim=(-2,2), fontsize=fs)
hv.save(img*(overlay), f'schottky_slice_hard_qand.png', fmt='png')
