from bella import slices, farey,cayley
from mpmath import mp
import holoviews as hv
import pandas as pd
hv.extension('bokeh')
from functools import reduce
from holoviews.operation.datashader import datashade
from holoviews.operation.resample import ResampleOperation2D
ResampleOperation2D.width=1600
ResampleOperation2D.height=1600


ray = [-0.773301 + 1.46771j, -0.715145 + 1.32335j, -0.57683 + 1.01168j, mp.exp(mp.pi*2j/3)]

t = 0
for ρ in ray:
  t += 1
  print(t)
  if t != 3:
    continue

  X = mp.matrix([[1,1],[0,1]])
  Y = mp.matrix([[1,0],[ρ,1]])
  G = cayley.GroupCache([X,Y],names=['X','Y'])
  A = G[G.fancyword_to_word('yxYXYxyXY')]
  B = G[G.fancyword_to_word('XyxYXYx')]
  C = G[G.fancyword_to_word('yXY')]
  G = cayley.GroupCache([X,Y,A,B,C],names=['X','Y','A','B','C'])


  limit_points = G.coloured_limit_set_fast(10**7)
  stuff = [limit_points]
  for lat in [-2,-1,0,1,2]:
        df = limit_points.copy()
        df['z'] = df['x']+1j*df['y']
        df['z'] = df['z'].apply(lambda z: complex(z+lat))
        df['x'] = df['z'].apply(lambda z: z.real)
        df['y'] = df['z'].apply(lambda z: z.imag)
        del df['z']
        stuff.append(df)
  limit_points = pd.concat(stuff, axis=0)

  scatter = datashade(hv.Scatter(limit_points, kdims = ['x'], vdims = ['y','colour'])\
              .opts(marker = "dot", size=1,\
                    data_aspect=1, frame_width=600,\
                    color = 'gray')\
              .redim(x=hv.Dimension('x', range=(-1.5,1.5)),\
                      y=hv.Dimension('y', range=(-1.5,1.5))), width=800, cmap=['gray'], min_alpha=0)

  splittings = farey.standard_peripheral_generators(3,5)

  splittings = [('X','yxYXYxyXY','k'), ('XyxYXYx','yXY','k')]
  splitting_scatters = []
  commcircles = [G.isometric_circle(G.fancyword_to_word('XyxYXYxyXY')),G.isometric_circle(G.fancyword_to_word('yxYXyxyXYx'))]
  splitting_scatters += [hv.Ellipse(float(centre.real), float(centre.imag), float(radius)*2).opts(color='black') for (centre, radius) in commcircles if radius != mp.inf]
  for u,v,colour in splittings:
      # print(''.join(u), ''.join(v))
      H = G.subgroup([G.fancyword_to_word(u), G.fancyword_to_word(v)])
      limit_points_2 = H.coloured_limit_set_fast(10**5)
      splitting_scatters.append(hv.Scatter(limit_points_2, kdims = ['x'], vdims = ['y']).opts(marker = "dot", size = 5,  color = 'black'))
      isocircles = [H.isometric_circle(w) for w in H.free_cayley_graph_bfs(1)]
      splitting_scatters += [hv.Ellipse(float(centre.real), float(centre.imag), float(radius)*2).opts(color=colour) for (centre, radius) in isocircles if radius != mp.inf]


  for s in splitting_scatters:
    scatter *= s

  hv.save(scatter.opts(width=800, data_aspect=1), f'fig8path_{t}.png', fmt='png')

