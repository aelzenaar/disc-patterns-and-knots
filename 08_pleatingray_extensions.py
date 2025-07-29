from bella import slices, farey,cayley
from mpmath import mp
import holoviews as hv
import pandas as pd
hv.extension('bokeh')
from functools import reduce

depth = 30

def reflect_list(X):
  return X + [-x for x in X]

def complex_list_to_points(X):
  return hv.Points([(float(x.real), float(x.imag)) for x in X])

fig8knot = farey.approximate_pleating_ray(3,5,mp.inf,mp.inf,R=20,N=100,end_at_cusp = False)[-1]

roots1 = complex_list_to_points(reflect_list([fig8knot])).opts(marker="+", size=20, frame_width=500,frame_height=300, data_aspect=1, color='red')\
          .redim(x=hv.Dimension('x', range=(-5,5)),y=hv.Dimension('y', range=(-3, 3)))  * hv.Ellipse(0,0,2).opts(color='gray')

groups = []
for r,s in farey.walk_tree_bfs(11):
  groups.append(farey.approximate_pleating_ray(r,s,mp.inf,mp.inf, R=50, N=100, end_at_cusp=False)[-1])

roots2 = complex_list_to_points(reflect_list(groups)).opts(marker="+", size=15, frame_width=500,frame_height=300, data_aspect=1, color='black')\
          .redim(x=hv.Dimension('x', range=(-5,5)),y=hv.Dimension('y', range=(-3, 3))) * roots1


df = slices.elliptic_exterior(mp.inf, mp.inf, depth)
roots2 *= hv.Scatter(df,  kdims=['x'],vdims=['y']).opts(marker="dot", size=2, frame_width=500,frame_height=300, data_aspect=1, color='black')\
          .redim(x=hv.Dimension('x', range=(-5,5)),y=hv.Dimension('y', range=(-3, 3)))


ray = farey.approximate_pleating_ray(3,5,mp.inf,mp.inf,R=20,N=100,end_at_cusp = False)

t = 0
for ρ in ray:
  if mp.fabs(ρ) > 3:
    print('(Skipping large ρ)')
    continue

  t += 1

  X = mp.matrix([[1,1],[0,1]])
  Y = mp.matrix([[1,0],[ρ,1]])
  G = cayley.GroupCache([X,Y],names=['X','Y'])
  tr = cayley.simple_tr(G[G.fancyword_to_word('XyxYXYxyXY')])
  print(f'{t:04d}\t\t{mp.nstr(ρ)}\t\t{mp.nstr(2*mp.acosh(tr/2))}')

  roots2 *= hv.Points([(float(ρ.real),float(ρ.imag))]).opts(marker="x",size=20,color="blue")

  if float(tr.real) < -2:
    continue

  hv.save(roots2, f'08/slice_{t:04d}.png', fmt='png')

  limit_points = pd.concat([G.coloured_limit_set_fast(10**6,seed=0),\
                            G.coloured_limit_set_fast(10**6,seed=1),\
                            G.coloured_limit_set_fast(10**6,seed=-1),\
                            G.coloured_limit_set_fast(10**6,seed=2),\
                            G.coloured_limit_set_fast(10**6,seed=-2)])
  scatter = hv.Scatter(limit_points, kdims = ['x'], vdims = ['y','colour'])\
              .opts(marker = "dot", size=1,\
                    data_aspect=1, frame_width=600,\
                    color = 'colour', cmap='glasbey_cool')\
              .redim(x=hv.Dimension('x', range=(-1.5,1.5)),\
                      y=hv.Dimension('y', range=(-1.5,1.5)))

  splittings = farey.standard_peripheral_generators(3,5)

  splittings = [('X','yxYXYxyXY','red'), ('XyxYXYx','yXY','green')]
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


      # fixed_points = [G.fixed_points(G.fancyword_to_word(u))[0]] + [G.fixed_points(G.fancyword_to_word(v))[0]] + [G.fixed_points(G.fancyword_to_word(u+v))[0]]
      # disc = cayley.circle_space_to_circle_or_line(mp.chop(cayley.circle_through_points(*fixed_points), tol=10**-15))
      # if(disc[-1] == True):
      #     pt1 = 100*disc[0] + (1-100)*disc[1]
      #     pt2 = -100*disc[0] + (1+100)*disc[1]
      #     splitting_scatters += [hv.Path([(float(pt1.real), float(pt1.imag)), (float(pt2.real), float(pt2.imag))]).opts(color=colour, alpha=0.5)]
      # else:
      #     splitting_scatters += [hv.Ellipse(float(disc[0].real), float(disc[0].imag), float(disc[1].real)*2).opts(color=colour, alpha=0.5)]

  for s in splitting_scatters:
    scatter *= s

  hv.save(scatter, f'08/limset_{t:04d}.png', fmt='png')

