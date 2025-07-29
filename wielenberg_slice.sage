import csv
import re
import itertools
import numpy as np
from math import isnan
from random import shuffle
import warnings
warnings.simplefilter("error")

import gc
gc.enable()

import multiprocessing

output_csv = "wielenberg_slice.csv"
# α = 2+2j
# β = 1-1j
α = 2*(1+I*sqrt(3))
β = -1+I*sqrt(3)

depth = 7

pixel_width=160
pixel_height=160
x_bounds = (-2,2)
y_bounds = (-2,2)

t = var('t')
P = Matrix([[1,α],[0,1]])
Q = Matrix([[1,β],[0,1]])
M = Matrix([[t,t**2-1],[1,t]])

generators = [(P,'P'),(Q,'Q'),(M,'M')]
generators += [(m**-1, l.swapcase()) for m,l in generators ]

def try_simplify_full(u):
    try:
        return u.expand().simplify()
    except:
        return u

list_of_words = []

last_row = [(Matrix([[1,0],[0,1]]), '')]
for k in range(0,depth):
    print(f"Computing words at height {k+1}/{depth}")
    next_row = []
    for m, l in last_row:
        for m2, l2 in generators:
            if len(l) == 0 or l2 != l[-1].swapcase():
                # print(l+l2)
                next_row.append((try_simplify_full(m*m2), l+l2))
    list_of_words += next_row
    last_row = next_row


print(f'Filtering out conjugates...')
length_before_filter = len(list_of_words)
def is_conjugate(w):
    if len(w) <= 1:
        return False
    if w[0] == w[-1].swapcase():
        return True
list_of_words = [x for x in list_of_words if not is_conjugate(x[1])]
print(f'Filtered out conjugates: from {length_before_filter} to {len(list_of_words)}, i.e. {len(list_of_words)/length_before_filter:.0%}')

print(f'Filtering out powers...')
length_before_filter = len(list_of_words)
p=re.compile('^(.+)\\1+$')
list_of_words = [x for x in list_of_words if not p.search(x[1])]
print(f'Filtered out powers: from {length_before_filter} to {len(list_of_words)}, i.e. {len(list_of_words)/length_before_filter:.0%}')


# Distance |tr^2 - 4|
print(f'Computing traces...')
list_of_traces = [try_simplify_full(x.trace()) for x,_ in list_of_words]

length_before_filter = len(list_of_traces)
# Filter words
print(f'Filtering out parabolics...')
list_of_traces = [x for x in list_of_traces if N(x.subs(t=-5+2j))**2-4 != 0]
print(f'Filtered out parabolics: from {length_before_filter} to {len(list_of_traces)}, i.e. {len(list_of_traces)/length_before_filter:.0%}')



def single_x_row(pixel):
    px_x,px_y = pixel
    perc = float((px_x*pixel_height + px_y)/(pixel_height*pixel_width))
    print(f'Computing heights, x={px_x} y={px_y}, approx perc = {perc:.0%}')
    x = float((px_x/pixel_width)*abs(x_bounds[1]-x_bounds[0]) + x_bounds[0])
    y = float((px_y/pixel_height)*abs(y_bounds[1]-y_bounds[0]) + y_bounds[0])
    def relen(tr):
        try:
            inp = N(tr.subs(t=(x+y*1j)))
            return float(N(abs(2*log(abs((inp + sqrt(inp**2-4))/2)))))
        except:
            return infinity
    height = min(filter(lambda x: not isnan(x), map(relen, list_of_traces)))
    return (float(x),float(y),float(height))

# single_x_row([0,3])

with open(output_csv, 'w') as f:
    c = csv.writer(f)
    with multiprocessing.Pool(maxtasksperchild=int(1)) as pool:
        prod = list(itertools.product(range(0,pixel_width),range(0,pixel_height)))
        # shuffle(prod)
        for height in pool.imap(single_x_row, prod):
            c.writerow(height)
            f.flush()

