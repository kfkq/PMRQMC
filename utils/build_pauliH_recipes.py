"""
This code was written in support of the experiments carried out
* Nic Ezzell, Lev Barash, Itay Hen, Exact and universal quantum Monte Carlo estimators for energy susceptibility and fidelity susceptibility, arXiv:2408.03924 (2024).
* Nic Ezzell and Itay Hen, Advanced measurement techniques in quantum Monte Carlo: The permutation matrix representation approach, arXiv:2504.07295 (2025).

Description: Used to build the 2 qubit PRL model and the TFIM, and XXZ models on various lattices.
"""

# %%
from pauli_manipulations import PauliTerm, PauliH

####################################
#-----------------------------------
# Nearest-neighbor grids
#-----------------------------------
#################################### 
# %%
#################################### 
# Nearest-neighbor rectangular grid
####################################
def form_nn_rect_grid(r, c, pbc=0, linidx=1):
    """
    Returns 2D [r] x [c] rectangular grid of nearest-neighbor 
    connections.
    - pbc = 1 means use periodic boundary conditions (torus)
    - linidx = 1 means use indexing [1,2, ..., rxc]
    """
    grid = {}
    for i in range(r):
        for j in range(c):
            rn = (i, (j + 1) % c) # right neighbor
            dn = ((i + 1) % r, j) # down neighbor
            if linidx == 0:
                # form (i, j) --> [(i + 1, j), (i, j + 1)] grid
                node = (i, j)
                right = rn
                down = dn
            else:
                # convert (i, j) to linear index (arithmetization)
                node = (c * i + j) + 1
                right = (c * rn[0] + rn[1]) + 1
                down = (c * dn[0] + dn[1]) + 1
            # handle perodic boundary conditions true or false
            if pbc == 1:
                neighbors = [right, down]
            else:
                if i == r - 1 and j == c - 1:
                    neighbors = []
                elif i == r - 1:
                    neighbors = [right]
                elif j == c - 1:
                    neighbors = [down]
                else:
                    neighbors = [right, down]
            grid[node] = neighbors

    return grid

#################################### 
# Nearest-neighbor triangular grid
####################################
def form_nn_tri_grid(r, c, pbc=0, linidx=1):
    """
    Returns 2D [r] x [c] triangular grid of nearest-neighbor 
    connections. PBC not supported for now.
    - pbc = 1 means use periodic boundary conditions (torus-style)
    - linidx = 1 means use indexing [1,2, ..., rxc]
    """
    # first form square grid without PBC, with (i,j) index
    grid = form_nn_rect_grid(r, c, pbc, 0)
    # connect diagonals to make triangles (think tilted square tile)
    for point in grid:
        i, j = point
        dn = ((i + 1) % r, (j + 1) % c) # diagonal neighbor
        if pbc == 1:
            grid[point].append(dn)
        else:
            if i + 1 < r and j + 1 < c:
                grid[point].append((i + 1, j + 1))
    # convert to linear indexing if wanted
    if linidx == 1:
        lin_grid = {}
        for point in grid:
            i, j = point
            lin_point = (c * i + j) + 1
            lin_grid[lin_point] = []
            for npt in grid[point]:
                l_npt = c * npt[0] + npt[1] + 1
                lin_grid[lin_point].append(l_npt)
        grid = lin_grid
    
    return grid
# %%

####################################
#-----------------------------------
# Isotropic TFIM Hamiltonians
#-----------------------------------
#################################### 
# %%
#################################### 
# 1D TFIM 
####################################
def build_1d_tfim_xx(n, pbc=0):
    """
    Builds XX term in 1D-TFIM.
    [pbc] = 1: uses periodic boundary conditions
    [pbc] = 0: uses open boundary conditions
    """
    h = PauliH(n)
    for j in range(n-1):
        h.add_term(PauliTerm(-1.0, [j+1, j+2], ['X', 'X'], n))
    if pbc:
        h.add_term(PauliTerm(-1.0, [n, 1], ['X', 'X'], n))

    return h

def build_1d_tfim_z(n):
    """
    Builds Z driving term in 1D-TFIM.
    """
    h = PauliH(n)
    for j in range(n):
        h.add_term(PauliTerm(-1.0, [j+1], ['Z'], n))

    return h

def build_1d_tfim_zz(n, pbc=0):
    """
    Builds ZZ term in 1D-TFIM.
    [pbc] = 1: uses periodic boundary conditions
    [pbc] = 0: uses open boundary conditions
    """
    h = PauliH(n)
    for j in range(n-1):
        h.add_term(PauliTerm(-1.0, [j+1, j+2], ['Z', 'Z'], n))
    if pbc:
        h.add_term(PauliTerm(-1.0, [n, 1], ['Z', 'Z'], n))

    return h

def build_1d_tfim_x(n):
    """
    Builds X driving term in 1D-TFIM.
    """
    h = PauliH(n)
    for j in range(n):
        h.add_term(PauliTerm(-1.0, [j+1], ['X'], n))

    return h

#################################### 
# 2D square lattice TFIM
####################################
def build_square_tfim_xx(l, pbc=0):
    """
    Builds lxl XX term for 2D TFIM on square
    lattice, i.e. \sum_<i,j> X_iX_j.
    - pbc: use periodic boundary conditions?
    """
    n = l**2
    h = PauliH(n)
    grid = form_nn_rect_grid(l, l, pbc=pbc)
    for i in grid:
        for j in grid[i]:
            h.add_term(PauliTerm(-1.0, [i, j], ['X', 'X'], n))
    h.simplify()
    
    return h

def build_square_tfim_z(l):
    """
    Builds lxl Z term for 2D TFIM on square
    lattice, i.e. \sum_i Z_i. 
    """
    n = l**2
    h = PauliH(n)
    for i in range(1, n+1):
        h.add_term(PauliTerm(-1.0, [i], ['Z'], n))
    
    return h

def build_square_tfim_zz(l, pbc=0):
    """
    Builds lxl ZZ term for 2D TFIM on square
    lattice, i.e. \sum_<i,j> X_iX_j.
    - pbc: use periodic boundary conditions?
    """
    n = l**2
    h = PauliH(n)
    grid = form_nn_rect_grid(l, l, pbc=pbc)
    for i in grid:
        for j in grid[i]:
            h.add_term(PauliTerm(-1.0, [i, j], ['Z', 'Z'], n))
    h.simplify()
    
    return h

def build_square_tfim_x(l):
    """
    Builds lxl Z term for 2D TFIM on square
    lattice, i.e. \sum_i Z_i. 
    """
    n = l**2
    h = PauliH(n)
    for i in range(1, n+1):
        h.add_term(PauliTerm(-1.0, [i], ['X'], n))
    
    return h

#################################### 
# 2D triangular lattice TFIM
####################################
def build_triangle_tfim_xx(l, pbc=0):
    """
    Builds lxl XX term for 2D TFIM on
    triangular lattice, i.e., \sum_<i,j> X_iX_j.
    """
    n = l**2
    h = PauliH(n)
    grid = form_nn_tri_grid(l, l, pbc)
    for i in grid:
        for j in grid[i]:
            h.add_term(PauliTerm(-1.0, [i, j], ['X', 'X'], n))
    h.simplify()
    
    return h

def build_triangle_tfim_z(l):
    """
    Builds lxl Z term for 2D TFIM on
    triangular lattice, i.e, \sum_i Z_i.
    """
    n = l**2
    h = PauliH(n)
    for i in range(1, n+1):
        h.add_term(PauliTerm(-1.0, [i], ['Z'], n))
    
    return h

def build_triangle_tfim_zz(l, pbc=0):
    """
    Builds lxl ZZ term for 2D TFIM on
    triangular lattice, i.e., \sum_<i,j> X_iX_j.
    """
    n = l**2
    h = PauliH(n)
    grid = form_nn_tri_grid(l, l, pbc)
    for i in grid:
        for j in grid[i]:
            h.add_term(PauliTerm(-1.0, [i, j], ['Z', 'Z'], n))
    h.simplify()
    
    return h

def build_triangle_tfim_x(l):
    """
    Builds lxl X term for 2D TFIM on
    triangular lattice, i.e, \sum_i Z_i.
    """
    n = l**2
    h = PauliH(n)
    for i in range(1, n+1):
        h.add_term(PauliTerm(-1.0, [i], ['X'], n))
    
    return h
# %%

####################################
#-----------------------------------
# Isotropic XXZ Hamiltonians
#-----------------------------------
#################################### 
# %%
#################################### 
# 1D XXZ 
####################################
def build_1d_xxz_xxyy(n, pbc=0, div2=True):
    """
    Builds XX + YY terms for 1D XXZ model.
    [pbc] = 1: uses periodic boundary conditions
    [pbc] = 0: uses open boundary conditions
    [div2] = True: uses Sx = (1/2)X operators
    """
    coeff = 0.25 if div2 == True else 1.0
    h = PauliH(n)
    for j in range(n-1):
        h.add_term(PauliTerm(coeff, [j+1, j+2], ['X', 'X'], n))
        h.add_term(PauliTerm(coeff, [j+1, j+2], ['Y', 'Y'], n))
    if pbc:
        h.add_term(PauliTerm(coeff, [n, 1], ['X', 'X'], n))
        h.add_term(PauliTerm(coeff, [n, 1], ['Y', 'Y'], n))

    return h

def build_1d_xxz_zz(n, pbc=0, div2=True):
    """
    Builds ZZ terms for 1D XXZ model.
    [pbc] = 1: uses periodic boundary conditions
    [pbc] = 0: uses open boundary conditions
    [div2] = True: uses Sz = (1/2)Z operators
    """
    coeff = 0.25 if div2 == True else 1.0
    h = PauliH(n)
    for j in range(n-1):
        h.add_term(PauliTerm(coeff, [j+1, j+2], ['Z', 'Z'], n))
    if pbc:
        h.add_term(PauliTerm(coeff, [n, 1], ['Z', 'Z'], n))

    return h

def build_1d_xxz_z(n, div2=True):
    """
    Builds ZZ terms for 1D XXZ model.
    [pbc] = 1: uses periodic boundary conditions
    [pbc] = 0: uses open boundary conditions
    [div2] = True: uses Sz = (1/2)Z operators
    """
    coeff = 0.5 if div2 == True else 1.0
    h = PauliH(n)
    for j in range(n-1):
        h.add_term(PauliTerm(-coeff, [j+1], ['Z'], n))

    return h

#################################### 
# 2D square lattice XXZ 
####################################
def build_square_xxz_xxyy(l, pbc=0, div2=True):
    """
    Builds lxl XX+YY terms for 2D XXZ on square
    lattice, i.e. \sum_<i,j> (X_iX_j + Y_iY_j).
    - pbc: use periodic boundary conditions?
    [div2] = True: uses Sx = (1/2)X operators
    """
    coeff = 0.25 if div2 == True else 1.0
    n = l**2
    h = PauliH(n)
    grid = form_nn_rect_grid(l, l, pbc=pbc)
    for i in grid:
        for j in grid[i]:
            h.add_term(PauliTerm(coeff, [i, j], ['X', 'X'], n))
            h.add_term(PauliTerm(coeff, [i, j], ['Y', 'Y'], n))
    h.simplify()
    
    return h

def build_square_xxz_zz(l, pbc=0, div2=True):
    """
    Builds lxl ZZ terms for 2D XXZ on square
    lattice, i.e. \sum_<i,j> (Z_iZ_j).
    - pbc: use periodic boundary conditions?
    [div2] = True: uses Sz = (1/2)Z operators
    """
    coeff = 0.25 if div2 == True else 1.0
    n = l**2
    h = PauliH(n)
    grid = form_nn_rect_grid(l, l, pbc=pbc)
    for i in grid:
        for j in grid[i]:
            h.add_term(PauliTerm(coeff, [i, j], ['Z', 'Z'], n))
    h.simplify()
    
    return h

def build_square_xxz_z(l, div2=True):
    """
    Builds lxl Z term for 2D XXZ on square
    lattice, i.e. \sum_i Z_i. 
    [div2] = True: uses Sz = (1/2)Z operators
    """
    coeff = 0.5 if div2 == True else 1.0
    n = l**2
    h = PauliH(n)
    for i in range(1, n+1):
        h.add_term(PauliTerm(-coeff, [i], ['Z'], n))
    
    return h


#################################### 
# 2D triangular lattice XXZ 
####################################
def build_triangle_xxz_xxyy(l, pbc=0, div2=True):
    """
    Builds lxl XX+YY terms for 2D XXZ on triangular
    lattice, i.e. \sum_<i,j> (X_iX_j + Y_iY_j).
    - pbc: use periodic boundary conditions?
    [div2] = True: uses Sx = (1/2)X operators
    """
    coeff = 0.25 if div2 == True else 1.0
    n = l**2
    h = PauliH(n)
    grid = form_nn_tri_grid(l, l, pbc)
    for i in grid:
        for j in grid[i]:
            h.add_term(PauliTerm(coeff, [i, j], ['X', 'X'], n))
            h.add_term(PauliTerm(coeff, [i, j], ['Y', 'Y'], n))
    h.simplify()
    
    return h

def build_triangle_xxz_zz(l, pbc=0, div2=True):
    """
    Builds lxl ZZ terms for 2D XXZ on triangular
    lattice, i.e. \sum_<i,j> (Z_iZ_j).
    - pbc: use periodic boundary conditions?
    [div2] = True: uses Sz = (1/2)Z operators
    """
    coeff = 0.25 if div2 == True else 1.0
    n = l**2
    h = PauliH(n)
    grid = form_nn_tri_grid(l, l, pbc)
    for i in grid:
        for j in grid[i]:
            h.add_term(PauliTerm(coeff, [i, j], ['Z', 'Z'], n))
    h.simplify()
    
    return h

def build_triangle_xxz_z(l, div2=True):
    """
    Builds lxl Z term for 2D XXZ on
    triangular lattice, i.e. \sum_i Z_i. 
    [div2] = True: uses Sz = (1/2)Z operators
    """
    coeff = 0.5 if div2 == True else 1.0
    n = l**2
    h = PauliH(n)
    for i in range(1, n+1):
        h.add_term(PauliTerm(-coeff, [i], ['Z'], n))
    
    return h
# %%
# %%
