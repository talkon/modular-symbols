'''
Very unoptimized implementation
'''

from sage.all import *

import time
import sys
from dataclasses import dataclass
from functools import cache
from typing import Any, Callable

def fraction_to_manin_comb(a: Integer, b: Integer, N: Integer) -> 'ManinSymbolCombination':
  if b == 0:
    return ManinSymbol(N, 1, 0).as_combination()

  out = ManinSymbolCombination.zero(N)
  stb = sym_to_basis(N)

  # not sure what's going on with signs, but this seems correct
  cf = continued_fraction(a / b)
  qs = [cf.q(i) for i in range(cf.length())]
  for j in range(1, cf.length()):
    if j % 2 == 1:
      # print(qs[j], qs[j-1])
      sym = ManinSymbol(N, qs[j], qs[j-1]).to_generator(omit_index = True)
      # print("=", stb[sym])
      out -= stb[sym]
    else:
      # print(-qs[j], qs[j-1])
      sym = ManinSymbol(N, -qs[j], qs[j-1]).to_generator(omit_index = True)
      # print("=", stb[sym])
      out -= stb[sym]

  return out

@dataclass
class ModularSymbol:
  a: Integer
  b: Integer
  c: Integer
  d: Integer

  def left_action_by(self, m):
    x, y, z, w = m
    a, b, c, d = self.a, self.b, self.c, self.d
    return ModularSymbol(
      x * a + y * c, x * b + y * d,
      z * a + w * c, z * b + w * d
    )

  def to_manin_comb(self, N):
    return fraction_to_manin_comb(self.b, self.d, N) - fraction_to_manin_comb(self.a, self.c, N)

@dataclass
class ManinSymbol:
  N: Integer
  c: Integer
  d: Integer

  def __init__(self, N, c, d):
    self.N = N
    self.c = c % N
    self.d = d % N

  def __str__(self):
    return f"({self.c}, {self.d})_{self.N}"

  def __eq__(self, other):
    assert isinstance(other, ManinSymbol)
    return self.N == other.N and ((self.c * other.d - self.d * other.c) % self.N) == 0

  def as_combination(self, coeff = 1) -> 'ManinSymbolCombination':
     return ManinSymbolCombination(self.N, {self : coeff})

  def to_generator(self, omit_index = True):
    return find_generator(self, omit_index)

  def to_modsym(self) -> ModularSymbol:
    g, x, y = xgcd(self.c, self.d)
    a, b = y, -x
    return ModularSymbol(a, b, self.c, self.d)

  def right_action_by(self, m) -> 'ManinSymbol':
    x, y, z, w = m
    return ManinSymbol(self.N, self.c * x + self.d * z, self.c * y + self.d * w).to_generator(omit_index = True)

  def modsym_left_action_by(self, m, M = None) -> 'ManinSymbolCombination':
    M = (self.N if M is None else M)
    return self.to_modsym().left_action_by(m).to_manin_comb(M)

  def __hash__(self):
     return hash((self.N, self.c, self.d))

  def tuple(self):
    return self.c, self.d

@dataclass
class ManinSymbolCombination:
  N: Integer
  components: dict[ManinSymbol, Any]

  @classmethod
  def zero(cls, N):
    return ManinSymbolCombination(N, dict())

  def map(self, f: Callable[[ManinSymbol],'ManinSymbolCombination'], M = None) -> 'ManinSymbolCombination':
    M = (self.N if M is None else M)
    out = ManinSymbolCombination.zero(M)
    for sym, coeff in self.components.items():
      out += coeff * f(sym)
    return out

  def __mul__(self, other):
    return ManinSymbolCombination(
      self.N,
      {k: other * v for k, v in self.components.items()}
    )

  def __rmul__(self, other):
    return self.__mul__(other)

  def __add__(self, other):
    assert isinstance(other, ManinSymbolCombination)
    # print(self, other)
    assert self.N == other.N
    components = {
      s : self.components.get(s, 0) + other.components.get(s, 0)
      for s in set(self.components) | set(other.components)
    }
    return ManinSymbolCombination(self.N, components)

  def __neg__(self):
    return ManinSymbolCombination(
      self.N,
      {k: -v for k, v in self.components.items()}
    )

  def __sub__(self, other):
    assert isinstance(other, ManinSymbolCombination)
    assert self.N == other.N
    return self + (-other)

  def __str__(self):
    out = ""
    for k, v in self.components.items():
      if v < -1:
        out += f"- {-v} * {k} "
      elif v == -1:
        out += f"- {k} "
      elif v == 0:
        pass
      elif v == 1:
        out += f"+ {k} "
      else:
        out += f"+ {v} * {k} "
    return out


'''
Manin symbol spaces
'''
@cache
def manin_generators(N : Integer) -> list[ManinSymbol]:
  l = divisors(N)
  out = [ManinSymbol(N, 0, 1)]
  for d in l[:-1]:
    cs = []
    M = N / d
    for x in range(N):
      if gcd(d, x) > 1:
        continue
      if any((x - c) % M == 0 for c in cs):
        continue
      cs.append(x)
    out += [ManinSymbol(N, d, c) for c in cs]
  return out

@cache
def find_generator(manin_sym : ManinSymbol, omit_index : bool = False) -> ManinSymbol:
  N = manin_sym.N
  gens = manin_generators(N)
  for i, gen in enumerate(gens):
    if manin_sym == gen:
      if omit_index:
        return gen
      else:
        return i, gen
  print("Could not find generator for", manin_sym)
  assert False

@cache
def ST_relation_matrix(N):
  gens = manin_generators(N)

  def find_generator_index(c, d):
    return find_generator(ManinSymbol(N, c, d))[0]

  # Compute P relations (eta relations / relations in Cremona Ch 2.5)
  done_P = set()
  filt_gens = []
  gen_index_to_filt_gen_index = {}
  P_equiv = {}
  i = 0
  for x, gen in enumerate(gens):
    c, d = gen.tuple()
    if x not in done_P:
      Px = find_generator_index(-c, d)
      gen_index_to_filt_gen_index[x] = i
      gen_index_to_filt_gen_index[Px] = i # if x = Px, this is fine
      i += 1
      filt_gens.append(gen)
      if x != Px:
        P_equiv[x] = Px
      done_P.update((x, Px))

  gen_to_filt_gen = {
    gens[k] : filt_gens[v] for k, v in gen_index_to_filt_gen_index.items()
  }

  print("# eta_gens:", len(filt_gens))

  # Compute S and T relations
  done_S = set()
  done_T = set()
  S_rels = []
  T_rels = []
  for x, gen in enumerate(gens):
    c, d = gen.tuple()
    # TODO: potential optimization: avoid creating duplicate rows
    if x not in done_S:
      Sx = find_generator_index(d, -c)
      rel = {x: 1}
      rel[Sx] = rel.get(Sx, 0) + 1
      S_rels.append(rel)
      done_S.update(rel)
    if x not in done_T:
      Tx = find_generator_index(c+d, -c)
      T2x = find_generator_index(d, -(c+d))
      rel = {x: 1}
      rel[Tx] = rel.get(Tx, 0) + 1
      rel[T2x] = rel.get(T2x, 0) + 1
      T_rels.append(rel)
      done_T.update(rel)

  def rel_to_filt(rel):
    out = {}
    for k, v in rel.items():
      fk = gen_index_to_filt_gen_index[k]
      out[fk] = out.get(fk, 0) + v
    return out

  S_rels = [rel_to_filt(rel) for rel in S_rels]
  T_rels = [rel_to_filt(rel) for rel in T_rels]

  def rel_to_row(rel):
    return tuple(rel.get(i, 0) for i in range(len(filt_gens)))

  rows = list(set(rel_to_row(rel) for rel in S_rels + T_rels))
  relation_matrix = Matrix(rows)
  # print(relation_matrix)

  return gen_to_filt_gen, filt_gens, relation_matrix

@cache
def manin_basis(N) -> list[ManinSymbol]:
  _, gens, mtx = ST_relation_matrix(N)
  nonpivots = mtx.nonpivots()
  # print(nonpivots)
  return [gens[i] for i in nonpivots]

@cache
def manin_basis_idx(N):
  return {manin_sym : i for i, manin_sym in enumerate(manin_basis(N))}

@cache
def sym_to_basis(N):
  gen_to_filt_gen, filt_gens, mtx = ST_relation_matrix(N)
  nonpivots = mtx.nonpivots()
  rref = mtx.rref()

  filt_gen_to_basis = {}
  for n in nonpivots:
     filt_gen_to_basis[filt_gens[n]] = filt_gens[n].as_combination()
  for row in rref:
    pivot = None
    val = ManinSymbolCombination.zero(N)
    for i, c in enumerate(row):
      if c != 0:
        if i in nonpivots:
          val += filt_gens[i].as_combination(-c)
        else:
          if pivot is None:
            pivot = i
          else:
            assert False
    if pivot is not None:
      filt_gen_to_basis[filt_gens[pivot]] = val
  assert len(filt_gen_to_basis) == len(filt_gens)

  gen_to_basis = {
    k : filt_gen_to_basis[v] for k, v in gen_to_filt_gen.items()
  }
  assert len(gen_to_basis) == len(gen_to_filt_gen)

  return gen_to_basis

@cache
def cusp_equivalent(C1, C2, N):
  # print ("cusp_equivalent called:", C1, C2, N)
  ((c1, d1), (c2, d2)) = C1, C2
  g1, s1, _ = xgcd(c1, d1)
  g2, s2, _ = xgcd(c2, d2)
  assert (g1 == 1) and (g2 == 1)
  g = gcd(d1 * d2, N)
  return (s1 * d2 - s2 * d1) % g == 0

@cache
def boundary_map(manin_sym : ManinSymbol):
  c, d = manin_sym.tuple()
  g, b, a = xgcd(c, d)
  assert (g == 1)
  return ((a, c), (-b, d))

@cache
def boundary_map_matrix(N):
  cusps = []

  def representative(cusp):
    c, d = cusp
    negated_cusp = (-c, d)
    for C in cusps:
      if cusp_equivalent(cusp, C, N):
        return C
      elif cusp_equivalent(negated_cusp, C, N):
        return C
    # new cusp
    cusps.append(cusp)
    return cusp

  basis = manin_basis(N)

  out = []
  for manin_sym in basis:
    c1, c2 = boundary_map(manin_sym)
    c1, c2 = representative(c1), representative(c2)
    if c1 == c2:
      out += [dict()]
    else:
      # Note: signs seem different from what Sage uses, but this should not matter since we only care about the kernel?
      '''
      Sage code:
      M = ModularSymbols(112)
      print(M.boundary_map().matrix().str())
      '''
      out += [{c1 : +1, c2 : -1}]

  def dict_to_row(d):
    return [d.get(cusp, 0) for cusp in cusps]

  mtx = Matrix([dict_to_row(d) for d in out])
  return basis, mtx

@cache
def boundary_map_kernel(N):
  basis, mtx = boundary_map_matrix(N)
  # print("Manin basis:\n" + str(basis))
  # print("Boundary map matrix:\n" + str(mtx))
  kernel = mtx.kernel().basis_matrix()
  return basis, kernel

@cache
def oldspace_map(manin_sym : ManinSymbol, d, M) -> ManinSymbolCombination:
  N = manin_sym.N
  assert N % (d * M) == 0
  out = manin_sym.modsym_left_action_by((d, 0, 0, 1), M)
  return out

def oldspace_map_comb(manin_comb : ManinSymbolCombination, d, M) -> ManinSymbolCombination:
  return manin_comb.map(lambda manin_sym : oldspace_map(manin_sym, d, M), M)

def rows_to_manin_combs(sym_basis : list[ManinSymbol], matrix):
  if len(sym_basis) == 0:
    return []

  N = sym_basis[0].N
  out = []
  for row in matrix:
    assert len(row) == len(sym_basis)
    out.append(ManinSymbolCombination(N, {k : v for k, v in zip(sym_basis, row)}))

  return out

def manin_combs_to_matrix(manin_combs : list[ManinSymbolCombination]):
  if len(manin_combs) == 0:
    return Matrix([[]])

  N = manin_combs[0].N
  return Matrix([[comb.components.get(manin_sym, 0) for manin_sym in manin_basis(N)] for comb in manin_combs])

@cache
def newspace(N):
  mbasis, current_basis = boundary_map_kernel(N)
  cuspidal_basis = current_basis
  print("Boundary map kernel:\n" + str(current_basis))
  for M in divisors(N):
    if M > 10 and M != N:
      for d in divisors(N/M):
        print(f"Taking oldspace map with M = {M} and d = {d}")
        combs = rows_to_manin_combs(mbasis, current_basis)
        map_output = [oldspace_map_comb(comb, d, M) for comb in combs]
        mtx = manin_combs_to_matrix(map_output)
        print("Output mtx:\n" + str(mtx))
        new_basis_in_old_basis = mtx.kernel().basis_matrix()
        print("Output mtx kernel:\n" + str(new_basis_in_old_basis))
        current_basis = new_basis_in_old_basis * current_basis
        print("New kernel:\n" + str(current_basis))
        print("New kernel dimension:", str(current_basis.nrows()))
  return mbasis, cuspidal_basis, current_basis

# from Cremona ch.2
@cache
def heilbronn_matrices(p):
  out = []
  out.append((1, 0, 0, p))
  for r in range(p):
    x1, x2 = p, -r
    y1, y2 = 0, 1
    a, b = -p, r
    out.append((x1, x2, y1, y2))
    while b != 0:
      q = Integer(round(a/b))
      c = a - b * q
      a = -b
      b = c
      x3 = q * x2 - x1
      x1 = x2
      x2 = x3
      y3 = q * y2 - y1
      y1 = y2
      y2 = y3
      out.append((x1, x2, y1, y2))
  return out

@cache
def hecke_action(manin_sym : ManinSymbol, p) -> ManinSymbolCombination:
  out = ManinSymbolCombination.zero(manin_sym.N)
  stb = sym_to_basis(manin_sym.N)
  for mtx in heilbronn_matrices(p):
    x = stb[manin_sym.right_action_by(mtx)]
    # print(manin_sym.right_action_by(mtx), " -> ", x)
    out += x
  return out

def hecke_action_comb(manin_comb : ManinSymbolCombination, p) -> ManinSymbolCombination:
  return manin_comb.map(lambda manin_sym : hecke_action(manin_sym, p))

# subspace basis should be in rref form
def hecke_action_matrix(mbasis, subspace_basis, p):
  combs = rows_to_manin_combs(mbasis, subspace_basis)
  map_output = [hecke_action_comb(comb, p) for comb in combs]
  mtx = manin_combs_to_matrix(map_output)
  pivots = subspace_basis.pivots()
  return mtx.matrix_from_columns(pivots)

# TODO: there's a lot of common code between actions of various operators, perhaps we could clean this up?
@cache
def atkin_lehner_matrix(N, Q):
  _, a, b = xgcd(Q, N/Q)
  return (Q, -b, N, a * Q)

@cache
def atkin_lehner_action(manin_sym : ManinSymbol, Q) -> ManinSymbolCombination:
  N = manin_sym.N
  return manin_sym.modsym_left_action_by(atkin_lehner_matrix(N, Q), N)

def atkin_lehner_action_comb(manin_comb : ManinSymbolCombination, Q) -> ManinSymbolCombination:
  return manin_comb.map(lambda manin_sym : atkin_lehner_action(manin_sym, Q))

def atkin_lehner_action_matrix(mbasis, subspace_basis, Q):
  combs = rows_to_manin_combs(mbasis, subspace_basis)
  map_output = [atkin_lehner_action_comb(comb, Q) for comb in combs]
  mtx = manin_combs_to_matrix(map_output)
  pivots = subspace_basis.pivots()
  return mtx.matrix_from_columns(pivots)

# decomposes a subspace with the given basis into simple (mtx)-modules
def decompose(mtx, basis):
  print("decompose() called")
  print("mtx:\n", mtx)
  print("basis:\n", basis)
  pivots = basis.pivots()
  mtx_on_basis = (basis * mtx).matrix_from_columns(pivots)
  print("mtx_on_basis:\n", basis * mtx, "\n", mtx_on_basis)
  out = []
  factors = factor(mtx_on_basis.minimal_polynomial())
  print("factors: ", factors)
  for poly, _ in factors:
    new_basis = poly(mtx_on_basis).kernel().basis_matrix()
    out += [(new_basis * basis, poly.degree() == new_basis.nrows())]
  print("decompose() output:\n", out)
  return out

@cache
def newform_subspaces(N, use_atkin_lehner = True):
  mbasis, _, newspace_basis = newspace(N)
  # print(mbasis)
  # print(newspace_basis)
  if newspace_basis.nrows() == 0:
    return []
  decomposition = []
  remaining = [identity_matrix(newspace_basis.nrows())]

  # Atkin-Lehner involutions
  # print("Atkin-Lehner involutions")
  if use_atkin_lehner:
    for p, l in factor(N):
      Q = p ** l
      mtx = atkin_lehner_action_matrix(mbasis, newspace_basis, Q)
      # print("Hecke matrix:", mtx)
      new_remaining = []
      for basis in remaining:
        # TODO: for Atkin-Lehner operators we know that the eigenvalues are +1 or -1,
        # so we shouldn't need to use a general decompose().
        for new_basis, _ in decompose(mtx, basis):
          if new_basis.nrows() > 0:
            new_remaining.append(new_basis)
      remaining = new_remaining

  # Hecke operators
  # print("Hecke operators")
  p = 0
  while True:
    p = next_prime(p)
    if len(remaining) == 0:
      break
    if N % p != 0:
      # print("p =", p)
      mtx = hecke_action_matrix(mbasis, newspace_basis, p)
      # print("Hecke matrix:", mtx)
      new_remaining = []
      for basis in remaining:
        print(basis.nrows(), basis.ncols())
        # print("Basis:\n" + str(basis))
        for new_basis, done in decompose(mtx, basis):
          if done:
            decomposition.append(new_basis)
          else:
            new_remaining.append(new_basis)
      remaining = new_remaining
  return decomposition

def perft(limit, use_atkin_lehner):
  out = dict()

  def print_out_item(cbasis, nbasis, t, decomp):
    decomp_dims = sorted(mtx.nrows() for mtx in decomp)
    print(f"{N:3} {cbasis.nrows():3} {nbasis.nrows():3} {round(t, 4):8.4f} {decomp_dims}", flush=True)


  abs_start_time = time.time()

  for N in range(1, limit+1):
    start_time = time.time()
    mbasis, cbasis, nbasis = newspace(N)
    decomposition = newform_subspaces(N, use_atkin_lehner)
    out[N] = (cbasis, nbasis, time.time() - start_time, decomposition)
    print_out_item(*out[N])

  print("Total time elapsed:", time.time() - abs_start_time)
