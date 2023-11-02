from sage.all import *

import time
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
      sym = ManinSymbol(N, qs[j], qs[j-1]).to_generator()[1]
      out -= stb[sym]
    else:
      sym = ManinSymbol(N, -qs[j], qs[j-1]).to_generator()[1]
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

  def to_generator(self) -> 'ManinSymbol':
    return find_generator(self)

  def to_modsym(self) -> ModularSymbol:
    g, x, y = xgcd(self.c, self.d)
    a, b = y, -x
    return ModularSymbol(a, b, self.c, self.d)

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

  def map(self, f: Callable[[ManinSymbol],'ManinSymbolCombination'], M) -> 'ManinSymbolCombination':
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
def find_generator(manin_sym : ManinSymbol) -> ManinSymbol:
  N = manin_sym.N
  gens = manin_generators(N)
  for i, gen in enumerate(gens):
    if manin_sym == gen:
      return i, gen
  print("Could not find generator for", manin_sym)
  assert False

@cache
def ST_relation_matrix(N):
  gens = manin_generators(N)
  n_gens = len(gens)

  def fgi(c, d):
    return find_generator(ManinSymbol(N, c, d))[0]

  S_rels = []
  T_rels = []

  done_S = set()
  done_T = set()
  for x in range(n_gens):
    c, d = gens[x].tuple()
    if x not in done_S:
      Sx = fgi(d, -c)
      rel = (x, Sx)
      S_rels.append(rel)
      done_S.update(rel)
    if x not in done_T:
      Tx = fgi(c+d, -c)
      T2x = fgi(d, -(c+d))
      rel = (x, Tx, T2x)
      T_rels.append(rel)
      done_T.update(rel)

  # print([tuple(gens[i] for i in rel) for rel in S_rels])
  # print([tuple(gens[i] for i in rel) for rel in T_rels])

  def rel_to_row(rel):
    return [1 if i in rel else 0 for i in range(n_gens)]

  return gens, Matrix([rel_to_row(rel) for rel in S_rels + T_rels])

@cache
def manin_basis(N) -> list[ManinSymbol]:
  gens, mtx = ST_relation_matrix(N)
  nonpivots = mtx.nonpivots()
  # print(nonpivots)
  return [gens[i] for i in nonpivots]

@cache
def manin_basis_idx(N):
  return {manin_sym : i for i, manin_sym in enumerate(manin_basis(N))}

@cache
def sym_to_basis(N):
  gens, mtx = ST_relation_matrix(N)
  nonpivots = mtx.nonpivots()
  rref = mtx.rref()
  out = {}
  for n in nonpivots:
     out[gens[n]] = gens[n].as_combination()
  for row in rref:
    pivot = None
    val = ManinSymbolCombination.zero(N)
    for i, c in enumerate(row):
      if c != 0:
        if i in nonpivots:
          val += gens[i].as_combination(-c)
        else:
          if pivot is None:
            pivot = i
          else:
            assert False
    if pivot is not None:
      out[gens[pivot]] = val
  assert len(out) == len(gens)
  return out

@cache
def cusp_equivalent(C1, C2, N):
  # print ("cusp_equivalent called:", C1, C2, N)
  ((c1, d1), (c2, d2)) = C1, C2
  g1, s1, _ = xgcd(c1, d1)
  g2, s2, _ = xgcd(c2, d2)
  assert (g1 == 1) and (g2 == 1)
  g = gcd(d1 * d2, N)
  return (s1 * d2 - s2 * d1) % g == 0

def boundary_map(manin_sym : ManinSymbol):
  c, d = manin_sym.tuple()
  g, b, a = xgcd(c, d)
  assert (g == 1)
  return ((a, c), (b, d))

@cache
def boundary_map_matrix(N):
  cusps = []

  def representative(cusp):
    for C in cusps:
      if cusp_equivalent(cusp, C, N):
        return C
    # new cusp
    cusps.append(cusp)
    return cusp

  basis = manin_basis(N)

  out = []
  for manin_sym in basis:
    # print(manin_sym)
    c1, c2 = boundary_map(manin_sym)
    c1, c2 = representative(c1), representative(c2)
    if c1 == c2:
      out += [dict()]
    else:
      out += [{representative(c1) : +1, representative(c2) : -1}]

  def dict_to_row(d):
    return [d.get(cusp, 0) for cusp in cusps]

  return basis, Matrix([dict_to_row(d) for d in out])

@cache
def boundary_map_kernel(N):
  basis, mtx = boundary_map_matrix(N)
  kernel = mtx.kernel().basis_matrix()
  return basis, kernel

@cache
def oldspace_map(manin_sym : ManinSymbol, d, M) -> ManinSymbolCombination:
  N = manin_sym.N
  assert N % (d * M) == 0
  return manin_sym.to_modsym().left_action_by((d, 0, 0, 1)).to_manin_comb(M)

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


def newspace(N):
  mbasis, current_basis = boundary_map_kernel(N)
  cuspidal_basis = current_basis
  # print("Manin basis:\n" + str(mbasis))
  # print("Boundary map kernel:\n" + str(current_basis))
  for M in divisors(N):
    if M > 10 and M != N:
      for d in divisors(N/M):
        # print(f"Taking oldspace map with M = {M} and d = {d}")
        combs = rows_to_manin_combs(mbasis, current_basis)
        map_output = [oldspace_map_comb(comb, d, M) for comb in combs]
        mtx = manin_combs_to_matrix(map_output)
        # print("Output mtx:\n" + str(mtx))
        new_basis_in_old_basis = mtx.kernel().basis_matrix()
        # print("Output mtx kernel:\n" + str(new_basis_in_old_basis))
        current_basis = new_basis_in_old_basis * current_basis
        # print("New kernel:\n" + str(current_basis))
  return mbasis, cuspidal_basis, current_basis

if __name__ == "__main__":
  ManinSymbol.__repr__ = ManinSymbol.__str__

  # N = 30
  # for sym in manin_basis(N):
  #   print(sym)
  #   print(sym.to_modsym())
  #   print(sym.to_modsym().to_manin_comb(N))


  limit = 200

  start_time = time.time()
  new_bases = dict()
  for N in range(1, limit + 1):
    #  print(f">> Starting computation for {N}")
    mbasis, cuspidal_basis, current_basis = newspace(N)
    # print(mbasis)
    # print(current_basis)
    new_bases[N] = (mbasis, cuspidal_basis, current_basis)

  end_time = time.time()

  for N in range(1, limit + 1):
    print(N, new_bases[N][1].ncols(), new_bases[N][1].nrows(), new_bases[N][2].nrows())

  print("Computation finished in ", (end_time - start_time), "seconds")

