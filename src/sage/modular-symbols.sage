from sage.all import *

from dataclasses import dataclass
from functools import cache
from typing import Any, Callable

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

  def map(self, f: Callable[[ManinSymbol],'ManinSymbolCombination']) -> 'ManinSymbolCombination':
    out = ManinSymbolCombination.zero()
    for sym, coeff in self.components.items():
      out += coeff * f(sym)
    return out

  def __mul__(self, other):
    return ManinSymbolCombination(
      self.N,
      {k: other * v for k, v in self.components.items()}
    )

  def __add__(self, other):
    assert isinstance(other, ManinSymbolCombination)
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
def ST_relation_matrix(N):
  gens = manin_generators(N)
  n_gens = len(gens)

  def fgc(c, d):
    s = ManinSymbol(N, c, d)
    for i, gen in enumerate(gens):
      if s == gen:
        return i
    print(s)
    assert False

  S_rels = []
  T_rels = []

  done_S = set()
  done_T = set()
  for x in range(n_gens):
    c, d = gens[x].tuple()
    if x not in done_S:
      Sx = fgc(d, -c)
      rel = (x, Sx)
      S_rels.append(rel)
      done_S.update(rel)
    if x not in done_T:
      Tx = fgc(c+d, -c)
      T2x = fgc(d, -(c+d))
      rel = (x, Tx, T2x)
      T_rels.append(rel)
      done_T.update(rel)

  # print([tuple(gens[i] for i in rel) for rel in S_rels])
  # print([tuple(gens[i] for i in rel) for rel in T_rels])

  def rel_to_row(rel):
    return [1 if i in rel else 0 for i in range(n_gens)]

  return gens, Matrix([rel_to_row(rel) for rel in S_rels + T_rels])

@cache
def manin_basis(N):
  gens, mtx = ST_relation_matrix(N)
  basis = mtx.nonpivots()
  rref = mtx.rref()
  return gens, basis, rref

@cache
def sym_to_basis(N):
  gens, nonpivots, rref = manin_basis(N)
  print(gens)
  print(nonpivots)
  print(rref)
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

if __name__ == "__main__":
  ManinSymbol.__repr__ = ManinSymbol.__str__
  out = sym_to_basis(12)
  for k, v in out.items():
    print(k, v)
