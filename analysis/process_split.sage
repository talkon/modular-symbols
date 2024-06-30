from sage.all import *

fname = "split"

def header_line() -> str:
  return "N,category,tau,omega,num_spaces,execution time (s),factors,split\n"

def process_line(line: str) -> str:
  n, split, time = line.split(":")
  n = Integer(n)
  tau = len(n.divisors())
  factors = list(factor(n))
  omega = len(factors)

  cat = "other"
  if n == 1:
    cat = "one"
  elif tau == 2:
    cat = "prime"
  elif all(x[1] == 1 for x in factors):
    cat = "squarefree composite"
  elif all(x[1] >= 2 for x in factors):
    cat = "powerful"

  num_spaces = len(eval(split))
  time = float(time)
  return f'{n},{cat},{tau},{omega},{num_spaces},{time},"{factors}","{split}"\n'

with open(f"{fname}.txt", "r") as i:
  out_lines = [header_line()] + [process_line(line) for line in i.readlines()]
  with open(f"{fname}.csv", "w") as o:
    o.writelines(out_lines)
