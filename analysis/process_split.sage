from sage.all import *

def header_line() -> str:
  return "n,tau,omega,num_spaces,time,split\n"

def process_line(line: str) -> str:
  n, split, time = line.split(":")
  n = Integer(n)
  tau = len(n.divisors())
  omega = len(n.prime_divisors())
  num_spaces = len(eval(split))
  time = float(time)
  return f'{n},{tau},{omega},{num_spaces},{time},"{split}"\n'

with open("split.txt", "r") as i:
  out_lines = [header_line()] + [process_line(line) for line in i.readlines()]
  with open("split.csv", "w") as o:
    o.writelines(out_lines)
