#!/usr/bin/env python

import subprocess
import sys
from dataclasses import dataclass
from statistics import geometric_mean
from colorama import Fore, Style

SPLIT_FILE = "./analysis/split.txt"
SPLIT_FILE_CALIBRATE = 13.4

BINARY = "./bin/main"
TESTS = [
  {
    "prime": (1009, 1123, 1523, 2053, 2531, 3343, 3719, 4133, 4889),
    "composite": (1846, 3587, 3589, 3716),
    "very_composite": (630, 840, 1050, 1260)
  },
  {
    "prime": (3607, 4201, 4801, 5623, 6427), # primes
    "semi_1": (43 * 67, 53 * 71, 41 * 113, 79 * 83, 73 * 137), # pq, with p ~ q
    "semi_2": (17 * 223, 23 * 191, 13 * 521, 19 * 419, 11 * 1009), # pq, with p << q
    "comp_1": (4 * 1019, 6 * 727, 8 * 601, 9 * 631, 12 * 479), # np, with small composite n
    "comp_2": (44 * 79, 58 * 61, 68 * 67, 69 * 71, 86 * 59), # np, with n ~ p
    "comp_3": (17 * 17 * 17, 13 * 13 * 23, 13 * 17 * 19, 11 * 23 * 29, 11 * 17 * 41), # pqr, with p ~ q ~ r
    "very_comp": (1680, 1890, 2100, 2310, 2520) # only small prime divisors
  }
]

@dataclass
class SplitInfo:
  level: int
  splits: list[int]
  time: float

  @classmethod
  def parse_line(cls, line: str) -> 'SplitInfo':
    parts = line.split(":")
    level = int(parts[0])
    splits = [int(x) for x in parts[1][1:-1].split(',') if x]
    time = float(parts[2])
    return SplitInfo(level, splits, time)

@dataclass
class TableLine:
  category: str
  level: int
  ref_t: float
  ratio: float
  time: float
  adj_t: float
  speedup: float
  sub_times: tuple[float, float, float]
  splits: list[int]

  @classmethod
  def from_test_result(cls, category: str, ratio: float, ref_split: SplitInfo, test_split: SplitInfo, sub_times: list[float]) -> 'TableLine':
    return TableLine(
      category=category,
      level=test_split.level,
      ref_t=ref_split.time,
      ratio=ratio,
      time=test_split.time,
      adj_t=test_split.time * ratio,
      speedup=(ref_split.time / ratio) / test_split.time,
      sub_times=tuple(sub_times),
      splits=ref_split.splits
    )

  @classmethod
  def pretty_print_header(cls) -> None:
    print(f'{"category":<9} | {"level":>5} | {"ref_t":>6} | {"time":>6} {"adj_t":>6} | {"speedup":>8} | {"basis":>6} {"nspace":>6} {"decomp":>6} | {"b_rel":>6} {"n_rel":>6} {"d_rel":>6} | {"splits"}')

  def pretty_print(self) -> None:
    adj_ref_t = self.ref_t / self.ratio
    print(f'{self.category[:9]:<9} | {self.level:>5} | {self.ref_t:6.2f} | {self.time:6.2f} {self.adj_t:6.2f} | {self.speedup:8.4f} | {self.sub_times[0]:6.2f} {self.sub_times[1]:6.2f} {self.sub_times[2]:6.2f} | {self.sub_times[0]/adj_ref_t:6.2f} {self.sub_times[1]/adj_ref_t:6.2f} {self.sub_times[2]/adj_ref_t:6.2f} | {self.splits}')

def load_splits(filename: str) -> dict[int, SplitInfo]:
  with open(filename, "r") as f:
    split_infos = [SplitInfo.parse_line(line) for line in f.readlines()]
    return {si.level : si for si in split_infos}

def run_test(command: list[str], level: int) -> tuple[SplitInfo, list[float]]:
  command = command + ["-n", str(level)]
  result = subprocess.run(command, timeout=100, capture_output=True)
  lines = result.stdout.splitlines()
  last_line = lines[-1].decode()
  c_sub_times = [float(l.decode().split()[0][5:-1]) for l in lines[:-1]]
  sub_times = [c_sub_times[2], c_sub_times[3] - c_sub_times[2], c_sub_times[4] - c_sub_times[3]]
  return SplitInfo.parse_line(last_line), sub_times

def primesieve() -> float:
  primesieve_result = subprocess.run(
    ["primesieve", "-t", "1", "1e11", "-q", "--time"],
    capture_output=True
  )
  time_line = primesieve_result.stdout.splitlines()[0].decode()
  return float(time_line.split(' ')[1])

def perft(test_sets: list[int], command: list[str]) -> None:
  print("[info] Loading reference split_file")
  ref_splits = load_splits(SPLIT_FILE)

  print("[info] Running primesieve for calibration")
  primesieve_time = primesieve()
  ratio = SPLIT_FILE_CALIBRATE / primesieve_time
  print(f"[info] Primesieve finished in: {primesieve_time:7.4}s")
  print(f"[info] Machine performance ratio: {ratio:6.4}")

  for i in test_sets:
    print(Fore.YELLOW + f"\n[info] Running test set {i}\n" + Style.RESET_ALL)

    category_scores = []
    for k, v in TESTS[i].items():
      print(Fore.YELLOW + f"[info] Running tests for category {k}" + Style.RESET_ALL)
      scores = []

      TableLine.pretty_print_header()

      for level in v:
        # print(f"[info] Testing level {level}")
        test_split, sub_times = run_test(command, level)
        correct = test_split.splits == ref_splits[level].splits
        if correct:
          table_line = TableLine.from_test_result(k, ratio, ref_splits[level], test_split, sub_times)
          table_line.pretty_print()
          scores.append(table_line.speedup)
        else:
          print(Fore.RED + f"[result] Output incorrect, aborting perft" + Style.RESET_ALL)
          return

      avg_score = geometric_mean(scores)
      category_scores.append(avg_score)
      print(Fore.GREEN + f"[result] Finished category {k}, average speedup: {avg_score:6.4f}\n" + Style.RESET_ALL)

    print(Fore.GREEN + "[result] Finished all categories" + Style.RESET_ALL)

    avg_cat_score = geometric_mean(category_scores)
    print(Fore.GREEN + f"[result] Overall speedup for set {i}: {avg_cat_score:6.4f}" + Style.RESET_ALL)

if __name__ == "__main__":
  test_sets = [int(i) for i in sys.argv[1].split(',')]
  command = ["./bin/main", "-v", "1"] + sys.argv[2:]
  perft(test_sets, command)











