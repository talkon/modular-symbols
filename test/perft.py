#!/usr/bin/env python

import subprocess
import sys
from dataclasses import dataclass
from statistics import geometric_mean
from colorama import Fore, Style

SPLIT_FILE = "./analysis/split.txt"
SPLIT_FILE_CALIBRATE = 13.4

BINARY = "./bin/main"
TESTS = {
  "prime": (1009, 1123, 1523, 2053, 2531, 3343, 3719, 4133, 4889),
  "composite": (1846, 3587, 3589, 3716),
  "very_composite": (630, 840, 1050, 1260)
}

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

def load_splits(filename: str) -> dict[int, SplitInfo]:
  with open(filename, "r") as f:
    split_infos = [SplitInfo.parse_line(line) for line in f.readlines()]
    return {si.level : si for si in split_infos}

def run_test(command: list[str], level: int) -> SplitInfo:
  command = command + ["-n", str(level)]
  result = subprocess.run(command, timeout=100, capture_output=True)
  last_line = result.stdout.splitlines()[-1].decode()
  return SplitInfo.parse_line(last_line)

def primesieve() -> float:
  primesieve_result = subprocess.run(
    ["primesieve", "-t", "1", "1e11", "-q", "--time"],
    capture_output=True
  )
  time_line = primesieve_result.stdout.splitlines()[0].decode()
  return float(time_line.split(' ')[1])

def perft(command: list[str]) -> None:
  print("[info] Loading reference split_file")
  ref_splits = load_splits(SPLIT_FILE)

  print("[info] Running primesieve for calibration")
  primesieve_time = primesieve()
  ratio = SPLIT_FILE_CALIBRATE / primesieve_time
  print(f"[info] Primesieve finished in: {primesieve_time:7.4}s")
  print(f"[info] Machine performance ratio: {ratio:6.4}\n")

  category_scores = []
  for k, v in TESTS.items():
    print(Fore.YELLOW + f"[info] Running tests for category {k}" + Style.RESET_ALL)
    scores = []

    for level in v:
      print(f"[info] Testing level {level}")
      test_split = run_test(command, level)
      correct = test_split.splits == ref_splits[level].splits
      if correct:
        rel_speed = ref_splits[level].time / (ratio * test_split.time)
        print(test_split.time, ratio, ref_splits[level].time)
        print(f"[result] Output correct, finished in {test_split.time:7.4f}s, speedup: {rel_speed:6.4f}")
        scores.append(rel_speed)
      else:
        print(Fore.RED + f"[result] Output incorrect, aborting perft" + Style.RESET_ALL)
        return

    avg_score = geometric_mean(scores)
    category_scores.append(avg_score)
    print(Fore.GREEN + f"[result] Finished category {k}, average speedup: {avg_score:6.4f}\n" + Style.RESET_ALL)

  print(Fore.GREEN + "[result] Finished all categories" + Style.RESET_ALL)

  avg_cat_score = geometric_mean(category_scores)
  print(Fore.GREEN + f"[result] Overall speedup: {avg_cat_score:6.4f}" + Style.RESET_ALL)

if __name__ == "__main__":
  command = ["./bin/main"] + sys.argv[1:]
  perft(command)











