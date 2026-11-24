"""Summarize paired ORFik timings and update the saved benchmark report."""
import csv
import random
import statistics
from pathlib import Path

root = Path("results/orfik_performance_2026-10-09")
source = root / "source_verification"
groups = {}
with (source / "timings.csv").open() as handle:
    for row in csv.DictReader(handle):
        key = row["function_name"], row["input"]
        groups.setdefault(key, {}).setdefault(int(row["iteration"]), {})[
            row["method"]
        ] = float(row["seconds"])

rng = random.Random(421)
summary = []
for (function, label), pairs in groups.items():
    before = [p["baseline"] for p in pairs.values()]
    after = [p["candidate"] for p in pairs.values()]
    speedup = statistics.median(before) / statistics.median(after)
    bootstrap = []
    for _ in range(5000):
        indices = rng.choices(range(len(before)), k=len(before))
        bootstrap.append(
            statistics.median(before[i] for i in indices)
            / statistics.median(after[i] for i in indices)
        )
    bootstrap.sort()
    summary.append(dict(
        function=function, input=label, paired_repetitions=len(before),
        before_seconds=statistics.median(before),
        after_seconds=statistics.median(after), speedup=speedup,
        time_saved_percent=100 * (1 - 1 / speedup),
        faster_pairs=sum(b > a for b, a in zip(before, after)),
        speedup_bootstrap_low=bootstrap[125],
        speedup_bootstrap_high=bootstrap[4874],
        before_min=min(before), before_max=max(before),
        after_min=min(after), after_max=max(after),
    ))
with (root / "final_paired_summary.csv").open("w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=list(summary[0]))
    writer.writeheader()
    writer.writerows(summary)

lines = [
    "Edited-source results; seconds per call. The four-library batch uses the",
    "same four actual columns for both implementations.", "",
    "| Function | Input | Before (s) | After (s) | Speedup | Paired repetitions |",
    "| --- | --- | ---: | ---: | ---: | ---: |",
]
for row in summary:
    lines.append(
        f"| {row['function']} | {row['input']} | {row['before_seconds']:.5f} | "
        f"{row['after_seconds']:.5f} | {row['speedup']:.2f}× | "
        f"{row['paired_repetitions']} |"
    )
lines.extend([
    "", "`final_paired_summary.csv` additionally retains observed ranges, the number",
    "of faster paired trials, and a paired bootstrap interval for the median",
    "speedup (5,000 resamples, fixed seed). These intervals describe this shared",
    "server run; they do not capture uncertainty across machines or cold disks.",
])
report = Path("md/orfik_function_performance_2026-10-09.md")
text = report.read_text()
begin = "<!-- BEGIN GENERATED TIMINGS -->"
end = "<!-- END GENERATED TIMINGS -->"
prefix, rest = text.split(begin, 1)
_, suffix = rest.split(end, 1)
report.write_text(prefix + begin + "\n" + "\n".join(lines) + "\n" + end + suffix)
for row in summary:
    print(f"{row['input']}: {row['before_seconds']:.5f} -> "
          f"{row['after_seconds']:.5f} s, {row['speedup']:.2f}x")
