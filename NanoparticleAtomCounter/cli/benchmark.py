"""
benchmark against the atomistic method
"""

from pathlib import Path
import sys
from ascii_colors import ASCIIColors
#import pytest
from atomcounter_benchmark.test_atomcount import main as benchmarking

def main() -> None:
    """
    Run the test_atomcount.py script, which benchmarks
    """
    print("This should take about 5 minutes,\n depending on how many cores your machine has . . .\n\n")
#    repo_root = Path(__file__).resolve().parents[2]
#    bnchk_dir = repo_root / "benchmark"
#    exit_code = pytest.main(["-s", str(bnchk_dir)])
    exit_code = benchmarking()
    text = "All benchmarking succeeded!" if exit_code == 0 else "Something went wrong! See .err files for details"
    color = ASCIIColors.color_green if exit_code == 0 else ASCIIColors.color_red
    ASCIIColors.print(
        text.upper(),
        color=color,
        style=ASCIIColors.style_bold,
        background=ASCIIColors.color_black,
        end="\n\n",
        flush=True,
        file=sys.stdout,
    )
    sys.exit(exit_code)  # sys.exit(pytest.main(["-s", str(tests_dir)]))


if __name__ == "__main__":
    main()
