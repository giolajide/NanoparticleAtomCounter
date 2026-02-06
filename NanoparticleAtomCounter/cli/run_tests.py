"""
Automatically test all functions
"""

from pathlib import Path
import sys
from ascii_colors import ASCIIColors
import pytest


def main() -> None:
    """
    Run every test, then exit with pytest’s status code.
    """
    repo_root = Path(__file__).resolve().parents[2]
    tests_dir = repo_root / "atomcounter_tests"
    exit_code = pytest.main(["-s", str(tests_dir)])
    text = "All tests passed!" if exit_code == 0 else "Some tests failed!"
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
