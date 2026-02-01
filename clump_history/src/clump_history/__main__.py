"""
Entry point for running clump_history as a module.

Usage:
    python -m clump_history          # Run CLI by default
    python -m clump_history cli      # Run CLI
    python -m clump_history gui      # Run GUI
"""

import sys
from .cli import main as cli_main
from .gui import main as gui_main


def main():
    """Main entry point that routes to CLI or GUI based on arguments"""
    if len(sys.argv) > 1 and sys.argv[1] == 'gui':
        # Remove 'gui' from args and run GUI
        sys.argv.pop(1)
        gui_main()
    else:
        # Run CLI by default
        cli_main()


if __name__ == "__main__":
    main()
