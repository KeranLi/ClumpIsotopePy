#!/usr/bin/env python
"""
Test script to launch the Clump History GUI.
Usage: python test_gui.py
"""

import sys
import os

# Add src to path
sys.path.insert(0, os.path.join('clump_history', 'src'))

from clump_history.gui import main

if __name__ == "__main__":
    print("Starting Clump History GUI...")
    print("If the GUI window doesn't appear, check for errors above.")
    main()
