"""
ClumpIsotope: Carbonate Clumped Isotope Geochemistry Modeling Framework

A comprehensive Python-based computational framework for carbonate clumped isotope 
(Δ47) forward/reordering calculations. The software implements solid-state isotope 
exchange kinetic models to simulate the temperature-dependent reordering of clumped 
isotopes in geological samples over time.

This enables reconstruction of thermal histories and paleotemperature estimates 
based on theoretical models from Stolper et al. (2015) and Hemingway & Henkes (2021).
"""

__all__ = ["__version__", "model", "fit", "plot", "io", "cli", "gui"]
__version__ = "0.1.1"