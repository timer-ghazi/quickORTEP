#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# 
# Wed Feb 19 11:56:58 CST 2025
# - Added `BOHR`
# 
# Prev:
# - Refactored to use split dictionaries from data_split.py plus color data

"""
elements_table.py (refactored with color data)

Provides a class `Elements` with classmethods for lookups:

  - mass(symbol, unit="u")
  - vdw_radius(symbol, unit="Ang")
  - covalent_radius(symbol, order="single", source="cordero", unit="Ang")
  - electronegativity(symbol, scale="pauling")
  - atomic_number(symbol)
  - symbol(atomic_number)
  - name(symbol)
  - period(symbol)
  - group(symbol)
  - classification(symbol)
  - color(symbol, palette="Jmol")  <-- new method
  - is_valid(symbol)
  - list_symbols()

Data references:
  - ELEMENTS_CORE, VDW_RADII, BOND_PARAMS, ELECTRONEGATIVITY (from data_split.py)
  - COLORS (from color_data.py or data_split.py, depending on your setup)
"""

import unittest

# #############################################################################
# 1) Dictionaries of element data 
# #############################################################################
#
# Sources:
# 1. Van der Waals Radii:
#    - Source: https://en.wikipedia.org/wiki/Van_der_Waals_radius
#    - All values in Ångströms (Å)
#
# 2. Electronegativity:
#    - Three different scales:
#      * Pauling: dimensionless
#      * Mulliken: electron volts (eV)
#      * Allen: electron volts (eV)
#
# 3. Covalent Radii:
#    - Pyykko Parameters:
#      * Source: Pyykkö, P. & Atsumi, M.
#      * "Molecular Double-Bond Covalent Radii for Elements Li–E112"
#      * Chemistry: A European Journal 15 (46): 12770–12779
#      * DOI: dx.doi.org/10.1002/chem.200901472
#    - Cordero Parameters:
#      * Source: Cordero, B. et al.
#      * "Covalent radii revisited"
#      * Dalton Trans. (21): 2832–2838
#      * DOI: http://dx.doi.org/10.1039/b801115j
#    - All values in Ångströms (Å)
#
# 4. Atomic Masses:
#    - All values in unified atomic mass units (u)
#
# 5. Element Metadata:
#    - Includes: full names, atomic numbers, groups, periods, and classifications
#    - Classification categories: Nonmetal, Noble Gas, Alkali Metal,
#      Alkaline Earth Metal, Metalloid, Post-transition Metal,
#      Halogen, Transition Metal
#
# Data Organization:
# - Elements from Period 1-6
# - Main Group Elements followed by Transition Metals (Period 4 and 5)
# - All measurements in standard units (Å, eV, u)
#
# Created by combining multiple parameter files into a unified structure
# Last Updated: December 23, 2024
# Combined elements data from all source files

ELEMENTS_CORE = { 'H': { 'period': 1,
         'mass': 1.008,
         'atomic_number': 1,
         'group': 1,
         'classification': 'nonmetal',
         'name': 'Hydrogen'},
  'He': { 'period': 1,
          'mass': 4.002602,
          'atomic_number': 2,
          'group': 18,
          'classification': 'noble gas',
          'name': 'Helium'},
  'Li': { 'period': 2,
          'mass': 6.94,
          'atomic_number': 3,
          'group': 1,
          'classification': 'alkali metal',
          'name': 'Lithium'},
  'Be': { 'period': 2,
          'mass': 9.0122,
          'atomic_number': 4,
          'group': 2,
          'classification': 'alkaline earth metal',
          'name': 'Beryllium'},
  'B': { 'period': 2,
         'mass': 10.81,
         'atomic_number': 5,
         'group': 13,
         'classification': 'metalloid',
         'name': 'Boron'},
  'C': { 'period': 2,
         'mass': 12.011,
         'atomic_number': 6,
         'group': 14,
         'classification': 'nonmetal',
         'name': 'Carbon'},
  'N': { 'period': 2,
         'mass': 14.007,
         'atomic_number': 7,
         'group': 15,
         'classification': 'nonmetal',
         'name': 'Nitrogen'},
  'O': { 'period': 2,
         'mass': 15.999,
         'atomic_number': 8,
         'group': 16,
         'classification': 'nonmetal',
         'name': 'Oxygen'},
  'F': { 'period': 2,
         'mass': 18.998403163,
         'atomic_number': 9,
         'group': 17,
         'classification': 'halogen',
         'name': 'Fluorine'},
  'Ne': { 'period': 2,
          'mass': 20.1797,
          'atomic_number': 10,
          'group': 18,
          'classification': 'noble gas',
          'name': 'Neon'},
  'Na': { 'period': 3,
          'mass': 22.98976928,
          'atomic_number': 11,
          'group': 1,
          'classification': 'alkali metal',
          'name': 'Sodium'},
  'Mg': { 'period': 3,
          'mass': 24.305,
          'atomic_number': 12,
          'group': 2,
          'classification': 'alkaline earth metal',
          'name': 'Magnesium'},
  'Al': { 'period': 3,
          'mass': 26.9815385,
          'atomic_number': 13,
          'group': 13,
          'classification': 'post-transition metal',
          'name': 'Aluminum'},
  'Si': { 'period': 3,
          'mass': 28.085,
          'atomic_number': 14,
          'group': 14,
          'classification': 'metalloid',
          'name': 'Silicon'},
  'P': { 'period': 3,
         'mass': 30.973761998,
         'atomic_number': 15,
         'group': 15,
         'classification': 'nonmetal',
         'name': 'Phosphorus'},
  'S': { 'period': 3,
         'mass': 32.06,
         'atomic_number': 16,
         'group': 16,
         'classification': 'nonmetal',
         'name': 'Sulfur'},
  'Cl': { 'period': 3,
          'mass': 35.45,
          'atomic_number': 17,
          'group': 17,
          'classification': 'halogen',
          'name': 'Chlorine'},
  'Ar': { 'period': 3,
          'mass': 39.948,
          'atomic_number': 18,
          'group': 18,
          'classification': 'noble gas',
          'name': 'Argon'},
  'K': { 'period': 4,
         'mass': 39.0983,
         'atomic_number': 19,
         'group': 1,
         'classification': 'alkali metal',
         'name': 'Potassium'},
  'Ca': { 'period': 4,
          'mass': 40.078,
          'atomic_number': 20,
          'group': 2,
          'classification': 'alkaline earth metal',
          'name': 'Calcium'},
  'Sc': { 'period': 4,
          'mass': 44.955908,
          'atomic_number': 21,
          'group': 3,
          'classification': 'transition metal',
          'name': 'Scandium'},
  'Ti': { 'period': 4,
          'mass': 47.867,
          'atomic_number': 22,
          'group': 4,
          'classification': 'transition metal',
          'name': 'Titanium'},
  'V': { 'period': 4,
         'mass': 50.9415,
         'atomic_number': 23,
         'group': 5,
         'classification': 'transition metal',
         'name': 'Vanadium'},
  'Cr': { 'period': 4,
          'mass': 51.9961,
          'atomic_number': 24,
          'group': 6,
          'classification': 'transition metal',
          'name': 'Chromium'},
  'Mn': { 'period': 4,
          'mass': 54.938044,
          'atomic_number': 25,
          'group': 7,
          'classification': 'transition metal',
          'name': 'Manganese'},
  'Fe': { 'period': 4,
          'mass': 55.845,
          'atomic_number': 26,
          'group': 8,
          'classification': 'transition metal',
          'name': 'Iron'},
  'Co': { 'period': 4,
          'mass': 58.933194,
          'atomic_number': 27,
          'group': 9,
          'classification': 'transition metal',
          'name': 'Cobalt'},
  'Ni': { 'period': 4,
          'mass': 58.6934,
          'atomic_number': 28,
          'group': 10,
          'classification': 'transition metal',
          'name': 'Nickel'},
  'Cu': { 'period': 4,
          'mass': 63.546,
          'atomic_number': 29,
          'group': 11,
          'classification': 'transition metal',
          'name': 'Copper'},
  'Zn': { 'period': 4,
          'mass': 65.38,
          'atomic_number': 30,
          'group': 12,
          'classification': 'transition metal',
          'name': 'Zinc'},
  'Ga': { 'period': 4,
          'mass': 69.723,
          'atomic_number': 31,
          'group': 13,
          'classification': 'post-transition metal',
          'name': 'Gallium'},
  'Ge': { 'period': 4,
          'mass': 72.63,
          'atomic_number': 32,
          'group': 14,
          'classification': 'metalloid',
          'name': 'Germanium'},
  'As': { 'period': 4,
          'mass': 74.921595,
          'atomic_number': 33,
          'group': 15,
          'classification': 'metalloid',
          'name': 'Arsenic'},
  'Se': { 'period': 4,
          'mass': 78.971,
          'atomic_number': 34,
          'group': 16,
          'classification': 'nonmetal',
          'name': 'Selenium'},
  'Br': { 'period': 4,
          'mass': 79.904,
          'atomic_number': 35,
          'group': 17,
          'classification': 'halogen',
          'name': 'Bromine'},
  'Kr': { 'period': 4,
          'mass': 83.798,
          'atomic_number': 36,
          'group': 18,
          'classification': 'noble gas',
          'name': 'Krypton'},
  'Rb': { 'period': 5,
          'mass': 85.4678,
          'atomic_number': 37,
          'group': 1,
          'classification': 'alkali metal',
          'name': 'Rubidium'},
  'Sr': { 'period': 5,
          'mass': 87.62,
          'atomic_number': 38,
          'group': 2,
          'classification': 'alkaline earth metal',
          'name': 'Strontium'},
  'Y': { 'period': 5,
         'mass': 88.90584,
         'atomic_number': 39,
         'group': 3,
         'classification': 'transition metal',
         'name': 'Yttrium'},
  'Zr': { 'period': 5,
          'mass': 91.224,
          'atomic_number': 40,
          'group': 4,
          'classification': 'transition metal',
          'name': 'Zirconium'},
  'Nb': { 'period': 5,
          'mass': 92.90637,
          'atomic_number': 41,
          'group': 5,
          'classification': 'transition metal',
          'name': 'Niobium'},
  'Mo': { 'period': 5,
          'mass': 95.95,
          'atomic_number': 42,
          'group': 6,
          'classification': 'transition metal',
          'name': 'Molybdenum'},
  'Tc': { 'period': 5,
          'mass': 98.0,
          'atomic_number': 43,
          'group': 7,
          'classification': 'transition metal',
          'name': 'Technetium'},
  'Ru': { 'period': 5,
          'mass': 101.07,
          'atomic_number': 44,
          'group': 8,
          'classification': 'transition metal',
          'name': 'Ruthenium'},
  'Rh': { 'period': 5,
          'mass': 102.9055,
          'atomic_number': 45,
          'group': 9,
          'classification': 'transition metal',
          'name': 'Rhodium'},
  'Pd': { 'period': 5,
          'mass': 106.42,
          'atomic_number': 46,
          'group': 10,
          'classification': 'transition metal',
          'name': 'Palladium'},
  'Ag': { 'period': 5,
          'mass': 107.8682,
          'atomic_number': 47,
          'group': 11,
          'classification': 'transition metal',
          'name': 'Silver'},
  'Cd': { 'period': 5,
          'mass': 112.414,
          'atomic_number': 48,
          'group': 12,
          'classification': 'transition metal',
          'name': 'Cadmium'},
  'In': { 'period': 5,
          'mass': 114.818,
          'atomic_number': 49,
          'group': 13,
          'classification': 'post-transition metal',
          'name': 'Indium'},
  'Sn': { 'period': 5,
          'mass': 118.71,
          'atomic_number': 50,
          'group': 14,
          'classification': 'post-transition metal',
          'name': 'Tin'},
  'Sb': { 'period': 5,
          'mass': 121.76,
          'atomic_number': 51,
          'group': 15,
          'classification': 'metalloid',
          'name': 'Antimony'},
  'Te': { 'period': 5,
          'mass': 127.6,
          'atomic_number': 52,
          'group': 16,
          'classification': 'metalloid',
          'name': 'Tellurium'},
  'I': { 'period': 5,
         'mass': 126.90447,
         'atomic_number': 53,
         'group': 17,
         'classification': 'halogen',
         'name': 'Iodine'},
  'Xe': { 'period': 5,
          'mass': 131.293,
          'atomic_number': 54,
          'group': 18,
          'classification': 'noble gas',
          'name': 'Xenon'},
  'Cs': { 'period': 6,
          'mass': 132.90545196,
          'atomic_number': 55,
          'group': 1,
          'classification': 'alkali metal',
          'name': 'Cesium'},
  'Ba': { 'period': 6,
          'mass': 137.327,
          'atomic_number': 56,
          'group': 2,
          'classification': 'alkaline earth metal',
          'name': 'Barium'},
  'La': { 'period': 6,
          'mass': 138.905,
          'atomic_number': 57,
          'group': 3,
          'classification': 'lanthanoid',
          'name': 'Lanthanum'},
  'Ce': { 'period': 6,
          'mass': 140.116,
          'atomic_number': 58,
          'group': None,
          'classification': 'lanthanoid',
          'name': 'Cerium'},
  'Pr': { 'period': 6,
          'mass': 140.908,
          'atomic_number': 59,
          'group': None,
          'classification': 'lanthanoid',
          'name': 'Praseodymium'},
  'Nd': { 'period': 6,
          'mass': 144.24,
          'atomic_number': 60,
          'classification': 'lanthanoid',
          'name': 'Neodymium'},
  'Pm': { 'period': 6,
          'mass': 145.0,
          'atomic_number': 61,
          'group': None,
          'classification': 'lanthanoid',
          'name': 'Promethium'},
  'Sm': { 'period': 6,
          'mass': 150.36,
          'atomic_number': 62,
          'group': None,
          'classification': 'lanthanoid',
          'name': 'Samarium'},
  'Eu': { 'period': 6,
          'mass': 151.964,
          'atomic_number': 63,
          'group': None,
          'classification': 'lanthanoid',
          'name': 'Europium'},
  'Gd': { 'period': 6,
          'mass': 157.25,
          'atomic_number': 64,
          'group': None,
          'classification': 'lanthanoid',
          'name': 'Gadolinium'},
  'Tb': { 'period': 6,
          'mass': 158.925,
          'atomic_number': 65,
          'group': None,
          'classification': 'lanthanoid',
          'name': 'Terbium'},
  'Dy': { 'period': 6,
          'mass': 162.5,
          'atomic_number': 66,
          'group': None,
          'classification': 'lanthanoid',
          'name': 'Dysprosium'},
  'Ho': { 'period': 6,
          'mass': 164.93,
          'atomic_number': 67,
          'group': None,
          'classification': 'lanthanoid',
          'name': 'Holmium'},
  'Er': { 'period': 6,
          'mass': 167.259,
          'atomic_number': 68,
          'group': None,
          'classification': 'lanthanoid',
          'name': 'Erbium'},
  'Tm': { 'period': 6,
          'mass': 168.934,
          'atomic_number': 69,
          'group': None,
          'classification': 'lanthanoid',
          'name': 'Thulium'},
  'Yb': { 'period': 6,
          'mass': 173.054,
          'atomic_number': 70,
          'group': None,
          'classification': 'lanthanoid',
          'name': 'Ytterbium'},
  'Lu': { 'period': 6,
          'mass': 174.967,
          'atomic_number': 71,
          'group': 3,
          'classification': 'lanthanoid',
          'name': 'Lutetium'},
  'Hf': { 'period': 6,
          'mass': 178.49,
          'atomic_number': 72,
          'group': 4,
          'classification': 'transition metal',
          'name': 'Hafnium'},
  'Ta': { 'period': 6,
          'mass': 180.948,
          'atomic_number': 73,
          'group': 5,
          'classification': 'transition metal',
          'name': 'Tantalum'},
  'W': { 'period': 6,
         'mass': 183.84,
         'atomic_number': 74,
         'group': 6,
         'classification': 'transition metal',
         'name': 'Tungsten'},
  'Re': { 'period': 6,
          'mass': 186.207,
          'atomic_number': 75,
          'group': 7,
          'classification': 'transition metal',
          'name': 'Rhenium'},
  'Os': { 'period': 6,
          'mass': 190.23,
          'atomic_number': 76,
          'group': 8,
          'classification': 'transition metal',
          'name': 'Osmium'},
  'Ir': { 'period': 6,
          'mass': 192.217,
          'atomic_number': 77,
          'group': 9,
          'classification': 'transition metal',
          'name': 'Iridium'},
  'Pt': { 'period': 6,
          'mass': 195.084,
          'atomic_number': 78,
          'group': 10,
          'classification': 'transition metal',
          'name': 'Platinum'},
  'Au': { 'period': 6,
          'mass': 196.967,
          'atomic_number': 79,
          'group': 11,
          'classification': 'transition metal',
          'name': 'Gold'},
  'Hg': { 'period': 6,
          'mass': 200.59,
          'atomic_number': 80,
          'group': 12,
          'classification': 'transition metal',
          'name': 'Mercury'},
  'Tl': { 'period': 6,
          'mass': 204.383,
          'atomic_number': 81,
          'group': 13,
          'classification': 'post-transition metal',
          'name': 'Thallium'},
  'Pb': { 'period': 6,
          'mass': 207.2,
          'atomic_number': 82,
          'group': 14,
          'classification': 'post-transition metal',
          'name': 'Lead'},
  'Bi': { 'period': 6,
          'mass': 208.98,
          'atomic_number': 83,
          'group': 15,
          'classification': 'post-transition metal',
          'name': 'Bismuth'},
  'Po': { 'period': 6,
          'mass': 209.0,
          'atomic_number': 84,
          'group': 16,
          'classification': 'post-transition metal',
          'name': 'Polonium'},
  'At': { 'period': 6,
          'mass': 210.0,
          'atomic_number': 85,
          'group': 17,
          'classification': 'halogen',
          'name': 'Astatine'},
  'Rn': { 'period': 6,
          'mass': 222.0,
          'atomic_number': 86,
          'group': 18,
          'classification': 'noble gas',
          'name': 'Radon'},
  'Fr': { 'period': 7,
          'mass': 223.0,
          'atomic_number': 87,
          'group': 1,
          'classification': 'alkali metal',
          'name': 'Francium'},
  'Ra': { 'period': 7,
          'mass': 226.0,
          'atomic_number': 88,
          'group': 2,
          'classification': 'alkaline earth metal',
          'name': 'Radium'},
  'Ac': { 'period': 7,
          'mass': 227.0,
          'atomic_number': 89,
          'group': 3,
          'classification': 'actinoid',
          'name': 'Actinium'},
  'Th': { 'period': 7,
          'mass': 232.038,
          'atomic_number': 90,
          'group': None,
          'classification': 'actinoid',
          'name': 'Thorium'},
  'Pa': { 'period': 7,
          'mass': 231.036,
          'atomic_number': 91,
          'group': None,
          'classification': 'actinoid',
          'name': 'Protactinium'},
  'U': { 'period': 7,
         'mass': 238.029,
         'atomic_number': 92,
         'group': None,
         'classification': 'actinoid',
         'name': 'Uranium'},
  'Np': { 'period': 7,
          'mass': 237.0,
          'atomic_number': 93,
          'group': None,
          'classification': 'actinoid',
          'name': 'Neptunium'},
  'Pu': { 'period': 7,
          'mass': 244.0,
          'atomic_number': 94,
          'group': None,
          'classification': 'actinoid',
          'name': 'Plutonium'},
  'Am': { 'period': 7,
          'mass': 243.0,
          'atomic_number': 95,
          'group': None,
          'classification': 'actinoid',
          'name': 'Americium'},
  'Cm': { 'period': 7,
          'mass': 247.0,
          'atomic_number': 96,
          'group': None,
          'classification': 'actinoid',
          'name': 'Curium'},
  'Bk': { 'period': 7,
          'mass': 247.0,
          'atomic_number': 97,
          'group': None,
          'classification': 'actinoid',
          'name': 'Berkelium'},
  'Cf': { 'period': 7,
          'mass': 251.0,
          'atomic_number': 98,
          'group': None,
          'classification': 'actinoid',
          'name': 'Californium'},
  'Es': { 'period': 7,
          'mass': 252.0,
          'atomic_number': 99,
          'group': None,
          'classification': 'actinoid',
          'name': 'Einsteinium'},
  'Fm': { 'period': 7,
          'mass': 257.0,
          'atomic_number': 100,
          'group': None,
          'classification': 'actinoid',
          'name': 'Fermium'},
  'Md': { 'period': 7,
          'mass': 258.0,
          'atomic_number': 101,
          'group': None,
          'classification': 'actinoid',
          'name': 'Mendelevium'},
  'No': { 'period': 7,
          'mass': 259.0,
          'atomic_number': 102,
          'group': None,
          'classification': 'actinoid',
          'name': 'Nobelium'},
  'Lr': { 'period': 7,
          'mass': 262.0,
          'atomic_number': 103,
          'group': 3,
          'classification': 'actinoid',
          'name': 'Lawrencium'},
  'Rf': { 'period': 7,
          'mass': 267.0,
          'atomic_number': 104,
          'group': 4,
          'classification': 'transition metal',
          'name': 'Rutherfordium'},
  'Db': { 'period': 7,
          'mass': 268.0,
          'atomic_number': 105,
          'group': 5,
          'classification': 'transition metal',
          'name': 'Dubnium'},
  'Sg': { 'period': 7,
          'mass': 269.0,
          'atomic_number': 106,
          'group': 6,
          'classification': 'transition metal',
          'name': 'Seaborgium'},
  'Bh': { 'period': 7,
          'mass': 270.0,
          'atomic_number': 107,
          'group': 7,
          'classification': 'transition metal',
          'name': 'Bohrium'},
  'Hs': { 'period': 7,
          'mass': 269.0,
          'atomic_number': 108,
          'group': 8,
          'classification': 'transition metal',
          'name': 'Hassium'},
  'Mt': { 'period': 7,
          'mass': 278.0,
          'atomic_number': 109,
          'group': 9,
          'classification': 'transition metal',
          'name': 'Meitnerium'},
  'Ds': { 'period': 7,
          'mass': 281.0,
          'atomic_number': 110,
          'group': 10,
          'classification': 'transition metal',
          'name': 'Darmstadtium'},
  'Rg': { 'period': 7,
          'mass': 282.0,
          'atomic_number': 111,
          'group': 11,
          'classification': 'transition metal',
          'name': 'Roentgenium'},
  'Cn': { 'period': 7,
          'mass': 285.0,
          'atomic_number': 112,
          'group': 12,
          'classification': 'transition metal',
          'name': 'Copernicium'},
  'Nh': { 'period': 7,
          'mass': 286.0,
          'atomic_number': 113,
          'group': 13,
          'classification': 'post-transition metal',
          'name': 'Nihonium'},
  'Fl': { 'period': 7,
          'mass': 289.0,
          'atomic_number': 114,
          'group': 14,
          'classification': 'post-transition metal',
          'name': 'Flerovium'},
  'Mc': { 'period': 7,
          'mass': 290.0,
          'atomic_number': 115,
          'group': 15,
          'classification': 'post-transition metal',
          'name': 'Moscovium'},
  'Lv': { 'period': 7,
          'mass': 293.0,
          'atomic_number': 116,
          'group': 16,
          'classification': 'post-transition metal',
          'name': 'Livermorium'},
  'Ts': { 'period': 7,
          'mass': 294.0,
          'atomic_number': 117,
          'group': 17,
          'classification': 'halogen',
          'name': 'Tennessine'},
  'Og': { 'period': 7,
          'mass': 294.0,
          'atomic_number': 118,
          'group': 18,
          'classification': 'noble gas',
          'name': 'Oganesson'}}

VDW_RADII = { 'H': 1.2,
  'He': 1.4,
  'Li': 1.82,
  'Be': 1.53,
  'B': 1.92,
  'C': 1.7,
  'N': 1.55,
  'O': 1.52,
  'F': 1.47,
  'Ne': 1.54,
  'Na': 2.27,
  'Mg': 1.73,
  'Al': 1.84,
  'Si': 2.1,
  'P': 1.8,
  'S': 1.8,
  'Cl': 1.75,
  'Ar': 1.88,
  'K': 2.75,
  'Ca': 2.31,
  'Sc': 2.11,
  'Ti': 2.0,
  'V': 2.0,
  'Cr': 2.0,
  'Mn': 2.0,
  'Fe': 2.0,
  'Co': 2.0,
  'Ni': 1.63,
  'Cu': 1.4,
  'Zn': 1.39,
  'Ga': 1.87,
  'Ge': 2.11,
  'As': 1.85,
  'Se': 1.9,
  'Br': 1.85,
  'Kr': 2.02,
  'Rb': 3.03,
  'Sr': 2.49,
  'Y': 2.19,
  'Zr': 2.06,
  'Nb': 1.98,
  'Mo': 1.9,
  'Tc': 1.83,
  'Ru': 1.78,
  'Rh': 1.73,
  'Pd': 1.63,
  'Ag': 1.72,
  'Cd': 1.58,
  'In': 1.93,
  'Sn': 2.17,
  'Sb': 2.06,
  'Te': 2.06,
  'I': 1.98,
  'Cs': 3.43,
  'Ba': 2.68,
  'La': 2.4,
  'Ce': 2.35,
  'Pr': 2.39,
  'Nd': 2.29,
  'Pm': 2.36,
  'Sm': 2.29,
  'Eu': 2.33,
  'Gd': 2.37,
  'Tb': 2.21,
  'Dy': 2.29,
  'Ho': 2.16,
  'Er': 2.35,
  'Tm': 2.27,
  'Yb': 2.42,
  'Lu': 2.21}

BOND_PARAMS = { 'H': {'cordero': {'single': 0.31}, 'pyykko': {'single': 0.32}},
  'He': {'cordero': {'single': 0.28}, 'pyykko': {'single': 0.46}},
  'Li': { 'cordero': {'single': 1.28},
          'pyykko': {'single': 1.33, 'double': 1.24}},
  'Be': { 'cordero': {'single': 0.96},
          'pyykko': {'single': 1.02, 'double': 0.9, 'triple': 0.85}},
  'B': { 'cordero': {'single': 0.84},
         'pyykko': {'single': 0.85, 'double': 0.78, 'triple': 0.73}},
  'C': { 'cordero': {'single': 0.76},
         'pyykko': {'single': 0.75, 'double': 0.67, 'triple': 0.6}},
  'N': { 'cordero': {'single': 0.71},
         'pyykko': {'single': 0.71, 'double': 0.6, 'triple': 0.54}},
  'O': { 'cordero': {'single': 0.66},
         'pyykko': {'single': 0.63, 'double': 0.57, 'triple': 0.53}},
  'F': { 'cordero': {'single': 0.57},
         'pyykko': {'single': 0.64, 'double': 0.59, 'triple': 0.53}},
  'Ne': { 'cordero': {'single': 0.58},
          'pyykko': {'single': 0.67, 'double': 0.96}},
  'Na': { 'cordero': {'single': 1.66},
          'pyykko': {'single': 1.55, 'double': 1.6}},
  'Mg': { 'cordero': {'single': 1.41},
          'pyykko': {'single': 1.39, 'double': 1.32, 'triple': 1.27}},
  'Al': { 'cordero': {'single': 1.21},
          'pyykko': {'single': 1.26, 'double': 1.13, 'triple': 1.11}},
  'Si': { 'cordero': {'single': 1.11},
          'pyykko': {'single': 1.16, 'double': 1.07, 'triple': 1.02}},
  'P': { 'cordero': {'single': 1.07},
         'pyykko': {'single': 1.11, 'double': 1.02, 'triple': 0.94}},
  'S': { 'cordero': {'single': 1.05},
         'pyykko': {'single': 1.03, 'double': 0.94, 'triple': 0.95}},
  'Cl': { 'cordero': {'single': 1.02},
          'pyykko': {'single': 0.99, 'double': 0.95, 'triple': 0.93}},
  'Ar': { 'cordero': {'single': 1.06},
          'pyykko': {'single': 0.96, 'double': 1.07, 'triple': 0.96}},
  'K': { 'cordero': {'single': 2.03},
         'pyykko': {'single': 1.96, 'double': 1.93}},
  'Ca': { 'cordero': {'single': 1.76},
          'pyykko': {'single': 1.71, 'double': 1.47, 'triple': 1.33}},
  'Sc': { 'cordero': {'single': 1.7},
          'pyykko': {'single': 1.48, 'double': 1.16, 'triple': 1.14}},
  'Ti': { 'cordero': {'single': 1.6},
          'pyykko': {'single': 1.36, 'double': 1.17, 'triple': 1.08}},
  'V': { 'cordero': {'single': 1.53},
         'pyykko': {'single': 1.34, 'double': 1.12, 'triple': 1.06}},
  'Cr': { 'cordero': {'single': 1.39},
          'pyykko': {'single': 1.22, 'double': 1.11, 'triple': 1.03}},
  'Mn': { 'cordero': {'single': 1.5},
          'pyykko': {'single': 1.19, 'double': 1.05, 'triple': 1.03}},
  'Fe': { 'cordero': {'single': 1.42},
          'pyykko': {'single': 1.16, 'double': 1.09, 'triple': 1.02}},
  'Co': { 'cordero': {'single': 1.38},
          'pyykko': {'single': 1.11, 'double': 1.03, 'triple': 0.96}},
  'Ni': { 'cordero': {'single': 1.24},
          'pyykko': {'single': 1.1, 'double': 1.01, 'triple': 1.01}},
  'Cu': { 'cordero': {'single': 1.32},
          'pyykko': {'single': 1.12, 'double': 1.15, 'triple': 1.2}},
  'Zn': { 'cordero': {'single': 1.22},
          'pyykko': {'single': 1.18, 'double': 1.2}},
  'Ga': { 'cordero': {'single': 1.22},
          'pyykko': {'single': 1.24, 'double': 1.17, 'triple': 1.21}},
  'Ge': { 'cordero': {'single': 1.2},
          'pyykko': {'single': 1.21, 'double': 1.11, 'triple': 1.14}},
  'As': { 'cordero': {'single': 1.19},
          'pyykko': {'single': 1.21, 'double': 1.14, 'triple': 1.06}},
  'Se': { 'cordero': {'single': 1.2},
          'pyykko': {'single': 1.16, 'double': 1.07, 'triple': 1.07}},
  'Br': { 'cordero': {'single': 1.2},
          'pyykko': {'single': 1.14, 'double': 1.09, 'triple': 1.1}},
  'Kr': { 'cordero': {'single': 1.16},
          'pyykko': {'single': 1.17, 'double': 1.21, 'triple': 1.08}},
  'Rb': {'cordero': {'single': 2.2}, 'pyykko': {'single': 2.1, 'double': 2.02}},
  'Sr': { 'cordero': {'single': 1.95},
          'pyykko': {'single': 1.85, 'double': 1.57, 'triple': 1.39}},
  'Y': { 'cordero': {'single': 1.9},
         'pyykko': {'single': 1.63, 'double': 1.3, 'triple': 1.24}},
  'Zr': { 'cordero': {'single': 1.75},
          'pyykko': {'single': 1.54, 'double': 1.27, 'triple': 1.21}},
  'Nb': { 'cordero': {'single': 1.64},
          'pyykko': {'single': 1.47, 'double': 1.25, 'triple': 1.16}},
  'Mo': { 'cordero': {'single': 1.54},
          'pyykko': {'single': 1.38, 'double': 1.21, 'triple': 1.13}},
  'Tc': { 'cordero': {'single': 1.47},
          'pyykko': {'single': 1.28, 'double': 1.2, 'triple': 1.1}},
  'Ru': { 'cordero': {'single': 1.46},
          'pyykko': {'single': 1.25, 'double': 1.14, 'triple': 1.03}},
  'Rh': { 'cordero': {'single': 1.42},
          'pyykko': {'single': 1.25, 'double': 1.1, 'triple': 1.06}},
  'Pd': { 'cordero': {'single': 1.39},
          'pyykko': {'single': 1.2, 'double': 1.17, 'triple': 1.12}},
  'Ag': { 'cordero': {'single': 1.45},
          'pyykko': {'single': 1.28, 'double': 1.39, 'triple': 1.37}},
  'Cd': { 'cordero': {'single': 1.44},
          'pyykko': {'single': 1.36, 'double': 1.44}},
  'In': { 'cordero': {'single': 1.42},
          'pyykko': {'single': 1.42, 'double': 1.36, 'triple': 1.46}},
  'Sn': { 'cordero': {'single': 1.39},
          'pyykko': {'single': 1.4, 'double': 1.3, 'triple': 1.32}},
  'Sb': { 'cordero': {'single': 1.39},
          'pyykko': {'single': 1.4, 'double': 1.33, 'triple': 1.27}},
  'Te': { 'cordero': {'single': 1.38},
          'pyykko': {'single': 1.36, 'double': 1.28, 'triple': 1.21}},
  'I': { 'cordero': {'single': 1.39},
         'pyykko': {'single': 1.33, 'double': 1.29, 'triple': 1.25}},
  'Xe': { 'cordero': {'single': 1.4},
          'pyykko': {'single': 1.31, 'double': 1.35, 'triple': 1.22}},
  'Cs': { 'cordero': {'single': 2.44},
          'pyykko': {'single': 2.32, 'double': 2.09}},
  'Ba': { 'cordero': {'single': 2.15},
          'pyykko': {'single': 1.96, 'double': 1.61, 'triple': 1.49}},
  'La': { 'cordero': {'single': 2.07},
          'pyykko': {'single': 1.8, 'double': 1.39, 'triple': 1.39}},
  'Ce': { 'cordero': {'single': 2.04},
          'pyykko': {'single': 1.63, 'double': 1.37, 'triple': 1.31}},
  'Pr': { 'cordero': {'single': 2.03},
          'pyykko': {'single': 1.76, 'double': 1.38, 'triple': 1.28}},
  'Gd': { 'cordero': {'single': 1.96},
          'pyykko': {'single': 1.69, 'double': 1.35, 'triple': 1.32}},
  'Lu': { 'cordero': {'single': 1.87},
          'pyykko': {'single': 1.62, 'double': 1.31, 'triple': 1.31}},
  'Hf': { 'cordero': {'single': 1.75},
          'pyykko': {'single': 1.52, 'double': 1.28, 'triple': 1.22}},
  'Ta': { 'cordero': {'single': 1.7},
          'pyykko': {'single': 1.46, 'double': 1.26, 'triple': 1.19}},
  'W': { 'cordero': {'single': 1.62},
         'pyykko': {'single': 1.37, 'double': 1.2, 'triple': 1.15}},
  'Re': { 'cordero': {'single': 1.51},
          'pyykko': {'single': 1.31, 'double': 1.19, 'triple': 1.1}},
  'Os': { 'cordero': {'single': 1.44},
          'pyykko': {'single': 1.29, 'double': 1.16, 'triple': 1.09}},
  'Ir': { 'cordero': {'single': 1.41},
          'pyykko': {'single': 1.22, 'double': 1.15, 'triple': 1.07}},
  'Pt': { 'cordero': {'single': 1.36},
          'pyykko': {'single': 1.23, 'double': 1.12, 'triple': 1.1}},
  'Au': { 'cordero': {'single': 1.36},
          'pyykko': {'single': 1.24, 'double': 1.21, 'triple': 1.23}},
  'Tl': { 'cordero': {'single': 1.45},
          'pyykko': {'single': 1.44, 'double': 1.42, 'triple': 1.5}},
  'Pb': { 'cordero': {'single': 1.46},
          'pyykko': {'single': 1.44, 'double': 1.35, 'triple': 1.37}},
  'Bi': { 'cordero': {'single': 1.48},
          'pyykko': {'single': 1.51, 'double': 1.41, 'triple': 1.35}},
  'Po': { 'cordero': {'single': 1.4},
          'pyykko': {'single': 1.45, 'double': 1.35, 'triple': 1.29}},
  'At': { 'cordero': {'single': 1.5},
          'pyykko': {'single': 1.47, 'double': 1.38, 'triple': 1.38}},
  'Rn': { 'cordero': {'single': 1.5},
          'pyykko': {'single': 1.42, 'double': 1.45, 'triple': 1.33}},
  'Ra': { 'cordero': {'single': 2.21},
          'pyykko': {'single': 2.01, 'double': 1.73, 'triple': 1.59}},
  'Ac': { 'cordero': {'single': 2.15},
          'pyykko': {'single': 1.86, 'double': 1.53, 'triple': 1.4}},
  'Th': { 'cordero': {'single': 2.06},
          'pyykko': {'single': 1.75, 'double': 1.43, 'triple': 1.36}},
  'Pa': { 'cordero': {'single': 2.0},
          'pyykko': {'single': 1.69, 'double': 1.38, 'triple': 1.29}},
  'U': { 'cordero': {'single': 1.96},
         'pyykko': {'single': 1.7, 'double': 1.34, 'triple': 1.18}},
  'Np': { 'cordero': {'single': 1.9},
          'pyykko': {'single': 1.71, 'double': 1.36, 'triple': 1.16}}}

ELECTRONEGATIVITY = { 'H': {'pauling': 2.2, 'allen': 2.3},
  'He': {'allen': 4.16},
  'Li': {'pauling': 0.98, 'allen': 0.912},
  'Be': {'pauling': 1.57, 'allen': 1.576},
  'B': {'pauling': 2.04, 'allen': 2.051},
  'C': {'pauling': 2.55, 'allen': 2.544},
  'N': {'pauling': 3.04, 'allen': 3.066},
  'O': {'pauling': 3.44, 'allen': 3.61},
  'F': {'pauling': 3.98, 'allen': 4.193},
  'Ne': {'allen': 4.787},
  'Na': {'pauling': 0.93, 'allen': 0.869},
  'Mg': {'pauling': 1.31, 'allen': 1.293},
  'Al': {'pauling': 1.61, 'allen': 1.613},
  'Si': {'pauling': 1.9, 'allen': 1.916},
  'P': {'pauling': 2.19, 'allen': 2.253},
  'S': {'pauling': 2.58, 'allen': 2.589},
  'Cl': {'pauling': 3.16, 'allen': 2.869},
  'Ar': {'allen': 3.242},
  'K': {'pauling': 0.82, 'allen': 0.734},
  'Ca': {'pauling': 1.0, 'allen': 1.034},
  'Sc': {'pauling': 1.36, 'allen': 1.19},
  'Ti': {'pauling': 1.54, 'allen': 1.38},
  'V': {'pauling': 1.63, 'allen': 1.53},
  'Cr': {'pauling': 1.66, 'allen': 1.65},
  'Mn': {'pauling': 1.55, 'allen': 1.75},
  'Fe': {'pauling': 1.83, 'allen': 1.8},
  'Co': {'pauling': 1.88, 'allen': 1.84},
  'Ni': {'pauling': 1.91, 'allen': 1.88},
  'Cu': {'pauling': 1.9, 'allen': 1.85},
  'Zn': {'pauling': 1.65, 'allen': 1.588},
  'Ga': {'pauling': 1.81, 'allen': 1.756},
  'Ge': {'pauling': 2.01, 'allen': 1.994},
  'As': {'pauling': 2.18, 'allen': 2.211},
  'Se': {'pauling': 2.55, 'allen': 2.424},
  'Br': {'pauling': 2.96, 'allen': 2.685},
  'Kr': {'pauling': 3.0, 'allen': 2.966},
  'Rb': {'pauling': 0.82, 'allen': 0.706},
  'Sr': {'pauling': 0.95, 'allen': 0.963},
  'Y': {'pauling': 1.22, 'allen': 1.12},
  'Zr': {'pauling': 1.33, 'allen': 1.32},
  'Nb': {'pauling': 1.6, 'allen': 1.41},
  'Mo': {'pauling': 2.16, 'allen': 1.47},
  'Tc': {'pauling': 1.9, 'allen': 1.51},
  'Ru': {'pauling': 2.2, 'allen': 1.54},
  'Rh': {'pauling': 2.28, 'allen': 1.56},
  'Pd': {'pauling': 2.2, 'allen': 1.58},
  'Ag': {'pauling': 1.93, 'allen': 1.87},
  'Cd': {'pauling': 1.69, 'allen': 1.521},
  'In': {'pauling': 1.78, 'allen': 1.656},
  'Sn': {'pauling': 1.96, 'allen': 1.824},
  'Sb': {'pauling': 2.05, 'allen': 1.984},
  'Te': {'pauling': 2.1, 'allen': 2.158},
  'I': {'pauling': 2.66, 'allen': 2.359},
  'Xe': {'allen': 2.582},
  'Cs': {'pauling': 0.79, 'allen': 0.659},
  'Ba': {'pauling': 0.89, 'allen': 0.881},
  'La': {'pauling': 1.1},
  'Ce': {'pauling': 1.12},
  'Pr': {'pauling': 1.13},
  'Nd': {'pauling': 1.14},
  'Pm': {'pauling': 1.13},
  'Sm': {'pauling': 1.17},
  'Eu': {'pauling': 1.2},
  'Gd': {'pauling': 1.2},
  'Tb': {'pauling': 1.2},
  'Dy': {'pauling': 1.22},
  'Ho': {'pauling': 1.23},
  'Er': {'pauling': 1.24},
  'Tm': {'pauling': 1.25},
  'Yb': {'pauling': 1.1},
  'Lu': {'pauling': 1.27, 'allen': 1.09},
  'Hf': {'pauling': 1.3, 'allen': 1.16},
  'Ta': {'pauling': 1.5, 'allen': 1.34},
  'W': {'pauling': 2.36, 'allen': 1.47},
  'Re': {'pauling': 1.9, 'allen': 1.6},
  'Os': {'pauling': 2.2, 'allen': 1.65},
  'Ir': {'pauling': 2.2, 'allen': 1.68},
  'Pt': {'pauling': 2.28, 'allen': 1.72},
  'Au': {'pauling': 2.54, 'allen': 1.92},
  'Hg': {'pauling': 2.0, 'allen': 1.765},
  'Tl': {'pauling': 1.62, 'allen': 1.789},
  'Pb': {'pauling': 2.33, 'allen': 1.854},
  'Bi': {'pauling': 2.02, 'allen': 2.01},
  'Po': {'pauling': 2.0, 'allen': 2.19},
  'At': {'pauling': 2.2, 'allen': 2.39},
  'Rn': {'allen': 2.6},
  'Fr': {'pauling': 0.7, 'allen': 0.67},
  'Ra': {'pauling': 0.9, 'allen': 0.89},
  'Ac': {'pauling': 1.1},
  'Th': {'pauling': 1.3},
  'Pa': {'pauling': 1.5},
  'U': {'pauling': 1.38},
  'Np': {'pauling': 1.36},
  'Pu': {'pauling': 1.28}}

COLORS = { 'H': {'Corey': '#ffffff', 'Koltun': '#ffffff', 'Jmol': '#ffffff', 'Rasmol': '#ffffff', 'PubChem': '#638c8c'},
  'He': {'Jmol': '#d9ffff', 'Rasmol': '#ffc0cb', 'PubChem': '#D593A1'},
  'Li': {'Jmol': '#cc80ff', 'Rasmol': '#b22222', 'PubChem': '#D56632'},
  'Be': {'Jmol': '#c2ff00', 'Rasmol': '#ff1493', 'PubChem': '#D5BAD5'},
  'B': {'Jmol': '#ffb5b5', 'Rasmol': '#00ff00', 'PubChem': '#2AD52A'},
  'C': {'Corey': '#202020', 'Koltun': '#202020', 'Jmol': '#909090', 'Rasmol': '#c8c8c8', 'PubChem': '#274A4A'},
  'N': {'Corey': '#2060ff', 'Koltun': '#2060ff', 'Jmol': '#3050f8', 'Rasmol': '#8f8fff', 'PubChem': '#0000FF'},
  'O': {'Corey': '#ee2010', 'Koltun': '#ee2010', 'Jmol': '#ff0d0d', 'Rasmol': '#f00000', 'PubChem': '#FF0000'},
  'F': {'Koltun': '#00ff00', 'Jmol': '#90e050', 'Rasmol': '#daa520', 'PubChem': '#D52092'},
  'Ne': {'Jmol': '#b3e3f5', 'Rasmol': '#ff1493', 'PubChem': '#FF00FF'},
  'Na': {'Jmol': '#ab5cf2', 'Rasmol': '#0000ff', 'PubChem': '#0E73D5'},
  'Mg': {'Jmol': '#8aff00', 'Rasmol': '#228b22', 'PubChem': '#198C19'},
  'Al': {'Jmol': '#bfa6a6', 'Rasmol': '#808090', 'PubChem': '#838C8C'},
  'Si': {'Jmol': '#f0c8a0', 'Rasmol': '#daa520', 'PubChem': '#D59E13'},
  'P': {'Koltun': '#8020ff', 'Jmol': '#ff8000', 'Rasmol': '#ffa500', 'PubChem': '#D58600'},
  'S': {'Koltun': '#ffff00', 'Jmol': '#ffff30', 'Rasmol': '#ffc832', 'PubChem': '#D5D500'},
  'Cl': {'Koltun': '#00bb00', 'Jmol': '#1ff01f', 'Rasmol': '#00ff00', 'PubChem': '#2AD52A'},
  'Ar': {'Jmol': '#80d1e3', 'Rasmol': '#ff1493', 'PubChem': '#FF00FF'},
  'K': {'Jmol': '#8f40d4', 'Rasmol': '#ff1493', 'PubChem': '#D50575'},
  'Ca': {'Jmol': '#3dff00', 'Rasmol': '#808090', 'PubChem': '#838C8C'},
  'Sc': {'Jmol': '#e6e6e6', 'Rasmol': '#ff1493', 'PubChem': '#838C8C'},
  'Ti': {'Jmol': '#bfc2c7', 'Rasmol': '#808090', 'PubChem': '#838C8C'},
  'V': {'Jmol': '#a6a6ab', 'Rasmol': '#ff1493', 'PubChem': '#838C8C'},
  'Cr': {'Jmol': '#8a99c7', 'Rasmol': '#808090', 'PubChem': '#838C8C'},
  'Mn': {'Jmol': '#9c7ac7', 'Rasmol': '#808090', 'PubChem': '#838C8C'},
  'Fe': {'Koltun': '#d0d0d0', 'Jmol': '#e06633', 'Rasmol': '#ffa500', 'PubChem': '#FFA900'},
  'Co': {'Koltun': '#d0d0d0', 'Jmol': '#f090a0', 'Rasmol': '#ff1493', 'PubChem': '#838C8C'},
  'Ni': {'Koltun': '#d0d0d0', 'Jmol': '#50d050', 'Rasmol': '#a52a2a', 'PubChem': '#838C8C'},
  'Cu': {'Koltun': '#d0d0d0', 'Jmol': '#c88033', 'Rasmol': '#a52a2a', 'PubChem': '#838C8C'},
  'Zn': {'Jmol': '#7d80b0', 'Rasmol': '#a52a2a', 'PubChem': '#838C8C'},
  'Ga': {'Jmol': '#c28f8f', 'Rasmol': '#ff1493', 'PubChem': '#D5CD72'},
  'Ge': {'Jmol': '#668f8f', 'Rasmol': '#ff1493', 'PubChem': '#D5CD72'},
  'As': {'Jmol': '#bd80e3', 'Rasmol': '#ff1493', 'PubChem': '#D56632'},
  'Se': {'Jmol': '#ffa100', 'Rasmol': '#ff1493', 'PubChem': '#D5CD72'},
  'Br': {'Koltun': '#008800', 'Jmol': '#a62929', 'Rasmol': '#a52a2a', 'PubChem': '#D58639'},
  'Kr': {'Jmol': '#5cb8d1', 'Rasmol': '#ff1493', 'PubChem': '#FF00FF'},
  'Rb': {'Jmol': '#702eb0', 'Rasmol': '#ff1493', 'PubChem': '#FF00FF'},
  'Sr': {'Jmol': '#00ff00', 'Rasmol': '#ff1493', 'PubChem': '#FF0000'},
  'Y': {'Jmol': '#94ffff', 'Rasmol': '#ff1493', 'PubChem': '#838C8C'},
  'Zr': {'Jmol': '#94e0e0', 'Rasmol': '#ff1493', 'PubChem': '#838C8C'},
  'Nb': {'Jmol': '#73c2c9', 'Rasmol': '#ff1493', 'PubChem': '#838C8C'},
  'Mo': {'Jmol': '#54b5b5', 'Rasmol': '#ff1493', 'PubChem': '#838C8C'},
  'Tc': {'Jmol': '#3b9e9e', 'Rasmol': '#ff1493', 'PubChem': '#838C8C'},
  'Ru': {'Jmol': '#248f8f', 'Rasmol': '#ff1493', 'PubChem': '#D5B100'},
  'Rh': {'Jmol': '#0a7d8c', 'Rasmol': '#ff1493', 'PubChem': '#D5B100'},
  'Pd': {'Jmol': '#006985', 'Rasmol': '#ff1493', 'PubChem': '#D5B100'},
  'Ag': {'Jmol': '#c0c0c0', 'Rasmol': '#808090', 'PubChem': '#838C8C'},
  'Cd': {'Jmol': '#ffd98f', 'Rasmol': '#ff1493', 'PubChem': '#838C8C'},
  'In': {'Jmol': '#a67573', 'Rasmol': '#ff1493', 'PubChem': '#838C8C'},
  'Sn': {'Jmol': '#668080', 'Rasmol': '#ff1493', 'PubChem': '#838C8C'},
  'Sb': {'Jmol': '#9e63b5', 'Rasmol': '#ff1493', 'PubChem': '#838C8C'},
  'Te': {'Jmol': '#d47a00', 'Rasmol': '#ff1493', 'PubChem': '#D5CD72'},
  'I': {'Koltun': '#005500', 'Jmol': '#940094', 'Rasmol': '#a020f0', 'PubChem': '#FF00FF'},
  'Xe': {'Jmol': '#429eb0', 'Rasmol': '#ff1493', 'PubChem': '#FF00FF'},
  'Cs': {'Jmol': '#57178f', 'Rasmol': '#ff1493', 'PubChem': '#D598D5'},
  'Ba': {'Jmol': '#00c900', 'Rasmol': '#ffa500', 'PubChem': '#2AD52A'},
  'La': {'Jmol': '#70d4ff', 'Rasmol': '#ff1493', 'PubChem': '#00CCD5'},
  'Ce': {'Jmol': '#ffffc7', 'Rasmol': '#ff1493', 'PubChem': '#00CCD5'},
  'Pr': {'Jmol': '#d9ffc7', 'Rasmol': '#ff1493', 'PubChem': '#00CCD5'},
  'Nd': {'Jmol': '#c7ffc7', 'Rasmol': '#ff1493', 'PubChem': '#00CCD5'},
  'Pm': {'Jmol': '#a3ffc7', 'Rasmol': '#ff1493', 'PubChem': '#00CCD5'},
  'Sm': {'Jmol': '#8fffc7', 'Rasmol': '#ff1493', 'PubChem': '#00CCD5'},
  'Eu': {'Jmol': '#61ffc7', 'Rasmol': '#ff1493', 'PubChem': '#00CCD5'},
  'Gd': {'Jmol': '#45ffc7', 'Rasmol': '#ff1493', 'PubChem': '#00CCD5'},
  'Tb': {'Jmol': '#30ffc7', 'Rasmol': '#ff1493', 'PubChem': '#00CCD5'},
  'Dy': {'Jmol': '#1fffc7', 'Rasmol': '#ff1493', 'PubChem': '#00CCD5'},
  'Ho': {'Jmol': '#00ff9c', 'Rasmol': '#ff1493', 'PubChem': '#00CCD5'},
  'Er': {'Jmol': '#00e675', 'Rasmol': '#ff1493', 'PubChem': '#00CCD5'},
  'Tm': {'Jmol': '#00d452', 'Rasmol': '#ff1493', 'PubChem': '#00CCD5'},
  'Yb': {'Jmol': '#00bf38', 'Rasmol': '#ff1493', 'PubChem': '#00CCD5'},
  'Lu': {'Jmol': '#00ab24', 'Rasmol': '#ff1493', 'PubChem': '#00CCD5'},
  'Hf': {'Jmol': '#4dc2ff', 'Rasmol': '#ff1493', 'PubChem': '#00CCD5'},
  'Ta': {'Jmol': '#4da6ff', 'Rasmol': '#ff1493', 'PubChem': '#838C8C'},
  'W': {'Jmol': '#2194d6', 'Rasmol': '#ff1493', 'PubChem': '#838C8C'},
  'Re': {'Jmol': '#267dab', 'Rasmol': '#ff1493', 'PubChem': '#D5B100'},
  'Os': {'Jmol': '#266696', 'Rasmol': '#ff1493', 'PubChem': '#838C8C'},
  'Ir': {'Jmol': '#175487', 'Rasmol': '#ff1493', 'PubChem': '#D5B100'},
  'Pt': {'Jmol': '#d0d0e0', 'Rasmol': '#ff1493', 'PubChem': '#D5B100'},
  'Au': {'Jmol': '#ffd123', 'Rasmol': '#daa520', 'PubChem': '#D5B100'},
  'Hg': {'Jmol': '#b8b8d0', 'Rasmol': '#ff1493', 'PubChem': '#6E8092'},
  'Tl': {'Jmol': '#a6544d', 'Rasmol': '#ff1493', 'PubChem': '#D593A1'},
  'Pb': {'Jmol': '#575961', 'Rasmol': '#ff1493', 'PubChem': '#838C8C'},
  'Bi': {'Jmol': '#9e4fb5', 'Rasmol': '#ff1493', 'PubChem': '#D593A1'},
  'Po': {'Jmol': '#ab5c00', 'Rasmol': '#ff1493', 'PubChem': '#D593A1'},
  'At': {'Jmol': '#754f45', 'Rasmol': '#ff1493', 'PubChem': '#FF00FF'},
  'Rn': {'Jmol': '#428296', 'Rasmol': '#ff1493', 'PubChem': '#D598D5'},
  'Fr': {'Jmol': '#420066', 'Rasmol': '#ff1493', 'PubChem': '#D598D5'},
  'Ra': {'Jmol': '#007d00', 'Rasmol': '#ff1493', 'PubChem': '#2AD52A'},
  'Ac': {'Jmol': '#70abfa', 'Rasmol': '#ff1493', 'PubChem': '#00CCD5'},
  'Th': {'Jmol': '#00baff', 'Rasmol': '#ff1493', 'PubChem': '#00CCD5'},
  'Pa': {'Jmol': '#00a1ff', 'Rasmol': '#ff1493', 'PubChem': '#00CCD5'},
  'U': {'Jmol': '#008fff', 'Rasmol': '#ff1493', 'PubChem': '#00CCD5'},
  'Np': {'Jmol': '#0080ff', 'Rasmol': '#ff1493', 'PubChem': '#00CCD5'},
  'Pu': {'Jmol': '#006bff', 'Rasmol': '#ff1493', 'PubChem': '#00CCD5'},
  'Am': {'Jmol': '#545cf2', 'Rasmol': '#ff1493', 'PubChem': '#00CCD5'},
  'Cm': {'Jmol': '#785ce3', 'Rasmol': '#ff1493', 'PubChem': '#00CCD5'},
  'Bk': {'Jmol': '#8a4fe3', 'Rasmol': '#ff1493', 'PubChem': '#00CCD5'},
  'Cf': {'Jmol': '#a136d4', 'Rasmol': '#ff1493', 'PubChem': '#00CCD5'},
  'Es': {'Jmol': '#b31fd4', 'Rasmol': '#ff1493', 'PubChem': '#00CCD5'},
  'Fm': {'Jmol': '#b31fba', 'Rasmol': '#ff1493', 'PubChem': '#00CCD5'},
  'Md': {'Jmol': '#b30da6', 'Rasmol': '#ff1493', 'PubChem': '#00CCD5'},
  'No': {'Jmol': '#bd0d87', 'Rasmol': '#ff1493', 'PubChem': '#00CCD5'},
  'Lr': {'Jmol': '#c70066', 'Rasmol': '#ff1493', 'PubChem': '#00CCD5'},
  'Rf': {'Jmol': '#cc0059', 'Rasmol': '#ff1493', 'PubChem': '#00CCD5'},
  'Db': {'Jmol': '#d1004f', 'Rasmol': '#ff1493', 'PubChem': '#00CCD5'},
  'Sg': {'Jmol': '#d90045', 'Rasmol': '#ff1493', 'PubChem': '#00CCD5'},
  'Bh': {'Jmol': '#e00038', 'Rasmol': '#ff1493', 'PubChem': '#00CCD5'},
  'Hs': {'Jmol': '#e6002e', 'Rasmol': '#ff1493', 'PubChem': '#00CCD5'},
  'Mt': {'Jmol': '#eb0026', 'Rasmol': '#ff1493', 'PubChem': '#00CCD5'},
  'Ds': {'Rasmol': '#ff1493', 'PubChem': '#00CCD5'},
  'Rg': {'Rasmol': '#ff1493', 'PubChem': '#00CCD5'},
  'Cn': {'Rasmol': '#ff1493', 'PubChem': '#00CCD5'},
  'Nh': {'Rasmol': '#ff1493'},
  'Fl': {'Rasmol': '#ff1493'},
  'Mc': {'Rasmol': '#ff1493'},
  'Lv': {'Rasmol': '#ff1493'},
  'Ts': {'Rasmol': '#ff1493'},
  'Og': {'Rasmol': '#ff1493'}}

PROPERTY_REGISTRY = {
    'vdw_radius': VDW_RADII,
    'bond_params': BOND_PARAMS,
    'electronegativity': ELECTRONEGATIVITY,
    'color': COLORS 
}


# --------------------------------------------------------------------
# 2) Unit normalization dictionaries
# --------------------------------------------------------------------
_DISTANCE_UNIT_ALIASES = {
    "a": ("Ang", 1.0),
    "ang": ("Ang", 1.0),
    "angstrom": ("Ang", 1.0),
    "angstroms": ("Ang", 1.0),
    "ångström": ("Ang", 1.0),
    "å": ("Ang", 1.0),
    "Å": ("Ang", 1.0),

    "pm": ("pm", 100.0),
    "picometer": ("pm", 100.0),
    "picometre": ("pm", 100.0),

    "nm": ("nm", 0.1),
    "nanometer": ("nm", 0.1),
    "nanometre": ("nm", 0.1),

    "bohr": ("bohr", 1.889725989),
    "BOHR": ("bohr", 1.889725989),
    "a0": ("bohr", 1.889725989),
    "au": ("bohr", 1.889725989),
    "atomic_unit": ("bohr", 1.889725989),
}

_MASS_UNIT_ALIASES = {
    "u": ("u", 1.0),
    "amu": ("u", 1.0),
    "g/mol": ("g/mol", 1.0),
    "grams/mol": ("g/mol", 1.0),
}

# --------------------------------------------------------------------
# 3) Helpers for normalizing symbol and units
# --------------------------------------------------------------------
def _normalize_symbol(symbol: str) -> str:
    s = symbol.strip().capitalize()
    if s not in ELEMENTS_CORE:
        raise KeyError(f"Unknown element symbol: '{symbol}'")
    return s

def _normalize_distance_unit(unit: str):
    key = unit.strip().lower()
    if key not in _DISTANCE_UNIT_ALIASES:
        raise KeyError(f"Unknown distance unit: '{unit}'")
    return _DISTANCE_UNIT_ALIASES[key]

def _normalize_mass_unit(unit: str):
    key = unit.strip().lower()
    if key not in _MASS_UNIT_ALIASES:
        raise KeyError(f"Unknown mass unit: '{unit}'")
    return _MASS_UNIT_ALIASES[key]

def _get_symbol_by_atomic_number(atomic_num: int) -> str:
    for sym, data in ELEMENTS_CORE.items():
        if data["atomic_number"] == atomic_num:
            return sym
    raise KeyError(f"No element found for atomic number {atomic_num}")

# --------------------------------------------------------------------
# 4) Elements class
# --------------------------------------------------------------------
class Elements:
    @classmethod
    def is_valid(cls, symbol: str) -> bool:
        try:
            _normalize_symbol(symbol)
            return True
        except KeyError:
            return False

    @classmethod
    def list_symbols(cls):
        """Return a list of all element symbols, sorted by atomic number."""
        return sorted(ELEMENTS_CORE.keys(), key=lambda s: ELEMENTS_CORE[s]["atomic_number"])

    @classmethod
    def atomic_number(cls, symbol: str) -> int:
        sym = _normalize_symbol(symbol)
        return ELEMENTS_CORE[sym]["atomic_number"]

    @classmethod
    def symbol(cls, atomic_num: int) -> str:
        return _get_symbol_by_atomic_number(atomic_num)

    @classmethod
    def name(cls, symbol: str) -> str:
        sym = _normalize_symbol(symbol)
        return ELEMENTS_CORE[sym]["name"]

    @classmethod
    def mass(cls, symbol: str, unit: str = "u") -> float:
        sym = _normalize_symbol(symbol)
        mass_u = ELEMENTS_CORE[sym]["mass"]
        _, factor = _normalize_mass_unit(unit)
        return mass_u * factor

    @classmethod
    def vdw_radius(cls, symbol: str, unit: str = "Ang") -> float:
        sym = _normalize_symbol(symbol)
        val = VDW_RADII.get(sym)
        if val is None:
            raise KeyError(f"No van der Waals radius available for '{sym}'")
        _, factor = _normalize_distance_unit(unit)
        return val * factor

    @classmethod
    def covalent_radius(cls,
                        symbol: str,
                        order: str = "single",
                        source: str = "cordero",
                        unit: str = "Ang") -> float:
        sym = _normalize_symbol(symbol)
        bond_data = BOND_PARAMS.get(sym, {})
        src_data = bond_data.get(source)
        if src_data is None:
            raise KeyError(f"No covalent radius data for source='{source}' in '{sym}'")
        if order not in src_data:
            raise KeyError(
                f"No covalent radius for bond order='{order}' in source='{source}' for '{sym}'"
            )
        radius_A = src_data[order]
        _, factor = _normalize_distance_unit(unit)
        return radius_A * factor

    @classmethod
    def electronegativity(cls, symbol: str, scale: str = "pauling") -> float:
        sym = _normalize_symbol(symbol)
        e_data = ELECTRONEGATIVITY.get(sym, {})
        if scale not in e_data:
            raise KeyError(f"No electronegativity found for scale='{scale}' in '{sym}'")
        val = e_data[scale]
        if val is None:
            raise KeyError(f"No electronegativity value on '{scale}' scale for '{sym}' (data is None)")
        return val

    @classmethod
    def period(cls, symbol: str) -> int:
        sym = _normalize_symbol(symbol)
        return ELEMENTS_CORE[sym]["period"]

    @classmethod
    def group(cls, symbol: str) -> int:
        sym = _normalize_symbol(symbol)
        return ELEMENTS_CORE[sym]["group"]

    @classmethod
    def classification(cls, symbol: str) -> str:
        sym = _normalize_symbol(symbol)
        return ELEMENTS_CORE[sym]["classification"]

    # ----------------------------------------------------------------
    # *** NEW METHOD FOR COLOR LOOKUP ***
    # ----------------------------------------------------------------
    @classmethod
    def color(cls, symbol: str, palette: str = "Jmol") -> str:
        """
        Return the color (hex string, e.g. "#FFFFFF") for the given element symbol
        from the specified palette (Jmol, Rasmol, Corey, Koltun, PubChem, etc.).
        Raises KeyError if no color data is available for that palette.
        """
        sym = _normalize_symbol(symbol)
        color_dict = COLORS.get(sym, {})
        if palette not in color_dict:
            raise KeyError(f"No color defined for '{sym}' in palette='{palette}'")
        return color_dict[palette]

# --------------------------------------------------------------------
# 5) Tests (including a new one for color)
# --------------------------------------------------------------------
class TestElements(unittest.TestCase):
    def test_is_valid(self):
        self.assertTrue(Elements.is_valid("H"))
        self.assertTrue(Elements.is_valid("he"))
        self.assertFalse(Elements.is_valid("foo"))

    def test_atomic_number(self):
        self.assertEqual(Elements.atomic_number("H"), 1)
        self.assertEqual(Elements.atomic_number("he"), 2)
        with self.assertRaises(KeyError):
            Elements.atomic_number("fakeSymbol")

    def test_symbol(self):
        self.assertEqual(Elements.symbol(1), "H")
        self.assertEqual(Elements.symbol(2), "He")
        with self.assertRaises(KeyError):
            Elements.symbol(9999)

    def test_name(self):
        self.assertEqual(Elements.name("H"), "Hydrogen")
        self.assertEqual(Elements.name("hE"), "Helium")

    def test_mass(self):
        self.assertAlmostEqual(Elements.mass("H"), 1.008, places=3)
        self.assertAlmostEqual(Elements.mass("H", "g/mol"), 1.008, places=3)
        with self.assertRaises(KeyError):
            Elements.mass("H", "pounds")

    def test_vdw_radius(self):
        self.assertAlmostEqual(Elements.vdw_radius("H"), 1.2, places=2)
        self.assertAlmostEqual(Elements.vdw_radius("H", "pm"), 120.0, places=1)

    def test_covalent_radius(self):
        self.assertAlmostEqual(Elements.covalent_radius("H"), 0.31, places=2)
        self.assertAlmostEqual(Elements.covalent_radius("H", source="pyykko"), 0.32, places=2)
        with self.assertRaises(KeyError):
            Elements.covalent_radius("H", order="quadruple")

    def test_electronegativity(self):
        self.assertAlmostEqual(Elements.electronegativity("H", "pauling"), 2.2, places=2)
        with self.assertRaises(KeyError):
            Elements.electronegativity("H", "bogusScale")

    def test_period_and_group(self):
        self.assertEqual(Elements.period("Li"), 2)
        self.assertEqual(Elements.group("Li"), 1)
        self.assertEqual(Elements.period("He"), 1)
        self.assertEqual(Elements.group("He"), 18)

    def test_classification(self):
        self.assertEqual(Elements.classification("H"), "nonmetal")
        self.assertEqual(Elements.classification("he"), "noble gas")

    # *** NEW TEST for color lookup ***
    def test_color(self):
        # 'H' in Jmol palette -> '#ffffff' if you're using the snippet above
        self.assertEqual(Elements.color("H", "Jmol").lower(), "#ffffff")
        # If a user tries a nonexistent palette, we get KeyError
        with self.assertRaises(KeyError):
            Elements.color("H", "FooPalette")

if __name__ == "__main__":
    unittest.main()
