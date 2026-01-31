# -*- mode: python ; coding: utf-8 -*-
from PyInstaller.utils.hooks import collect_submodules, collect_data_files
from PyInstaller.utils.hooks import collect_all
import os

datas = []
binaries = []
hiddenimports = []

# Collect all necessary submodules
hiddenimports += collect_submodules('scipy')
hiddenimports += collect_submodules('numpy')
hiddenimports += collect_submodules('pandas')
hiddenimports += collect_submodules('matplotlib')
hiddenimports += collect_submodules('isotopylog')

# Collect all necessary data files
datas += collect_data_files('matplotlib')
datas += collect_data_files('isotopylog')

# Process matplotlib and isotopylog specifically
tmp_ret = collect_all('matplotlib')
datas += tmp_ret[0]; binaries += tmp_ret[1]; hiddenimports += tmp_ret[2]
tmp_ret = collect_all('isotopylog')
datas += tmp_ret[0]; binaries += tmp_ret[1]; hiddenimports += tmp_ret[2]

# Add clump_history package data
hiddenimports += collect_submodules('clump_history.src.clump_history')

block_cipher = None


a = Analysis(
    ['run_gui.py'],
    pathex=['../src'],
    binaries=binaries,
    datas=datas,
    hiddenimports=hiddenimports,
    hookspath=[],
    hooksconfig={},
    runtime_hooks=[],
    excludes=[],
    win_no_prefer_redirects=False,
    win_private_assemblies=False,
    cipher=block_cipher,
    noarchive=False,
)

# Exclude unnecessary modules to reduce size
exclude_modules = [
    'tkinter.test', 'tkinter.test.test_tkinter', 'tkinter.test.test_widgets',
    'matplotlib.tests', 'pandas.tests', 'numpy.f2py.tests', 'scipy.sparse.csgraph.tests',
    'scipy.spatial.distance.tests', 'scipy.special._precompute.testfunctions'
]
for mod in exclude_modules:
    try:
        if mod in a.scripts:
            a.scripts.remove(mod)
    except:
        pass

pyz = PYZ(a.pure, a.zipped_data, cipher=block_cipher)

exe = EXE(
    pyz,
    a.scripts,
    [],
    exclude_binaries=True,
    name='clump-history-gui',
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=True,
    console=False,
    disable_windowed_traceback=False,
    argv_emulation=False,
    target_arch=None,
    codesign_identity=None,
    entitlements_file=None,
    icon=None  # Add an icon if available
)

coll = COLLECT(
    exe,
    a.binaries,
    a.zipfiles,
    a.datas,
    strip=False,
    upx=True,
    upx_exclude=[],
    name='clump-history-gui',
)