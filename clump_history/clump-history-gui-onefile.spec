# -*- mode: python ; coding: utf-8 -*-
# Single-file executable configuration
# Use this for simpler distribution (one .exe file)

from PyInstaller.utils.hooks import collect_submodules, collect_data_files
from PyInstaller.utils.hooks import collect_all

datas = []
binaries = []
hiddenimports = []

# Collect all necessary submodules
hiddenimports += collect_submodules('scipy')
hiddenimports += collect_submodules('numpy')
hiddenimports += collect_submodules('pandas')
hiddenimports += collect_submodules('matplotlib')
hiddenimports += collect_submodules('isotopylog')
hiddenimports += [
    'matplotlib.backends.backend_tkagg',
    'matplotlib.backends.backend_agg',
]

# Collect data files
datas += collect_data_files('matplotlib')
datas += collect_data_files('isotopylog')

tmp_ret = collect_all('matplotlib')
datas += tmp_ret[0]; binaries += tmp_ret[1]; hiddenimports += tmp_ret[2]
tmp_ret = collect_all('isotopylog')
datas += tmp_ret[0]; binaries += tmp_ret[1]; hiddenimports += tmp_ret[2]

# Add clump_history package
hiddenimports += [
    'clump_history',
    'clump_history.model',
    'clump_history.fit',
    'clump_history.io',
    'clump_history.plot',
    'clump_history.cli',
    'clump_history.gui',
]

block_cipher = None

a = Analysis(
    ['run_gui.py'],
    pathex=['src'],
    binaries=binaries,
    datas=datas,
    hiddenimports=hiddenimports,
    hookspath=[],
    hooksconfig={},
    runtime_hooks=[],
    excludes=[
        'tkinter.test',
        'matplotlib.tests',
        'pandas.tests',
        'numpy.f2py.tests',
        'scipy.tests',
        'unittest',
        'pytest',
        '_pytest',
        'pylint',
        'mypy',
    ],
    win_no_prefer_redirects=False,
    win_private_assemblies=False,
    cipher=block_cipher,
    noarchive=False,
)

pyz = PYZ(a.pure, a.zipped_data, cipher=block_cipher)

# Single-file executable (onefile mode)
exe = EXE(
    pyz,
    a.scripts,
    a.binaries,
    a.zipfiles,
    a.datas,
    [],
    name='ClumpHistoryGUI',
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=True,
    upx_exclude=[],
    runtime_tmpdir=None,
    console=False,
    disable_windowed_traceback=False,
    argv_emulation=False,
    target_arch=None,
    codesign_identity=None,
    entitlements_file=None,
    icon=None,
)
