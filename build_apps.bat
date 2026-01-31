@echo off
echo Building ClumpIsotope applications...

REM 设置环境变量
set PYTHONPATH=%~dp0clump_history\src;%PYTHONPATH%

REM 检查是否安装了pyinstaller
python -c "import PyInstaller" 2>nul
if errorlevel 1 (
    echo Installing PyInstaller...
    pip install pyinstaller
)

echo Building CLI application...
pyinstaller clump-history.spec || echo Error building CLI application

echo Building GUI application...
pyinstaller clump-history-gui.spec || echo Error building GUI application

echo Build process completed. Check the dist/ folder for executables.
pause