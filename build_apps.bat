@echo off
REM ClumpIsotope Application Build Script
REM This script rebuilds both CLI and GUI applications using PyInstaller

echo Starting ClumpIsotope application build process...

REM Change to the clump_history directory
cd /d "%~dp0\clump_history"

echo Building GUI application...
pyinstaller clump-history-gui.spec --noconfirm --clean
if %ERRORLEVEL% NEQ 0 (
    echo Error building GUI application
    pause
    exit /b %ERRORLEVEL%
)

echo Building CLI application...
pyinstaller clump-history.spec --noconfirm --clean
if %ERRORLEVEL% NEQ 0 (
    echo Error building CLI application
    pause
    exit /b %ERRORLEVEL%
)

echo Build process completed successfully!
echo.
echo GUI application located at: dist\clump-history-gui\clump-history-gui.exe
echo CLI application located at: dist\clump-history\clump-history.exe
echo.
echo Press any key to exit...
pause > nul