@echo off
chcp 65001 >nul
REM Build script for Clump History GUI Executable
REM Author: Keran Li
REM Date: 2026-02-01

echo ============================================
echo Clump History GUI - Build Script
echo ============================================
echo.

REM Check if conda is available
where conda >nul 2>nul
if %errorlevel% neq 0 (
    echo [ERROR] Conda not found. Please install Anaconda or Miniconda.
    exit /b 1
)

echo [1/8] Activating conda environment...
call conda activate clump-isotope
if %errorlevel% neq 0 (
    echo [WARNING] clump-isotope environment not found. Using current environment.
    echo To create the environment, run: conda env create -f create_env.yaml
)

echo.
echo [2/8] Cleaning old build files...
cd clump_history
if exist "dist" (
    rmdir /s /q "dist"
    echo      Removed old dist directory
)
if exist "build" (
    rmdir /s /q "build"
    echo      Removed old build directory
)

echo.
echo [3/8] Installing pyinstaller...
pip install pyinstaller -q

echo.
echo [4/8] Building GUI executable...
echo      This may take 5-10 minutes...
pyinstaller clump-history-gui.spec --noconfirm --clean
if %errorlevel% neq 0 (
    echo [ERROR] Build failed!
    exit /b 1
)

echo.
echo [5/8] Creating distribution package...
set "DIST_DIR=dist\ClumpHistoryGUI"
if not exist "%DIST_DIR%" mkdir "%DIST_DIR%"

REM Copy executable
xcopy /s /e /i /y "dist\clump-history-gui\*" "%DIST_DIR%\" >nul

echo.
echo [6/8] Copying sample data...
if not exist "%DIST_DIR%\datasets" mkdir "%DIST_DIR%\datasets"
xcopy /y "..\datasets\Thermal_History_Hu.csv" "%DIST_DIR%\datasets\" >nul
xcopy /y "..\datasets\acutal_test_Hu.csv" "%DIST_DIR%\datasets\" >nul
xcopy /y "..\datasets\Thermal_History.csv" "%DIST_DIR%\datasets\" >nul
xcopy /y "..\datasets\acutal_test.csv" "%DIST_DIR%\datasets\" >nul

echo.
echo [7/8] Creating user documentation...
if not exist "%DIST_DIR%\docs" mkdir "%DIST_DIR%\docs"
xcopy /y "..\README.md" "%DIST_DIR%\docs\" >nul
xcopy /y "..\QUICKSTART_GUI.md" "%DIST_DIR%\docs\" >nul
xcopy /y "..\GUI_USAGE.md" "%DIST_DIR%\docs\" >nul 2>nul
xcopy /y "..\THEORY.md" "%DIST_DIR%\docs\" >nul 2>nul

echo.
echo [8/8] Creating start script...
(
echo @echo off
echo chcp 65001 ^>nul
echo echo Starting Clump History GUI...
echo echo.
echo echo 如果窗口无法显示完整，请最大化窗口
echo.
echo start "" "%%~dp0clump-history-gui.exe"
) > "%DIST_DIR%\启动程序.bat"

echo.
echo ============================================
echo Build completed successfully!
echo ============================================
echo.
echo Output location: clump_history\dist\ClumpHistoryGUI\
echo.
echo To distribute:
echo   1. Zip the folder: clump_history\dist\ClumpHistoryGUI
echo   2. Share the zip file with users
echo   3. Users can run: 启动程序.bat
echo.
pause
