@echo off
chcp 65001 >nul
REM One-file build script for easy distribution

echo ============================================
echo Building Single-File Executable
echo ============================================
echo.

REM Check conda
call conda activate clump-isotope 2>nul
if %errorlevel% neq 0 (
    echo [WARNING] Using current Python environment
)

echo [Step 1/4] Installing dependencies...
pip install pyinstaller -q

echo.
echo [Step 2/4] Cleaning old files...
cd clump_history
if exist "dist" rmdir /s /q "dist"
if exist "build" rmdir /s /q "build"

echo.
echo [Step 3/4] Building single-file executable...
echo This will take 5-10 minutes...
pyinstaller clump-history-gui-onefile.spec --noconfirm --clean
if %errorlevel% neq 0 (
    echo [ERROR] Build failed!
    pause
    exit /b 1
)

echo.
echo [Step 4/4] Creating release package...
set "RELEASE_DIR=..\release"
if not exist "%RELEASE_DIR%" mkdir "%RELEASE_DIR%"

REM Copy executable
copy "dist\ClumpHistoryGUI.exe" "%RELEASE_DIR%\" >nul

REM Create sample data folder
if not exist "%RELEASE_DIR%\datasets" mkdir "%RELEASE_DIR%\datasets"
copy "..\datasets\Thermal_History_Hu.csv" "%RELEASE_DIR%\datasets\" >nul
copy "..\datasets\acutal_test_Hu.csv" "%RELEASE_DIR%\datasets\" >nul

REM Copy documentation
copy "..\RELEASE_NOTES.md" "%RELEASE_DIR%\README.txt" >nul
copy "..\QUICKSTART_GUI.md" "%RELEASE_DIR%\使用说明.txt" >nul

REM Create start script
(
echo @echo off
echo echo ========================================
echo echo  Clump History GUI v0.2.0
echo echo ========================================
echo echo.
echo echo Starting program...
echo echo 首次启动可能需要1-2分钟解压资源
echo echo.
echo start "" "ClumpHistoryGUI.exe"
echo timeout /t 3
echo exit
) > "%RELEASE_DIR%\启动程序.bat"

echo.
echo ============================================
echo Build completed successfully!
echo ============================================
echo.
echo Output: release\ClumpHistoryGUI.exe
echo.
echo To distribute:
echo   Zip the 'release' folder and share
echo.
cd ..
pause
