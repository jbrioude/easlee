@echo off
rem ============================================================================
rem  EASLEE - one-click installer for Windows
rem  Double-click this file to install Julia (if needed) and all the libraries
rem  the model and the web GUI require.  Run it once.
rem ============================================================================
setlocal enabledelayedexpansion
cd /d "%~dp0"

echo ============================================================
echo    EASLEE - installing required libraries
echo ============================================================
echo.

rem ---- 1. Julia ---------------------------------------------------------------
where julia >nul 2>nul
if %errorlevel%==0 (
    echo [1/3] Julia is already installed:
    julia --version
) else (
    echo [1/3] Julia not found - installing with winget...
    winget install --id 9NJNWW8PVKMN -e --accept-source-agreements --accept-package-agreements
    if !errorlevel! neq 0 (
        echo.
        echo  Could not install Julia automatically.
        echo  Please install it manually from:  https://julialang.org/downloads/
        echo  then double-click this file again.
        echo.
        pause
        exit /b 1
    )
    echo.
    echo  Julia installed. Please CLOSE this window and double-click
    echo  install_easlee.bat again so Windows picks up the new PATH.
    echo.
    pause
    exit /b 0
)
echo.

rem ---- 2. model dependencies --------------------------------------------------
echo [2/3] Installing EASLEE model libraries (this can take 10-20 minutes)...
julia --project=easlee -e "import Pkg; Pkg.instantiate(); Pkg.precompile()"
if %errorlevel% neq 0 (
    echo  ERROR while installing the model libraries.
    pause
    exit /b 1
)
echo.

rem ---- 3. web GUI dependencies ------------------------------------------------
echo [3/3] Installing web GUI libraries...
julia --project=webgui -e "import Pkg; Pkg.instantiate(); Pkg.precompile()"
if %errorlevel% neq 0 (
    echo  ERROR while installing the web GUI libraries.
    pause
    exit /b 1
)

echo.
echo ============================================================
echo    Installation complete.
echo    Double-click  run_gui.bat  to start EASLEE.
echo ============================================================
echo.
pause
