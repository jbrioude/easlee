@echo off
rem ============================================================================
rem  EASLEE - build a fast-start sysimage (Windows)
rem  Run this ONCE, after install_easlee.bat.  It compiles the model and GUI
rem  libraries into  easlee_sysimage.dll  so the GUI then starts in seconds.
rem  Takes 10-25 minutes.  Double-click and wait.
rem ============================================================================
setlocal
cd /d "%~dp0"

where julia >nul 2>nul
if %errorlevel% neq 0 (
    echo Julia is not installed yet.
    echo Please double-click  install_easlee.bat  first.
    echo.
    pause
    exit /b 1
)

echo ============================================================
echo    Building the EASLEE sysimage (10-25 min, one time)...
echo    You can keep using your computer while this runs.
echo ============================================================
echo.

julia --project=webgui webgui\build_sysimage.jl
if %errorlevel% neq 0 (
    echo.
    echo  Sysimage build FAILED.  The GUI still works without it,
    echo  just slower to start.  See WINDOWS.md for help.
    echo.
    pause
    exit /b 1
)

echo.
echo ============================================================
echo    Done.  Double-click  run_gui.bat  (now starts fast),
echo    or  EASLEE.vbs  to launch with no console window.
echo ============================================================
echo.
pause
