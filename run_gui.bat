@echo off
rem ============================================================================
rem  EASLEE - start the web GUI (Windows)
rem  Double-click this file.  Your browser opens at http://localhost:8080.
rem  Close this window to stop the model.
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
echo    Starting the EASLEE web GUI...
echo    A browser tab will open at http://localhost:8080
echo    Keep this window open while you use EASLEE.
echo    Close it to stop the model.
echo ============================================================
echo.

rem Use the pre-baked sysimage if present (fast start; build it once with
rem build_sysimage.bat).  Otherwise fall back to a plain launch.
if exist "easlee_sysimage.dll" (
    echo Using fast-start sysimage easlee_sysimage.dll
    julia --threads=auto --sysimage easlee_sysimage.dll --project=webgui webgui\easlee_app.jl
) else (
    julia --threads=auto --project=webgui webgui\easlee_app.jl
)

echo.
echo EASLEE has stopped.
pause
