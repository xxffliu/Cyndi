@echo off
rem  Everything below is mirrored into build_log.txt, so the result can be read
rem  afterwards without scrolling a console window: the script re-invokes itself
rem  once with the output redirected, then prints the log.
if "%CYNDI_LOGGED%"=="1" goto :main
set "CYNDI_LOGGED=1"
cd /d "%~dp0"
call "%~f0" %* > "%~dp0build_log.txt" 2>&1
type "%~dp0build_log.txt"
set "CYNDI_LOGGED="
if not "%CYNDI_NOPAUSE%"=="1" pause
exit /b

:main
setlocal enabledelayedexpansion
rem ===========================================================================
rem  Cyndi -- one-shot CMake build on Windows with the MinGW-w64 toolchain.
rem
rem  Double-click this file, or run it from any shell:
rem      build_windows.bat            Release build + ctest
rem      build_windows.bat debug      Debug build
rem      build_windows.bat linear     Release build, pre-v1 linear torsion encoding
rem
rem  It installs cmake and ninja with winget if they are not already on PATH,
rem  then configures, builds and tests through the presets in CMakePresets.json.
rem  Nothing here touches the source tree.
rem ===========================================================================

cd /d "%~dp0"

set "PRESET=mingw-release"
if /i "%~1"=="debug"  set "PRESET=mingw-debug"
if /i "%~1"=="linear" set "PRESET=mingw-torsion-linear"

echo.
echo === Cyndi CMake build (preset: %PRESET%) ===
echo.

rem --- the compiler ----------------------------------------------------------
rem  Looked for in this order:
rem    1. %CYNDI_MINGW%       -- set this to override
rem    2. g++ already on PATH
rem    3. %USERPROFILE%\tools\mingw64 and a couple of common install locations
if defined CYNDI_MINGW if exist "%CYNDI_MINGW%\bin\g++.exe" set "MINGW=%CYNDI_MINGW%"
if not defined MINGW (
    for /f "delims=" %%p in ('where g++ 2^>nul') do if not defined MINGW (
        for %%d in ("%%~dpp..") do set "MINGW=%%~fd"
    )
)
if not defined MINGW (
    for %%c in ("%USERPROFILE%\tools\mingw64" "C:\msys64\mingw64" "C:\mingw64") do (
        if not defined MINGW if exist "%%~c\bin\g++.exe" set "MINGW=%%~c"
    )
)
if not defined MINGW (
    echo ERROR: no MinGW-w64 g++ found.
    echo        Put its bin directory on PATH, or set CYNDI_MINGW to the
    echo        toolchain root, e.g.:
    echo            set CYNDI_MINGW=C:\msys64\mingw64
    goto :fail
)
set "PATH=%MINGW%\bin;%PATH%"
for /f "tokens=*" %%v in ('g++ -dumpversion 2^>nul') do set "GXXVER=%%v"
echo   g++      %GXXVER%   (%MINGW%\bin)

rem --- cmake and ninja -------------------------------------------------------
rem winget updates the PATH of NEW shells only, so also add the two standard
rem install locations to this session's PATH after installing.
call :ensure cmake Kitware.CMake         || goto :fail
call :ensure ninja Ninja-build.Ninja     || goto :fail

for /f "tokens=*" %%v in ('cmake --version 2^>nul ^| findstr /r "^cmake"') do set "CMAKEVER=%%v"
echo   %CMAKEVER%
for /f "tokens=*" %%v in ('ninja --version 2^>nul') do set "NINJAVER=%%v"
echo   ninja    %NINJAVER%
echo.

rem --- configure, build, test ------------------------------------------------
echo --- configure ---
cmake --preset %PRESET%
if errorlevel 1 goto :fail

echo.
echo --- build ---
cmake --build --preset %PRESET%
if errorlevel 1 goto :fail

echo.
echo --- test ---
rem Only the Release presets have a matching test preset; Debug/linear just build.
if /i "%PRESET%"=="mingw-release" (
    ctest --preset mingw-release
    if errorlevel 1 goto :fail
)

echo.
echo === OK ===
echo Executable and the two .parm files are in the build directory for this
echo preset (build-mingw for Release). Run it from there, or from any folder
echo that has MMFF94.parm and TAFF.parm in it:
echo.
echo     build-mingw\cyndi.exe -input Cyndi\3ert.mol2 -output out.mol2 -parm parameter\CyndiParam.in
echo.
goto :done

rem ---------------------------------------------------------------------------
:ensure
rem  %1 = executable name, %2 = winget package id
where %1 >nul 2>nul && exit /b 0
echo   %1 not found -- installing %2 with winget...
where winget >nul 2>nul
if errorlevel 1 (
    echo.
    echo ERROR: %1 is missing and winget is not available either.
    echo        Install %1 by hand, or grab the portable zip:
    echo            cmake  https://cmake.org/download/
    echo            ninja  https://github.com/ninja-build/ninja/releases
    echo        then put its bin directory on PATH and re-run this script.
    exit /b 1
)
winget install --id %2 -e --accept-package-agreements --accept-source-agreements --disable-interactivity
rem winget only edits the PATH of future shells, so add the usual spots now.
set "PATH=%ProgramFiles%\CMake\bin;%LOCALAPPDATA%\Microsoft\WinGet\Links;%PATH%"
where %1 >nul 2>nul && exit /b 0
echo.
echo %1 installed but not yet on this shell's PATH.
echo Close this window, open a new terminal and run build_windows.bat again.
exit /b 1

rem ---------------------------------------------------------------------------
:fail
echo.
echo === FAILED ===
echo Scroll up for the first error. Paste it back if you want help reading it.
:done
endlocal
exit /b
