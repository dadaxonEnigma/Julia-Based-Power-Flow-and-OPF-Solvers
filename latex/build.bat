@echo off
REM ── Build thesis.pdf locally with MiKTeX (no Overleaf, no timeout) ──
setlocal

set TEX=C:\Users\User\AppData\Local\Programs\MiKTeX\miktex\bin\x64
set PATH=%TEX%;%PATH%

cd /d "%~dp0"

echo [1/4] xelatex (first pass)...
xelatex -interaction=nonstopmode thesis.tex >nul

echo [2/4] biber (bibliography)...
biber thesis >nul

echo [3/4] xelatex (second pass)...
xelatex -interaction=nonstopmode thesis.tex >nul

echo [4/4] xelatex (final pass)...
xelatex -interaction=nonstopmode thesis.tex >nul

echo.
if exist thesis.pdf (
    echo Done. thesis.pdf is ready.
    start "" thesis.pdf
) else (
    echo Build FAILED. Check thesis.log for errors.
)

endlocal
