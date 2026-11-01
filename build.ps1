# One-folder build of Scheme Molecules on Windows (PowerShell).
# Run it inside the python env that has rdkit and pyside6, from any folder that
# also contains scheme_gui.py.  If scripts are blocked, run:
#   powershell -ExecutionPolicy Bypass -File .\build.ps1
$ErrorActionPreference = "Stop"
Set-Location $PSScriptRoot

python -m pip install --quiet pyinstaller
if ($LASTEXITCODE -ne 0) { throw "could not install pyinstaller" }

foreach ($p in @("build", "dist", "SchemeMolecules.spec")) {
    if (Test-Path $p) { Remove-Item -Recurse -Force $p }
}

python -m PyInstaller --noconfirm --windowed --name SchemeMolecules scheme_gui.py
if ($LASTEXITCODE -ne 0) { throw "PyInstaller failed" }

# Smoke test. The windowed exe has no console, so we read its exit code
# (0 = OK) instead of its printed output.
Write-Host "Smoke test of the built program..."
$env:QT_QPA_PLATFORM = "offscreen"
$proc = Start-Process -FilePath ".\dist\SchemeMolecules\SchemeMolecules.exe" `
        -ArgumentList "--selftest" -Wait -PassThru
Remove-Item Env:QT_QPA_PLATFORM
if ($proc.ExitCode -eq 0) {
    Write-Host "selftest: OK"
} else {
    Write-Warning ("selftest exit code " + $proc.ExitCode + " - open dist\SchemeMolecules\SchemeMolecules.exe by hand to check it")
}

Compress-Archive -Path ".\dist\SchemeMolecules" -DestinationPath ".\dist\SchemeMolecules-windows.zip" -Force
Write-Host "Done: dist\SchemeMolecules-windows.zip   (run dist\SchemeMolecules\SchemeMolecules.exe)"
