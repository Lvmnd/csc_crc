#!/usr/bin/env pwsh
# =============================================================================
# Run the full bioinformatics pipeline (Windows/PowerShell)
# Usage:
#   .\scripts\run_pipeline.ps1
#
# Pipeline stages:
#   1. Install R packages (00_install_packages.R)
#   2. Download TCGA data (01_download_tcga.R)
#   3. Download GEO validation data (02_download_geo_validation.R)
#   4. QC discovery cohort (03_qc_discovery.R)
#   5. QC validation cohort (03_qc_validation.R)
#   6. DEG analysis (04_deg_analysis.R)
#   7. CSC marker focus (05_csc_marker_focus.R)
#   8. Signature construction (06_signature_construction.R)
#   9. Validation (07_validation.R)
#  10. Enrichment (08_enrichment.R)
#  11. Immune infiltration (09_immune_infiltration.R)
#  12. Robustness check (10_robustness_check.R)
# =============================================================================

$ErrorActionPreference = "Stop"
$startTime = Get-Date

# ── Workspace ─────────────────────────────────────────────────────────────
$ProjectRoot = Split-Path -Parent (Split-Path -Parent $PSCommandPath)
$ScriptsDir = Join-Path $ProjectRoot "scripts"
$LogsDir = Join-Path $ProjectRoot "logs"
$null = New-Item -ItemType Directory -Force -Path $LogsDir

$CondaExe = "D:\Miniforge3\Scripts\conda.exe"
$EnvName = "bioinformatics"

# ── Helper ─────────────────────────────────────────────────────────────────
function Run-Step {
    param([string]$Name, [string]$Script, [int]$TimeoutSeconds = 86400)
    $logFile = Join-Path $LogsDir "$([System.IO.Path]::GetFileNameWithoutExtension($Script)).log"
    $scriptPath = Join-Path $ScriptsDir $Script

    Write-Host "`n========================================" -ForegroundColor Cyan
    Write-Host "  [$Name] Running $Script ..." -ForegroundColor Cyan
    Write-Host "  Log → $logFile" -ForegroundColor Cyan
    Write-Host "========================================" -ForegroundColor Cyan

    $stepStart = Get-Date
    $cmd = "& '$CondaExe' run -n $EnvName Rscript '$scriptPath' 2>&1"
    # Run in separate process with timeout
    $proc = Start-Process -NoNewWindow -FilePath "powershell" -ArgumentList "-NoProfile -Command $cmd" -RedirectStandardOutput $logFile -RedirectStandardError "${logFile}.err" -PassThru

    $timeout = $TimeoutSeconds * 1000
    $completed = $proc.WaitForExit($timeout)
    if (-not $completed) {
        $proc.Kill()
        Write-Host "  ❌ [$Name] TIMEOUT after ${TimeoutSeconds}s" -ForegroundColor Red
        Write-Host "  Last 20 lines of log:" -ForegroundColor Red
        Get-Content $logFile -Tail 20
        throw "Step [$Name] timed out"
    }

    $elapsed = [math]::Round(((Get-Date) - $stepStart).TotalMinutes, 1)
    if ($proc.ExitCode -eq 0) {
        Write-Host "  ✅ [$Name] completed in ${elapsed}min" -ForegroundColor Green
    } else {
        Write-Host "  ❌ [$Name] FAILED (exit code $($proc.ExitCode)) after ${elapsed}min" -ForegroundColor Red
        Write-Host "  Last 30 lines of log:" -ForegroundColor Red
        Get-Content $logFile -Tail 30
        throw "Step [$Name] failed with exit code $($proc.ExitCode). Check $logFile"
    }
}

# ── Run Steps ──────────────────────────────────────────────────────────────
try {
    Write-Host "╔════════════════════════════════════════════════╗" -ForegroundColor Magenta
    Write-Host "║   Bioinformatics Pipeline - Full Run          ║" -ForegroundColor Magenta
    Write-Host "║   $(Get-Date -Format 'yyyy-MM-dd HH:mm')                   ║" -ForegroundColor Magenta
    Write-Host "╚════════════════════════════════════════════════╝" -ForegroundColor Magenta
    Write-Host "Project root: $ProjectRoot"
    Write-Host "Conda env:    $EnvName"
    Write-Host ""

    # ── Step 0: Install R packages ──────────────────────────────────────────
    Run-Step -Name "pkg-install" -Script "00_install_packages.R" -TimeoutSeconds 1800

    # ── Step 1-2: Downloads (parallel) ──────────────────────────────────────
    $null = New-Item -ItemType Directory -Force -Path (Join-Path $ProjectRoot "data")
    Run-Step -Name "download-tcga" -Script "01_download_tcga.R" -TimeoutSeconds 7200
    Run-Step -Name "download-geo" -Script "02_download_geo_validation.R" -TimeoutSeconds 7200

    # ── Step 3: QC (parallel) ───────────────────────────────────────────────
    Run-Step -Name "qc-discovery" -Script "03_qc_discovery.R" -TimeoutSeconds 3600
    Run-Step -Name "qc-validation" -Script "03_qc_validation.R" -TimeoutSeconds 3600

    # ── Step 4: DEG analysis ────────────────────────────────────────────────
    Run-Step -Name "deg" -Script "04_deg_analysis.R" -TimeoutSeconds 3600

    # ── Step 5: CSC marker focus ────────────────────────────────────────────
    Run-Step -Name "csc-markers" -Script "05_csc_marker_focus.R" -TimeoutSeconds 3600

    # ── Step 6: Signature construction ──────────────────────────────────────
    Run-Step -Name "signature" -Script "06_signature_construction.R" -TimeoutSeconds 7200

    # ── Step 7: Validation ──────────────────────────────────────────────────
    Run-Step -Name "validation" -Script "07_validation.R" -TimeoutSeconds 7200

    # ── Step 8: Enrichment ──────────────────────────────────────────────────
    Run-Step -Name "enrichment" -Script "08_enrichment.R" -TimeoutSeconds 3600

    # ── Step 9: Immune infiltration ─────────────────────────────────────────
    Run-Step -Name "immune" -Script "09_immune_infiltration.R" -TimeoutSeconds 3600

    # ── Step 10: Robustness check ───────────────────────────────────────────
    Run-Step -Name "robustness" -Script "10_robustness_check.R" -TimeoutSeconds 3600

    # ── Done ────────────────────────────────────────────────────────────────
    $totalMin = [math]::Round(((Get-Date) - $startTime).TotalMinutes, 1)
    Write-Host "`n╔════════════════════════════════════════════════╗" -ForegroundColor Green
    Write-Host "║   ✅ PIPELINE COMPLETE                         ║" -ForegroundColor Green
    Write-Host "║   Total time: ${totalMin} minutes                ║" -ForegroundColor Green
    Write-Host "╚════════════════════════════════════════════════╝" -ForegroundColor Green
}
catch {
    Write-Host "`n❌ PIPELINE FAILED: $_" -ForegroundColor Red
    Write-Host "See logs in: $LogsDir" -ForegroundColor Red
    exit 1
}
