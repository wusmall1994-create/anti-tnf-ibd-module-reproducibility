param(
  [string]$Tar = ""
)

$ErrorActionPreference = "Stop"

$Root = Split-Path -Parent (Split-Path -Parent $MyInvocation.MyCommand.Path)
Set-Location $Root

$Python = Join-Path $Root ".venv\Scripts\python.exe"

if (-not $Tar) {
  throw "Please provide the downloaded GSE282122 processed tar path, for example: powershell -ExecutionPolicy Bypass -File scripts\run_geo_10x_full_gate.ps1 -Tar <LOCAL_DATA_DIR>\GSE282122_filtered_processed_data.tar.gz"
}

if (-not (Test-Path $Python)) {
  throw "Missing project Python: $Python"
}

if (-not (Test-Path $Tar)) {
  throw "Missing GEO processed tar: $Tar"
}

Write-Host "Checking download status..."
& $Python "scripts\check_download_status.py"

$size = (Get-Item -LiteralPath $Tar).Length
$expected = 3027066520
if ($size -lt $expected) {
  throw "Tar file appears incomplete: $size bytes. Wait for download to finish, then rerun."
}

Write-Host "Auditing tar structure..."
& $Python "scripts\audit_gse282122_processed_tar.py"

Write-Host "Scoring whole-biopsy 10x modules..."
& $Python "scripts\score_10x_h5_modules_from_tar.py" --label geo_10x_full

Write-Host "Computing recovery statistics..."
& $Python "scripts\compute_recovery_synchrony.py" `
  --score-csv "results\module_scores\geo_10x_full_sample_module_scores.csv" `
  --label geo_10x_full_siteaware `
  --pair-by-site

Write-Host "Plotting recovery figure..."
& $Python "scripts\plot_recovery_synchrony.py" `
  --recovery-csv "results\recovery\geo_10x_full_siteaware_patient_lineage_module_recovery_for_tests.csv" `
  --synchrony-csv "results\recovery\geo_10x_full_siteaware_patient_synchrony.csv" `
  --label geo_10x_full_siteaware

Write-Host "Summarizing gate decision..."
& $Python "scripts\summarize_recovery_gate.py" `
  --tests "results\recovery\geo_10x_full_siteaware_recovery_synchrony_tests.csv" `
  --label geo_10x_full_siteaware

Write-Host "Done. Key outputs:"
Write-Host "results\module_scores\geo_10x_full_sample_module_scores.csv"
Write-Host "results\recovery\geo_10x_full_siteaware_recovery_synchrony_tests.csv"
Write-Host "results\recovery\geo_10x_full_siteaware_gate_decision.json"
Write-Host "results\figures\fig1_geo_10x_full_siteaware_recovery_synchrony.png"
