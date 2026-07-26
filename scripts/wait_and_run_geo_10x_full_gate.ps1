param(
  [string]$Tar = ""
)

$ErrorActionPreference = "Stop"

$Root = Split-Path -Parent (Split-Path -Parent $MyInvocation.MyCommand.Path)
Set-Location $Root

if (-not $Tar) {
  throw "Please provide the downloaded GSE282122 processed tar path, for example: powershell -ExecutionPolicy Bypass -File scripts\wait_and_run_geo_10x_full_gate.ps1 -Tar <LOCAL_DATA_DIR>\GSE282122_filtered_processed_data.tar.gz"
}

$Expected = 3027066520
$LogDir = Join-Path $Root "results\logs"
New-Item -ItemType Directory -Force -Path $LogDir | Out-Null
$Log = Join-Path $LogDir "wait_and_run_geo_10x_full_gate.log"

function Write-Log {
  param([string]$Message)
  $stamp = Get-Date -Format "yyyy-MM-dd HH:mm:ss"
  "$stamp`t$Message" | Tee-Object -FilePath $Log -Append
}

Write-Log "Monitor started. Waiting for $Tar to reach $Expected bytes."

while ($true) {
  if (Test-Path $Tar) {
    $size = (Get-Item -LiteralPath $Tar).Length
    $pct = [math]::Round(($size / $Expected) * 100, 2)
    Write-Log "Current size: $size bytes ($pct%)."
    if ($size -ge $Expected) {
      break
    }
  } else {
    Write-Log "Tar file not found yet."
  }
  Start-Sleep -Seconds 60
}

Write-Log "Tar appears complete. Running full GEO 10x gate."
try {
  & (Join-Path $Root "scripts\run_geo_10x_full_gate.ps1") -Tar $Tar *>&1 | Tee-Object -FilePath $Log -Append
  Write-Log "Full GEO 10x gate completed."
} catch {
  Write-Log "ERROR: $($_.Exception.Message)"
  throw
}
