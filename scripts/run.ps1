param(
    [Parameter(Position=0)]
    [string]$Command = "doctor",
    [Parameter(ValueFromRemainingArguments=$true)]
    [string[]]$RemainingArgs
)

$ErrorActionPreference = "Stop"
$ScriptDir = Split-Path -Parent $MyInvocation.MyCommand.Path
$SkillRoot = Split-Path -Parent $ScriptDir
$RunRoot = 'D:\ChatGPT\geneprotein-dock-runs'
$TempRoot = Join-Path $RunRoot 'tmp'
New-Item -ItemType Directory -Force -Path $TempRoot | Out-Null
$env:TEMP = $TempRoot
$env:TMP = $TempRoot
$env:TMPDIR = $TempRoot
$env:PYTHONDONTWRITEBYTECODE = '1'
$env:PYMOL_SCRIPTS = Join-Path $RunRoot 'pymol'
$env:MPLCONFIGDIR = Join-Path $RunRoot 'mplconfig'
$env:PATH = 'D:\GeneProteinDocking\env\Library\bin;D:\GeneProteinDocking\env\Scripts;' + $env:PATH

function Find-Python {
    $candidates = @(
        $env:DOCKINGSKILL_PYTHON,
        "D:\GeneProteinDocking\env\python.exe",
        "$env:CONDA_PREFIX\python.exe",
        "$env:USERPROFILE\miniforge3\envs\docking\python.exe",
        "$env:USERPROFILE\miniconda3\envs\docking\python.exe",
        "$env:LOCALAPPDATA\Programs\Python\Python312\python.exe"
    )
    foreach ($p in $candidates) {
        if ($p -and (Test-Path $p)) {
            return $p
        }
    }
    $pyCmd = Get-Command python.exe -ErrorAction SilentlyContinue
    if ($pyCmd) {
        return $pyCmd.Source
    }
    throw "No suitable Python environment found. Expected D:\GeneProteinDocking\env\python.exe"
}

$Python = Find-Python
$env:PYTHONUTF8 = "1"
$PythonRoot = Split-Path -Parent $Python
$env:PATH = (Join-Path $PythonRoot 'Library\bin') + ';' + (Join-Path $PythonRoot 'Scripts') + ';' + $env:PATH
if (Test-Path "D:\GeneProteinDocking\env\Library\bin\java.exe") {
    $env:JAVA_HOME = "D:\GeneProteinDocking\env\Library"
}

$NeedsPreflight = $Command.ToLower() -in @('dock', 'dock-single', 'batch', 'dock-batch')
$HelpRequested = ($RemainingArgs -contains '--help') -or ($RemainingArgs -contains '-h')
if ($NeedsPreflight -and -not $HelpRequested) {
    & $Python (Join-Path $ScriptDir 'doctor.py') '--local-only'
    if ($LASTEXITCODE -ne 0) {
        Write-Host 'Preflight failed: required local tools are missing or unusable. No docking job was started.'
        exit 1
    }
}

switch ($Command.ToLower()) {
    "doctor" {
        & $Python (Join-Path $ScriptDir "doctor.py") @RemainingArgs
    }
    { $_ -in "dock", "dock-single" } {
        & $Python (Join-Path $ScriptDir "dock_pipeline.py") @RemainingArgs
    }
    { $_ -in "batch", "dock-batch" } {
        & $Python (Join-Path $ScriptDir "dock_pipeline.py") "--batch" @RemainingArgs
    }
    "analyze" {
        & $Python (Join-Path $ScriptDir "analyze_results.py") @RemainingArgs
    }
    "export-pose" {
        & $Python (Join-Path $ScriptDir "export_pose.py") @RemainingArgs
    }
    "render" {
        & $Python (Join-Path $ScriptDir "render_results.py") @RemainingArgs
    }
    default {
        Write-Host "Unknown command: $Command"
        Write-Host "Supported commands: doctor, dock, batch, analyze, export-pose, render"
        exit 1
    }
}
exit $LASTEXITCODE
