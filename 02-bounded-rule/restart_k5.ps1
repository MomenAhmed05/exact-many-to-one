$killed = @()
# 1) the S1.5 chain shells
Get-CimInstance Win32_Process | Where-Object { $_.CommandLine -like '*W26_ENV_KS=5*' -and $_.Name -eq 'bash.exe' } | ForEach-Object {
    Write-Output ("kill chain shell " + $_.ProcessId)
    Stop-Process -Id $_.ProcessId -Force -ErrorAction SilentlyContinue
}
# 2) the S1.5 parent + effdf parent (by script name)
foreach ($t in @('w5_env', 'ar1_effdf')) {
    Get-CimInstance Win32_Process -Filter "name = 'python.exe'" | ForEach-Object {
        if ($_.CommandLine -and $_.CommandLine -like "*$t*") {
            Write-Output ("kill parent " + $_.ProcessId + " (" + $t + ")")
            Stop-Process -Id $_.ProcessId -Force -ErrorAction SilentlyContinue
        }
    }
}
# 3) orphaned heavy workers: no script name in argv, identified by CPU burn
Get-CimInstance Win32_Process -Filter "name = 'python.exe'" | ForEach-Object {
    $pr = Get-Process -Id $_.ProcessId -ErrorAction SilentlyContinue
    if ($pr -and $pr.CPU -gt 1200) {
        Write-Output ("kill heavy worker " + $_.ProcessId + " cpu=" + [int]$pr.CPU)
        Stop-Process -Id $_.ProcessId -Force -ErrorAction SilentlyContinue
    }
}
Start-Sleep -Seconds 3
$left = @(Get-CimInstance Win32_Process -Filter "name = 'python.exe'")
Write-Output ("python processes remaining: " + $left.Count)
foreach ($p in $left) {
    $pr = Get-Process -Id $p.ProcessId -ErrorAction SilentlyContinue
    $m = [regex]::Match($p.CommandLine, '[a-z0-9_]+\.py')
    $name = if ($m.Success) { $m.Value } else { '(worker)' }
    Write-Output ("  {0,7} {1,-16} cpu={2:N0}s" -f $p.ProcessId, $name, $pr.CPU)
}
