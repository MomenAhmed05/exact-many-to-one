$procs = Get-CimInstance Win32_Process -Filter "name = 'python.exe'"
foreach ($p in $procs) {
    $cmd = ""
    if ($p.CommandLine) { $cmd = $p.CommandLine }
    if ($cmd -like "*w5_env*") {
        Write-Output ("killing pid " + $p.ProcessId)
        Stop-Process -Id $p.ProcessId -Force -ErrorAction SilentlyContinue
    }
}
Start-Sleep -Seconds 2
$left = Get-CimInstance Win32_Process -Filter "name = 'python.exe'" |
    Where-Object { $_.CommandLine -like "*w5_env*" }
Write-Output ("remaining w5_env python processes: " + @($left).Count)
