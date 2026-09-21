$targets = @('w6_ar2', 'ar1_effdf')
foreach ($t in $targets) {
    $procs = Get-CimInstance Win32_Process -Filter "name = 'python.exe'"
    foreach ($p in $procs) {
        $cmd = ""
        if ($p.CommandLine) { $cmd = $p.CommandLine }
        if ($cmd -like "*$t*") {
            Write-Output ("killing pid " + $p.ProcessId + " (" + $t + ")")
            Stop-Process -Id $p.ProcessId -Force -ErrorAction SilentlyContinue
        }
    }
}
Start-Sleep -Seconds 2
$left = Get-CimInstance Win32_Process -Filter "name = 'python.exe'" |
    Where-Object { $_.CommandLine -like "*w6_ar2*" -or $_.CommandLine -like "*ar1_effdf*" }
Write-Output ("remaining w6_ar2/ar1_effdf processes: " + @($left).Count)
