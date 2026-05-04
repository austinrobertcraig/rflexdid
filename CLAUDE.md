# rflexdid — Claude Code Notes

## Running R on Windows

`Rscript` is typically not on the PATH in Windows shell environments. To run R commands, resolve the executable dynamically:

```powershell
$rscript = (Get-ChildItem "C:\Program Files\R" | Sort-Object Name -Descending | Select-Object -First 1).FullName + "\bin\Rscript.exe"
& $rscript -e "devtools::test()"
```

This picks the highest-installed R version regardless of which version is current.

## Running Tests

```powershell
$rscript = (Get-ChildItem "C:\Program Files\R" | Sort-Object Name -Descending | Select-Object -First 1).FullName + "\bin\Rscript.exe"
& $rscript -e "devtools::test()"
```

All 56 tests should pass with `FAIL 0 | WARN 0`.
