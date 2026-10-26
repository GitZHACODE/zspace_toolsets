param(
  [string]$EmsdkDir = "$env:USERPROFILE/emsdk",
  [string]$CoreDir = "$PSScriptRoot/../../../zspace_core",
  [string]$TetgenDir = "$PSScriptRoot/../cache/tetgen",
  [int]$Jobs = 6
)
$ErrorActionPreference = 'Stop'
$root = (Resolve-Path "$PSScriptRoot/../..").Path
$pythonDir = Get-ChildItem (Join-Path $EmsdkDir 'python') -Directory | Where-Object { Test-Path (Join-Path $_.FullName 'python.exe') } | Sort-Object Name -Descending | Select-Object -First 1
$python = if ($pythonDir) { Join-Path $pythonDir.FullName 'python.exe' } else { (Get-Command python -ErrorAction Stop).Source }
$emcmake = Join-Path $EmsdkDir 'upstream/emscripten/emcmake.py'
if (!(Test-Path -LiteralPath $emcmake) -or !(Test-Path -LiteralPath $python)) {
  throw "Install and activate the official Emscripten SDK first (emsdk install latest; emsdk activate latest)."
}
$cmake = (Get-Command cmake -ErrorAction SilentlyContinue).Source
if (!$cmake) { $cmake = 'C:/Program Files/Microsoft Visual Studio/2022/Professional/Common7/IDE/CommonExtensions/Microsoft/CMake/CMake/bin/cmake.exe' }
$ninja = (Get-Command ninja -ErrorAction SilentlyContinue).Source
if (!$ninja) { $ninja = 'C:/Program Files/Microsoft Visual Studio/2022/Professional/Common7/IDE/CommonExtensions/Microsoft/CMake/Ninja/ninja.exe' }
$env:EM_CONFIG = Join-Path $EmsdkDir '.emscripten'
$env:EMSDK_PYTHON = $python
$env:EM_CACHE = "$root/wasm/cache/emscripten"
$env:EMCC_TEMP_DIR = "$root/wasm/tmp"
New-Item -ItemType Directory -Path $env:EMCC_TEMP_DIR -Force | Out-Null
$env:PATH = "$(Split-Path $cmake);$(Split-Path $ninja);$EmsdkDir/upstream/emscripten;$env:PATH"
& $python $emcmake $cmake -S "$root/wasm/slicing" -B "$root/wasm/build/slicing" -G Ninja "-DCMAKE_MAKE_PROGRAM=$ninja" "-DZSPACE_CORE_DIR=$CoreDir" "-DZSPACE_TETGEN_SOURCE_DIR=$TetgenDir"
if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }
& $cmake --build "$root/wasm/build/slicing" --parallel $Jobs
exit $LASTEXITCODE
