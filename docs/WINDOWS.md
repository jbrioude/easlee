# EASLEE on Windows — click-to-run

Everything is a double-click; no command line needed.

## 1. Install (once)

Double-click **`install_easlee.bat`**.

It will:
1. install **Julia** if it isn't already (via the Microsoft Store / `winget`);
2. download and pre-compile all the libraries the model needs;
3. do the same for the web interface.

This takes 10–20 minutes the first time (it downloads a lot).  If it installs
Julia, it will ask you to close the window and double-click the file once more
so Windows can find the new `julia` command.

If automatic Julia installation fails, install it by hand from
<https://julialang.org/downloads/> and run `install_easlee.bat` again.

## 2. (Optional but recommended) Make it start fast

Double-click **`build_sysimage.bat`**.

This compiles the model + GUI into `easlee_sysimage.dll` so the interface then
launches in a few **seconds** instead of ~a minute.  It runs once and takes
10–25 minutes.  `run_gui.bat` and `EASLEE.vbs` use the sysimage automatically
once it exists (and work fine without it, just slower to start).

## 3. Run

Either:

- **`EASLEE.vbs`** — double-click.  No terminal in your face; the server runs
  minimised on the taskbar and your browser opens at <http://localhost:8080>.
  **To stop:** click the minimised window on the taskbar and close it.
- **`run_gui.bat`** — same thing, but the black console window stays visible
  (close it to stop).  Use this if you want to see messages.

For a real app feel: right-click `EASLEE.vbs` → **Send to → Desktop (create
shortcut)**, rename it `EASLEE`, and optionally set an icon via the shortcut's
Properties → Change Icon.

## Using the interface

- **Config** tab — every option from the TOML configuration, generated from the
  model schema and grouped (Run, Meteorology, Species, Release, Physics, Output,
  Output grid).  Edit, **Save config**, then **RUN**.
- **Load TOML** — type the path to any `configs\*.toml` and load it to edit.
- **Download TOML** — save your edited configuration to a `.toml` file to re-run
  later or share.
- **Figures** tab — live column tracer, vertical cross-section, particle
  positions, and raw met fields, refreshed as the run steps.
- **Console** tab — a Julia REPL bound to the running model.

## Notes

- An NVIDIA GPU is optional; without one the model runs on the CPU
  automatically.
- Met data and emission inventories are read from the paths you set in the
  Config tab — make sure those folders exist on the machine.

## Why isn't there a single prebuilt `EASLEE.exe`?

Julia's `PackageCompiler` cannot **cross-compile**: a Linux machine builds only a
Linux image, a Windows machine only a Windows `.dll`/`.exe`.  EASLEE is developed
on Linux, so instead of shipping a (impossible-to-build-here) Windows binary we
ship `build_sysimage.bat`, which produces the Windows-native fast-start image
*on your machine* — the same payload a bundled `.exe` would carry.

**Advanced (on Windows):** for a folder that runs without Julia on `PATH`, build
one with PackageCompiler's `create_app`.  That needs the GUI wrapped as a package
with a `julia_main` entry point (not included by default); ask if you want it —
it must be built on a Windows host.

Prefer the command line?  `julia --threads=auto --project=webgui
webgui\easlee_app.jl` starts the GUI, and `julia --project=easlee
scripts\easlee.jl my_run.toml` runs a config (see [CONFIG.md](CONFIG.md)).
