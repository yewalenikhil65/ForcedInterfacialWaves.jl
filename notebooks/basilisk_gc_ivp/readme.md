   ## Key Files

   | File | Purpose |
   |------|---------|
   | `FreeSurface.dat` + `velocity_interpolated_below.dat`   | Needed for initialisation in `run_gc_ivp.c`  |
   | `run_gc_ivp.c` | Basilisk IVP solver; writes dumps + interface facets |
   | `extract.c` | Post-process: dumps → VTU |
   | `vtu_to_interface_dat.py` | Post-process: VTU → interface geometry |
   | `../gc_comoving_ivp_readable.ipynb` | Notebook: solve IVP, overlay with Basilisk, generate plots |

   ## Running the Simulation and extracting the interface data for comparison
   Following commands will work in terminal provided Basilisk and preview is installed
   ```bash
   # Compile and run serially Basilisk CFD simulation (generates dumpfile/dump-*, interface_data/interface-*.dat)
   qcc -O2 -o run_gc_ivp run_gc_ivp.c -lm
   ./run_gc_ivp

   # Post-process: dumps → VTU
   qcc -O2 -o extract extract.c -lm
   ./extract

   # Post-process (needs preview installation): VTU → interface DAT (optional, for visualization)
    /Applications/ParaView-5.13.1.app/Contents/bin/pvpython \
           vtu_to_interface_dat.py \
           vtufiles/series.pvd \
           --output-dir interface_data \
           --overwrite
