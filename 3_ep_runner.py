import os
from time import strftime, localtime, time
import sys
import platform
import multiprocessing as mp
import shutil
from eppy.runner.run_functions import run


ROOT_DIR = os.path.abspath(os.path.dirname(__file__))
EP_DIR = os.path.join(ROOT_DIR, 'EnergyPlus')
IDF_DIR = os.path.join(ROOT_DIR, 'idf_files_overhang')
EPW_DIR = os.path.join(ROOT_DIR, 'weather_data')
OUT_DIR = os.path.join(ROOT_DIR, 'ep_outputs')

weather_file = 'IND_Ahmedabad.426470_IWEC.epw'

# Set path to E+ idd file based on the computer's operating system
system = platform.system().lower()
if system in ['windows', 'linux', 'darwin']:
    iddfile = os.path.join(EP_DIR, 'ep8.9_{}/Energy+.idd'.format(system))

# Number of available threads for simulation
# Preset custom value N for maximum allowed number of parallel processes
# (usually to leave some processing power available for other pc activities)
N = int(8)
available_nodes = min(mp.cpu_count(), N)  # Number of available nodes.


def main():
    start = time()
    print(strftime('%H:%M:%S', localtime()),
          '- EnergyPlus simulations start time')
    sys.stdout.flush()

    epwfile = os.path.join(EPW_DIR, weather_file)

    sim_inputs = list()
    with os.scandir(IDF_DIR) as d:
        for entry in d:
            idfname = entry.path
            outfolder = os.path.join(OUT_DIR, os.path.splitext(entry.name)[0])
            sim_inputs.append([idfname, epwfile, outfolder])

    # Parallel processing
    p = mp.Pool(available_nodes)
    p.starmap_async(ep_runner, sim_inputs)
    p.close()
    p.join()

    print(strftime('%H:%M:%S', localtime()),
          '- EnergyPlus simulations finish time')
    sys.stdout.flush()
    pt('##### All EnergyPlus simulations completed in:', start)

# END OF MAIN  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


def pt(printout, pst):
    pft = time()
    process_time = pft - pst
    if process_time <= 60:
        unit = 'sec'
    elif process_time <= 3600:
        process_time = process_time / 60
        unit = 'min'
    else:
        process_time = process_time / 3600
        unit = 'hr'
    loctime = strftime('%H:%M:%S', localtime())
    print('{0} - {1} {2:.2f} {3}'.format(loctime, printout, process_time, unit))
    sys.stdout.flush()
    return

# END OF FUNCTION  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


def ep_runner(idfname, epwfile, outfolder):
    pst = time()
    # Delete output folder if exists

    # if output_dir exists then delete it
    if os.path.isdir(outfolder):
        shutil.rmtree(outfolder)

    filename = os.path.splitext(os.path.basename(idfname))[0]
    print(strftime('%H:%M:%S', localtime()),
          '- {} simulation started'.format(filename))
    sys.stdout.flush()

    run(idf=idfname, weather=epwfile, output_directory=outfolder, idd=iddfile,
        expandobjects=True, readvars=True, verbose='q', ep_version=8.9)

    pt('{} simulation completed in:'.format(filename), pst)
    return

# END OF FUNCTION  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


if __name__ == '__main__':
    main()
