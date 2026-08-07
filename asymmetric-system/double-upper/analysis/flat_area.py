import numpy as np
from concurrent.futures import ProcessPoolExecutor
from functools import partial
import os
import flat_area_helpers as helpers
import membpy

frame_path = r"./gro_snapshots/frame_{frame_num}.gro"
outpath = "./area"
first_frame = 0
last_frame = 400
frame_skip = 1
frame_spacing = 2.5

def make_outdir(outpath, dirname=None):
    """
    Create an output directory if it does not already exist.

    Parameters:
        outpath (str): The path where the output directory should be created.
        dirname (str): The name of the output directory to create.

    Returns:
        str: The full path to the created or existing output directory.
    """
    if dirname is None:
        outdir = outpath
        if not os.path.exists(outdir):
            os.makedirs(outdir)
    else:
        outdir = os.path.join(outpath, dirname)
        if not os.path.exists(outdir):
            os.makedirs(outdir)
    return outdir

frame_list = np.arange(first_frame, last_frame + 1, frame_skip)
time_list = frame_list * frame_spacing
idx_list = np.arange(len(frame_list))

areas = np.zeros(len(frame_list))
apls = np.zeros(len(frame_list))

worker = partial(helpers.process_frame, 
                 frame_path=frame_path)

cpus_available = os.cpu_count() - 1 if os.cpu_count() > 1 else 1
max_workers = min(cpus_available, len(frame_list))
with ProcessPoolExecutor(max_workers=max_workers) as executor:
    for result in executor.map(worker, idx_list):
        idx, area, apl = result
        areas[idx] = area
        apls[idx] = apl

_ = make_outdir(outpath)
area_time_series = membpy.plot_time_series({'Area': areas}, 
                                                 time_list, 
                                                 'Undulation Area Over Time',
                                                 'Time (ns)', 'Undulation Area (nm^2)')
area_time_series.write_html(os.path.join(outpath, 'undulation_area_time_series.html'))
apl_time_series = membpy.plot_time_series({'Area per Lipid': apls},
                                          time_list,
                                          'Area per Lipid Over Time',
                                          'Time (ns)', 'Area per Lipid (nm^2)')
apl_time_series.write_html(os.path.join(outpath, 'area_per_lipid_time_series.html'))

area_boxplot = membpy.box_plot({'Area': areas},
                               'Distribution of Undulation Area Measurements',
                               xlabel=None, ylabel='Undulation Area (nm^2)')
area_boxplot.write_html(os.path.join(outpath, 'undulation_area_boxplot.html'))
apl_pdf = membpy.plot_pdf({'Area per Lipid': apls},
                          'Distribution of Area per Lipid Measurements',
                          xlabel=None, ylabel='Area per Lipid (nm^2)')
apl_pdf.write_html(os.path.join(outpath, 'area_per_lipid_pdf.html'))