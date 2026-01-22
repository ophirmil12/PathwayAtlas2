import os
import sys

def boot_plot_folder():
    back_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
    sys.path.append(back_dir)