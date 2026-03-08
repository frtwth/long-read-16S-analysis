##### Code by S. Jang
##### Last updated: 2026-02-09
##### Code definition: Running Nancomp

import os
import glob
from src.utils.mods_sjang import current_time, job_count_sleep

def make_names(input_files, suffixes):
    names = []
    for fp in input_files:
        name = os.path.splitext(os.path.basename(fp))[0]
        for suf in suffixes:
            if name.endswith(suf):
                name = name[:-len(suf)]
                break
        names.append(name)
    return names

def run_nanocomp_qc(dir_in, dir_out, log_base, job_limit, threads, input_ext, label, plot, run_cfg=None):
    start_time = current_time()
    os.makedirs(dir_out, exist_ok=True)

    dir_log = os.path.join(log_base, "00.qc")
    os.makedirs(dir_log, exist_ok=True)

    input_files = sorted(glob.glob(os.path.join(dir_in, f"*.{input_ext}")))
    if not input_files:
        raise FileNotFoundError(f"No *.{input_ext} files found in: {dir_in}")

    suffixes = []
    if run_cfg and "name_strip_suffixes" in run_cfg:
        s = run_cfg["name_strip_suffixes"]
        suffixes = [s] if isinstance(s, str) else list(s)

    names = make_names(input_files, suffixes)
    
    log_name = f"{label}.nanocomp.{start_time}.log"
    log_path = os.path.join(dir_log, log_name)

    cmd_args = [
        "NanoComp",
        "--fastq", *input_files,
        "--names", *names,
        "--outdir", dir_out,
        "--threads", str(int(threads)),
        "--plot", str(plot),
        #"--format", str(formats),
    ]

    print(f"[RUN] {label} -> {os.path.basename(dir_out)} | log: {log_name}")
    job_count_sleep(cmd_args, int(job_limit), log_path=log_path)

