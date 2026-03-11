#!/usr/bin/python3

##### Code by S. Jang
##### Last updated: 2026-01-11
##### Code definition: Running chopper filter

import os
import glob
from src.utils.mods_sjang import current_time, job_count_sleep

def run_chopper_filter(dir_in, dir_out, log_base, job_limit, threads, quality, minlength, maxlength, input_ext,
    label="chopper",
):
    start_time = current_time()
    os.makedirs(dir_out, exist_ok=True)

    dir_log = os.path.join(log_base, "02.chopper")
    os.makedirs(dir_log, exist_ok=True)

    input_files = sorted(glob.glob(os.path.join(dir_in, f"*.{input_ext}")))
    if not input_files:
        raise FileNotFoundError(f"No *.{input_ext} files found in: {dir_in}")

    for infile in input_files:
        base = os.path.basename(infile)
        stem, ext = os.path.splitext(base)

        outfile = os.path.join(dir_out, f"{stem}_filtered{ext}")

        log_name = f"{stem}.chopper.{start_time}.log"
        log_path = os.path.join(dir_log, log_name)

        cmd_args = [
            "chopper",
            "-t", str(int(threads)),
            "-q", str(int(quality)),
            "--minlength", str(int(minlength)),
            "--maxlength", str(int(maxlength)),
            "-i", infile,
        ]

        print(f"[RUN] {base} -> {os.path.basename(outfile)} | log: {log_name}")

        job_count_sleep(cmd_args, int(job_limit), log_path=log_path, stdout_path=outfile)

