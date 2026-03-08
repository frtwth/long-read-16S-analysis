##### Code by S. Jang
##### Last updated: 2026-01-12
##### Code definition: Running chopper filter

import os
import glob
import shlex

from src.utils.mods_sjang import current_time, job_count_sleep

def run_emu_profile(dir_in, dir_out, log_base, job_limit, threads, db, min_abundance, extra, input_ext,):
    start_time = current_time()
    os.makedirs(dir_out, exist_ok=True)

    dir_log = os.path.join(log_base, "03.emu")
    os.makedirs(dir_log, exist_ok=True)

    input_files = sorted(glob.glob(os.path.join(dir_in, f"*.{input_ext}")))
    if not input_files:
        raise FileNotFoundError(f"No *.{input_ext} files found in: {dir_in}")

    extra_args = []
    if extra is not None and str(extra).strip() != "":
        extra_args = shlex.split(str(extra))

    for infile in input_files:
        base = os.path.basename(infile)
        stem, _ = os.path.splitext(base)

        sample_outdir = os.path.join(dir_out, stem)
        os.makedirs(sample_outdir, exist_ok=True)

        log_name = f"{stem}.emu.{start_time}.log"
        log_path = os.path.join(dir_log, log_name)

        cmd = [
            "emu", "abundance",
            infile,
            "--db", db,
            "--threads", str(int(threads)),
            "--output-dir", sample_outdir,
            "--output-basename", stem,
            "--min-abundance", str(min_abundance),
            *extra_args,
        ]


        print(f"[RUN] {base} -> {stem} | log: {log_name}")
        job_count_sleep(cmd, int(job_limit), log_path=log_path)

