#!/usr/bin/python3

##### Code by S. Jang
##### Last updated: 2026-01-11
##### Code definition: Run nanocomp_qc using parameters from a config file

import argparse
import yaml

from src.step02_chopper_filter import run_chopper_filter
from src.utils.mods_sjang import wait_all

def load_config(path):
    with open(path, "r") as f:
        return yaml.safe_load(f) or {}

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--config", default="config/params.yaml")
    args = p.parse_args()

    config = load_config(args.config)
    root = config["project"]["root"].rstrip("/")
    defaults = config["defaults"]
    cfg = config["chopper_filter"]

    run_chopper_filter(
        dir_in=f"{root}/{cfg['dir_in_rel']}",
        dir_out=f"{root}/{cfg['dir_out_rel']}",
        log_base=f"{root}/{defaults['log_base_rel']}",
        job_limit=cfg["job_limit"],
        threads=cfg["threads"],
        quality=cfg["quality"],
        minlength=cfg["minlength"],
        maxlength=cfg["maxlength"], 
        input_ext=cfg["input_ext"],
    )

    wait_all()

if __name__ == "__main__":
    main()

