#!/usr/bin/python3

##### Code by S. Jang
##### Last updated: 2026-01-12
##### Code definition: Run EMU pipeline using parameters from a config file

import argparse
import yaml

from src.step03_emu_profile import run_emu_profile
from src.utils.mods_sjang import wait_all

def load_config(path: str) -> dict:
    with open(path, "r", encoding="utf-8") as f:
        return yaml.safe_load(f) or {}

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--config", default="config/params.yaml")
    args = p.parse_args()

    config = load_config(args.config)
    root = config["project"]["root"].rstrip("/")
    defaults = config["defaults"]
    cfg = config["emu_profile"]

    run_emu_profile(
        dir_in=f"{root}/{cfg['dir_in_rel']}",
        dir_out=f"{root}/{cfg['dir_out_rel']}",
        log_base=f"{root}/{defaults['log_base_rel']}",
        job_limit=cfg["job_limit"],
        threads=cfg["threads"],
        db=f"{root}/{cfg['db_rel']}",
        min_abundance=cfg["min_abundance"],
        extra=cfg.get("extra", ""),
        input_ext=cfg["input_ext"],
    )

    wait_all()

if __name__ == "__main__":
    main()

