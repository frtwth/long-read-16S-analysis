##### Code by S. Jang
##### Last updated: 2026-02-09
##### Description: Run nanocomp_qc using parameters from a config file

import argparse
import yaml

from src.step00_nanocomp_qc import run_nanocomp_qc
from src.utils.mods_sjang import wait_all

def load_config(path):
    with open(path, "r") as f:
        return yaml.safe_load(f) or {}

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--config", default="config/params.yaml")
    p.add_argument("--run", choices=["raw", "trimmed", "filtered"], required=True)
    args = p.parse_args()

    config = load_config(args.config)
    root = config["project"]["root"].rstrip("/")
    defaults = config["defaults"]
    qc_runs = config["qc_runs"]
    r = qc_runs[args.run]
    cfg = config["nanocomp_qc"]
    
    run_nanocomp_qc(
        dir_in=f"{root}/{r['dir_in_rel']}",
        dir_out=f"{root}/{r['dir_out_rel']}",
        log_base=f"{root}/{defaults['log_base_rel']}",
        job_limit=cfg["job_limit"],
        threads=cfg["threads"],
        input_ext=cfg["input_ext"],
        plot = cfg["plot"], #cfg.get
        label=r.get("label"),
        run_cfg=r,
        #formats=cfg.get["formats"],
    )

    wait_all()

if __name__ == "__main__":
    main()
