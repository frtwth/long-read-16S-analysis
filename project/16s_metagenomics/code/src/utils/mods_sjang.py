#!/usr/bin/python3

##### Code by S. Jang
##### Last updated: 2026-01-08
##### Description: Custom modules by S. Jang

import time
import shlex
import subprocess

_RUNNING = []

def current_time():
    return time.strftime("%y%m%d_%H%M%S", time.localtime())

def job_count_sleep(cmd, job_limit, sleep_submit=5, sleep_poll=15, log_path=None, stdout_path=None):
    global _RUNNING

    while True:
        _RUNNING = [p for p in _RUNNING if p.poll() is None]
        if len(_RUNNING) < job_limit:
            break
        time.sleep(sleep_poll)

    log_fh = None
    if log_path:
        log_fh = open(log_path, "w", encoding="utf-8")
        cmd_str = (
            " ".join(shlex.quote(x) for x in cmd)
            if isinstance(cmd, (list, tuple))
            else str(cmd)
        )
        log_fh.write(f"[CMD] {cmd_str}\n\n")
        log_fh.flush()

    out_fh = open(stdout_path, "w", encoding="utf-8") if stdout_path else None

    p = subprocess.Popen(
        cmd,
        shell=isinstance(cmd, str),
        #stdout=log_fh,
        #stderr=subprocess.STDOUT if log_fh else None,
        stdout=out_fh if out_fh else log_fh,
        stderr=(log_fh if out_fh else (subprocess.STDOUT if log_fh else None)),
        text=True
    )

    p._log_fh = log_fh
    p._out_fh = out_fh
    _RUNNING.append(p)
    time.sleep(sleep_submit)
    return p

def wait_all():
    global _RUNNING
    for p in _RUNNING:
        p.wait()
        if getattr(p, "_out_fh", None):
            p._out_fh.close()
        if getattr(p, "_log_fh", None):
            p._log_fh.close()
    _RUNNING = []

