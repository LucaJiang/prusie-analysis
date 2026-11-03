"""Write regenerated numerical tables to build/tables without modifying inputs."""
from pathlib import Path
import os
import shutil
import subprocess
import sys
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'build/tables'
OUT.mkdir(parents=True,exist_ok=True)
for name in ('ld_class_census.tsv','thread_scaling.tsv'):
    shutil.copyfile(ROOT/'03_results'/name,OUT/name)
env=dict(os.environ,ANALYSIS_OUTPUT=str(OUT))
subprocess.run([sys.executable,str(Path(__file__).with_name('summarize_benchmark.py'))],env=env,check=True)
env['ANALYSIS_INPUT']=str(OUT)
subprocess.run([sys.executable,str(Path(__file__).with_name('derived_analyses.py'))],env=env,check=True)
subprocess.run([sys.executable,str(Path(__file__).with_name('manuscript_numbers.py')),'--output',str(OUT/'manuscript_numbers.json')],check=True)
print('Numerical tables:',OUT)
