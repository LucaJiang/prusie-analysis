"""Load paths from the explicitly supplied YAML configuration."""
import os
from pathlib import Path
import yaml

def load():
    path=Path(os.environ['ANALYSIS_CONFIG']).resolve()
    cfg=yaml.safe_load(path.read_text())['supplementary']
    return {key: str((path.parent/val).resolve()) for key,val in cfg.items()}
