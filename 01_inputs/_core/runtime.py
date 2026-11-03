"""Configuration and ordinary process I/O for input preparation."""
from pathlib import Path
import json, subprocess, time
import yaml
from .acquire import atomic_json, sha256

def config(path):
    path=Path(path).resolve()
    raw=path.read_text()
    try:data=json.loads(raw)
    except json.JSONDecodeError:data=yaml.safe_load(raw)
    # YAML 1.1 parsers treat exponents without a decimal point as strings.
    for key in ('tol','coverage','min_abs_corr','scaled_prior_variance','residual_variance',
                'null_weight','prior_tol','r_tol','check_null_threshold'):
        if key in data.get('scientific_settings',{}):
            data['scientific_settings'][key]=float(data['scientific_settings'][key])
    values=data['inputs']
    for key,value in list(values.items()):
        if isinstance(value,str) and key not in ('fasta_sha256',):
            p=Path(value).expanduser()
            values[key]=str(p if p.is_absolute() else path.parent/p)
    return data

def paths(config_path):
    cfg=config(config_path); inp=cfg['inputs']; root=Path(inp['scratch'])
    result=dict(cfg=cfg,run=root,root=root,source=Path(inp['source_pgen_prefix']),index=Path(inp['source_index']),
                samples=Path(inp['sample_dir']),fasta=Path(inp['fasta']),reference=Path(inp['reference_root']),
                inputs=root/'cases',reports=root/'reports',metadata=root/'metadata',
                plink1=Path(inp['plink1']),plink2=Path(inp['plink2']))
    for key in ('reports','metadata'): result[key].mkdir(parents=True,exist_ok=True)
    return result

def log_command(config_path,argv,log,**kwargs):
    log=Path(log);log.parent.mkdir(parents=True,exist_ok=True)
    with log.open('w') as stream:
        subprocess.run(list(map(str,argv)),stdout=stream,stderr=subprocess.STDOUT,check=True,**kwargs)

def fasta_identity(config_path,build):
    p=paths(config_path);fa=p['fasta']; st=fa.stat()
    return dict(path=str(fa),bytes=st.st_size,mtime_ns=st.st_mtime_ns,sha256=sha256(fa),build=build)
