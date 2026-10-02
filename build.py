"""Regenerate the algebra, checks, data, and images without touching the papers."""
import subprocess
import sys
from pathlib import Path

HERE=Path(__file__).resolve().parent
stages=[("derive.py","ordinary"),("derive.py","twisted"),("derive.py","tate"),("derive.py","mixed"),
        ("boundaries.py",),("verify.py",),("c4.py",),("render.py",),("physical/kawahara.py",)]
for stage in stages:
    label=Path(stage[0]).stem+("_"+stage[1] if len(stage)>1 else "")
    print("Running", " ".join(stage),flush=True)
    with (HERE/"data"/(label+".log")).open("w") as log:
        subprocess.run([sys.executable,"-u",str(HERE/stage[0]),*stage[1:]],cwd=HERE,stdout=log,check=True)
print("Open",HERE/"index.html")
