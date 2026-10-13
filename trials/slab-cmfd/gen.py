import sys
from pathlib import Path

sys.path.insert(0, sys.argv[1])
import axial_slab

for width, name in [(4.0, "4"), (2.0, "2"), (1.0, "1")]:
    axial_slab.WIDTH = width
    text, _ = axial_slab.slab("16")
    text = text.replace(
        'library = "reflector.msgpack"',
        'library = "../../examples/reflector/reflector.msgpack"',
    )
    text = text.replace(
        'medium = "medium-3A.msgpack"',
        'medium = "../../examples/reflector/medium-3A.msgpack"',
    )
    text = text.replace('output = "axial-16/', f'output = "out-{name}/')
    Path(f"slab-{name}.toml").write_text(text)
    cols = text.split('rows = ["')[1].split('"')[0].split()
    print(name, len(cols))
