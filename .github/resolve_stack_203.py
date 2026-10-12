"""#203 scratch helper: resolve every stack word in a gdb dump to
object+offset via `info proc mappings`, then addr2line it."""
import re, subprocess, sys

txt = open(sys.argv[1]).read()
maps = []
for m in re.finditer(r'^\s*(0x[0-9a-f]+)\s+(0x[0-9a-f]+)\s+0x[0-9a-f]+\s+(0x[0-9a-f]+)(?:\s+[rwxsp-]{4})?\s+(/\S+)', txt, re.M):
    maps.append((int(m.group(1), 16), int(m.group(2), 16), m.group(4)))
base = {}
for lo, hi, path in maps:
    base[path] = min(base.get(path, lo), lo)
stack = re.findall(r'0x[0-9a-f]{6,16}', txt.split('### STACK', 1)[-1].split('### END', 1)[0])
seen = set()
for w in stack:
    a = int(w, 16)
    if a in seen:
        continue
    seen.add(a)
    for lo, hi, path in maps:
        if lo <= a < hi and ('libnvf' not in path):
            pie = path.endswith('.so') or '.so.' in path
            off = a - base[path] if pie else a
            obj = sys.argv[2] if path.endswith('/stripped_vc') else path
            r = subprocess.run(['addr2line', '-f', '-i', '-C', '-e', obj, hex(off)],
                               capture_output=True, text=True).stdout.split()
            print(f'{w} {path.rsplit("/", 1)[-1]}+{hex(off)} -> {" ".join(r)}')
            break
