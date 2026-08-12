#!/bin/bash
# Tier-2 = pass count (>=20 certified) but miss survival (<15%). Usable-with-caveat; regenerable from the log.
cd /home/jovyan/tremorferometry
/home/jovyan/envs/tremorferometry/bin/python -c "
import re
print('# tier-2 survival-borderline (>=20 certified, <15% survival) -- usable-with-caveat, opt-in later')
seen=set()
for line in open('logs/fleet_broadband.log'):
    m=re.search(r'(\w+) RESULT: (\d+) certified / (\d+) densified \((\d+)% survival\) -> FLAG',line)
    if m:
        s,c,nd,v=m.group(1),int(m.group(2)),int(m.group(3)),int(m.group(4))
        if c>=20 and v<15 and s not in seen: seen.add(s); print(f'{s}\t{c} certified\t{v}% survival')
"
