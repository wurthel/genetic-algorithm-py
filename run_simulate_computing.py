import os
import time
import random

tmp_in = "tmp_out"
tmp_out = "tmp_inp"

N = 0
with open(tmp_in, 'r') as inf:
    with open(tmp_out, 'w') as outf:
        for _ in inf.readlines():
            v = random.random()
            print(f"{v}", file=outf)
            N += 1
print(f"Computed {N} proteins")