#!/usr/bin/env python3

"""This program populates the seed file used by sim3d, specified in
$SEEDS_FILENAME (usually seeds.dat)
"""

from datetime import datetime
import random
import os
import sys

SEEDS_COUNT = os.getenv("SEEDS_COUNT", 200)

def main():
    iseed = -1 * seconds_of_year()
    generator = random.Random(iseed)

    # get filename
    seeds_fn = os.getenv("SEEDS_FILENAME")
    if not seeds_fn:
        print("SEEDS_FILENAME not specified", file=sys.stderr)
        sys.exit(1)

    nran = []

    with open(seeds_fn, "w") as seeds_outfile:
        for _ in range(SEEDS_COUNT):
            while True:
                xran = generator.random()
                cur_ran = -xran * 2_147_483_647 - 1
                if cur_ran not in nran:
                    break
            nran.append(cur_ran)
            seeds_outfile.write(str(int(cur_ran)))
            seeds_outfile.write("\n")



def seconds_of_year() -> int:
    """ Return how many seconds have passed since the
        start of the current year"""
    cur_time = datetime.now()
    start_of_year = datetime(year=cur_time.year, month=1, day=1)

    since_start_of_year = cur_time - start_of_year

    return int(since_start_of_year.total_seconds())



if __name__ == '__main__':
    main()
