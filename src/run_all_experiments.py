#!/usr/bin/python
# -*- coding: utf-8 -*-

import itertools
import os
import multiprocessing

def launcher(_approach, _area, _total_of_items, _reloading_depth, _relocation_cost):
    if _reloading_depth == "INF":
        _reloading_depth = int(_total_of_items)
    os.system("./dtsppl --approach %s --pickuparea ../instances/%sp.tsp --deliveryarea ../instances/%sd.tsp --n %02d --l %02d --h %02d --outputsolution ../solutions/%s/%s_%02d_%02d_%02d" % (_approach, _area, _area, _total_of_items, _reloading_depth, _relocation_cost, _approach, _area, _total_of_items, _reloading_depth, _relocation_cost))

if __name__ == "__main__":

    approach = ["ILP1", "ILP2", "BRKGA", ]
    area = ["R05", "R06", "R07", "R08", "R09", ]
    total_of_items = [6, 8, 10, 12, 15, 20, ]
    reloading_depth = [1, 2, 3, 4, 5, "INF", ]
    relocation_cost = [0, 1, 2, 5, 10, 20, ]

    os.system("make clean")
    os.system("make")

    pool = multiprocessing.Pool(max(1, multiprocessing.cpu_count() - 2))

    for _product in itertools.product(approach, total_of_items, area, reloading_depth, relocation_cost):
        _approach, _total_of_items, _area, _reloading_depth, _relocation_cost = _product
        pool.apply_async(launcher, args=(_approach, _area, _total_of_items, _reloading_depth, _relocation_cost))
        
    pool.close()
    pool.join()
