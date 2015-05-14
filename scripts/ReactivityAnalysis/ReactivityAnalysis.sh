#!/bin/sh
./pnest bistable.seq > pnestOutput
python PijToPi.py compat.pairs > compat.pi
python PijToPi.py noncompat.pairs > noncompat.pi
