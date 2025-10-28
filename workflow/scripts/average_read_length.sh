#!/bin/bash

zcat $1 | head -400 | awk 'NR % 4 == 2 {sum += length($0); count++}END{print sum/count}'
