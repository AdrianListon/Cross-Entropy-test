#!/bin/sh

src_dir=../00_src

R CMD BATCH --quiet --vanilla $src_dir/analyze_flow_cytometry.r \
    ./analyze_flow_cytometry.out


