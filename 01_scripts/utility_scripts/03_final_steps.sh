#!/bin/bash
# This is basically the step (6) part of the GitHub repo when there is no ranking of 'best markers' but rather just the markers that we want. 

# Instead of combining neutral and adaptive amplicons, we will use all identified
cat 05_amplicons/completed_all_amplicons.fa > 06_output/all_amplicons.fa
