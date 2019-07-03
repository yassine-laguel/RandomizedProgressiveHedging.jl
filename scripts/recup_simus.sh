#!bin/bash
scp -r luke.ciment:/bettik/PROJECTS/pr-cvar/RPH_num_exps ./logdir
find ./logdir/ -mindepth 2 -maxdepth 2 -type d '!' -exec test -e "{}/full_logs.json" ';' -print | xargs rm -rf
