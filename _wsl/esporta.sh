#!/bin/bash
mkdir -p /mnt/i/AA_ENGINE/runs_cfd
cp -f "$1" "/mnt/i/AA_ENGINE/runs_cfd/$(basename "$1")"
ls -la "/mnt/i/AA_ENGINE/runs_cfd/$(basename "$1")"
which ffmpeg /root/zef/bin/ffmpeg
