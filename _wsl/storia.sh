#!/bin/bash
cp "/mnt/i/AA_ENGINE/_wsl/storia_git.py" /root/_storia_git.py
sed -i 's/\r$//' /root/_storia_git.py
/root/zef/bin/python /root/_storia_git.py
