#!/bin/bash
echo -e '\0033\0143'
echo -e "Run point-electrode model simulation for source configuration A"
python3.7 patcha_pointelec.py

echo -e '\0033\0143'
echo -e "Run insulating disc-electrode model with insulating layer simulation for source configuration A"
python3.7 patcha_discelecins.py

echo -e '\0033\0143'
echo -e "Run metal electrode model simulation for source configuration A"
python3.7 patcha_metalelec.py

echo -e '\0033\0143'
echo -e "Run metal electrode model with insulating layer simulation for source configuration A"
python3.7 patcha_metalelecins.py

echo -e '\0033\0143'
echo -e "Run point-electrode model simulation for source configuration B"
python3.7 patchb_pointelec.py

echo -e '\0033\0143'
echo -e "Run point-electrode model simulation for source configuration B divided in deep and superficial sources"
python3.7 patchb_pointelec_DS.py
