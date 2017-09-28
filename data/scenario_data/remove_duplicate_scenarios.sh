#!/bin/bash

awk '!seen[$0]++' case14_scenarios.txt
